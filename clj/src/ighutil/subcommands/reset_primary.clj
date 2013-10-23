(ns ighutil.subcommands.reset-primary
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory
            SAMFileHeader$SortOrder])
  (:require [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [ighutil.sam :refer [primary?
                                 alignment-score
                                 partition-by-name
                                 partition-by-name-segment
                                 ig-segment
                                 bam-writer]]))

(defn- set-supp-and-bases-per-gene [reads]
  "Set read bases and quals for each gene primary record"
  (let [^SAMRecord primary (first (filter primary? reads))
        ^bytes read-seq (.getReadBases primary)
        ^bytes read-qual (.getBaseQualities primary)
        update-first (fn [r]
                       (let [^SAMRecord f (first r)]
                         (doto f
                           (.setReadBases read-seq)
                           (.setBaseQualities read-qual)
                           (.setNotPrimaryAlignmentFlag false)
                           (.setSupplementaryAlignmentFlag
                            (not= (ig-segment f) \V)))
                         (cons f (rest r))))]
    (->> reads
         (partition-by ig-segment)
         (mapcat update-first))))

(defn- max-score-tiebreak [reads {:keys [rand-f] :or {rand-f rand-nth}}]
  "Break the tie between SAM records with identical maximum alignment
   score using rand-nth."
  (if (= 1 (count reads))
    reads
    (let [^SAMRecord primary (first (filter primary? reads))
          sorted (sort-by alignment-score #(compare %2 %1) reads)
          max-score (-> sorted
                        first
                        alignment-score)
          max-records (vec (take-while
                            (comp (partial = max-score) alignment-score)
                            sorted))
          ^SAMRecord selection (rand-f max-records)]
      (when (and primary (not= selection primary))
        (do
          ;; Swap out the primary read
          (doto selection
            (.setReadBases (.getReadBases primary))
            (.setBaseQualities (.getBaseQualities primary))
            (.setNotPrimaryAlignmentFlag false))
          (doto primary
            (.setNotPrimaryAlignmentFlag true)
            (.setReadBases SAMRecord/NULL_SEQUENCE)
            (.setBaseQualities SAMRecord/NULL_QUALS))))
      reads)))

(defn reset-primary-record [sam-records & {:keys [randomize]
                                           :or {randomize true}}]
  "Given a sequence of SAMRecord objects, sorted by name, assigns a
   new primary record for each partition (if :randomize is true
   [default]), copies bases and qualities to each primary segment
   record, sets the supplementary alignment flag for each non-V segment"
  (->> sam-records
       partition-by-name
       (mapcat set-supp-and-bases-per-gene)
       partition-by-name-segment
       (map vec)
       (mapcat #(max-score-tiebreak
                 %
                 :rand-f (if randomize rand-nth first)))))

(defcommand reset-primary
  "Reset primary alignment, breaking ties randomly"
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-o" "--out-file" "Destination path":required true
                :parse-fn io/file]
               ["--[no-]randomize" "Randomize primary record among top-scoring alignments."
                :default true]
               ["--[no-]compress" "Compress output?" :default true]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [header (.getFileHeader reader)
          read-iterator (->> reader
                             .iterator
                             iterator-seq)]
      (.setSortOrder header SAMFileHeader$SortOrder/unsorted)
      (with-open [writer (bam-writer
                          header
                          out-file
                          :compress compress)]
        (doseq [^SAMRecord read (reset-primary-record read-iterator
                                                      :randomize randomize)]
          (assert (not (nil? read)))
          (.addAlignment writer read))))))
