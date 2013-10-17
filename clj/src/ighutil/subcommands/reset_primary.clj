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
                                 partition-by-name-type
                                 gene-type
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
                           (.setSupplementaryAlignmentFlag true)
                           (.setNotPrimaryAlignmentFlag false))
                         (cons f (rest r))))]
    (->> reads
         (partition-by gene-type)
         (mapcat update-first))))

(defn- random-tiebreak [reads]
  (let [^SAMRecord primary (first (filter primary? reads))
        sorted (sort-by alignment-score #(compare %2 %1) reads)
        max-score (-> sorted
                      first
                      alignment-score)
        max-records (vec (take-while
                          (comp (partial = max-score) alignment-score)
                          sorted))
        ^SAMRecord selection (rand-nth max-records)]
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
    reads))

(defn- assign-primary-for-partition [select-fn reads]
  "Assign a primary read for a partition of mappings for a single read"
  (if (= 1 (count reads))
    reads
    (select-fn reads)))

(defcommand reset-primary
  "Reset primary alignment, breaking ties randomly"
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-o" "--out-file" "Destination path":required true
                :parse-fn io/file]
               ["--[no-]sorted" "Input values are sorted." :default false]
               ["--compression-level" "Compression level"
                :parse-fn #(Integer/valueOf ^String %) :default 9]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [header (.getFileHeader reader)
          read-iterator (->> reader
                             .iterator
                             iterator-seq)
          partitioned-reads (->> read-iterator
                                 partition-by-name
                                 (mapcat set-supp-and-bases-per-gene)
                                 partition-by-name-type
                                 (map vec)
                                 (mapcat #(assign-primary-for-partition
                                           random-tiebreak
                                           %))
                                 )]
      (.setSortOrder header SAMFileHeader$SortOrder/unsorted)
      (with-open [writer (.makeBAMWriter (SAMFileWriterFactory.)
                                         header
                                         true
                                         ^java.io.File out-file
                                         compression-level)]
        (doseq [^SAMRecord read partitioned-reads]
          (assert (not (nil? read)))
          (.addAlignment writer read))))))
