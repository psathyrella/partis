(ns ighutil.subcommands.reset-primary
  "Resets the primary record for each sequence,
  first adding sequence data to each segment (IGHV, IGHD, IGHJ),
  then finding alleles which may be present in the data,
  finally choosing randomly among best-scoring alignments."
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
                                 reference-name
                                 partition-by-name
                                 partition-by-name-segment
                                 ig-segment
                                 bam-writer]]
            [ighutil.imgt :refer [strip-allele]]
            [plumbing.core :refer [frequencies-fast]]))

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

(defn- unlikely-alleles [alignments & {:keys [min-prop] :or {min-prop 0.1}}]
  "List alleles present at less than 'min-prop' proportion of reads for a gene"
  (letfn [(filter-max-score
            [alignments]
            (let [sorted (sort-by alignment-score #(compare %2 %1) alignments)
                  max-score (-> sorted first alignment-score)
                  to-keep (vec (take-while #(= max-score (alignment-score %)) sorted))
                  n-keep (count to-keep)]
              (into
               {}
               (for [record to-keep]
                 [(reference-name record) (/ 1.0 n-keep)]))))
          (drop-group
            [freqs]
            (let [s (reduce + (map second freqs))
                  minimum (* min-prop (float s))]
              (->> freqs
                   (filter (comp #(<= % minimum) second))
                   (map first))))]
    (->> alignments
         partition-by-name-segment
         (map filter-max-score)
         (reduce #(merge-with + %1 %2))
         seq
         sort
         (partition-by (comp strip-allele first))
         (mapcat drop-group)
         (into #{}))))

(defn- max-score-tiebreak [alignments & {:keys [rand-f to-remove]
                                         :or {rand-f rand-nth
                                              to-remove #{}}}]
  "alignments should all have the same query

   Break the tie between SAM records with identical maximum alignment
   score using rand-nth.

   Drops any alignments to "
  (assert (set? to-remove))
  (if (= 1 (count alignments))
    alignments
    (let [^SAMRecord primary (first (filter primary? alignments))
          sorted (sort-by alignment-score #(compare %2 %1) alignments)
          dropped-unlikely (remove #(-> % reference-name to-remove) sorted)]
      (if (seq dropped-unlikely)
        (let [max-score (-> dropped-unlikely
                            first
                            alignment-score)
              max-records (vec (take-while
                                (comp #(= max-score %) alignment-score)
                                dropped-unlikely))
              ^SAMRecord selection (rand-f max-records)]
          (when (and primary (not= selection primary))
            (do
              ;; Swap out the primary read and the new selection
              (doto selection
                (.setReadBases (.getReadBases primary))
                (.setBaseQualities (.getBaseQualities primary))
                (.setNotPrimaryAlignmentFlag false)
                (.setSupplementaryAlignmentFlag
                 (.getSupplementaryAlignmentFlag primary)))
              (doto primary
                (.setNotPrimaryAlignmentFlag true)
                (.setReadBases SAMRecord/NULL_SEQUENCE)
                (.setBaseQualities SAMRecord/NULL_QUALS)
                (.setSupplementaryAlignmentFlag false))))
          alignments)
        ;; drop all alignments for the read
        []))))

(defn reset-primary-record [sam-records & {:keys [randomize to-remove]
                                           :or {randomize true
                                                to-remove #{}}}]
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
                 :rand-f (if randomize rand-nth first)
                 :to-remove to-remove))))

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
  (let [to-remove (with-open [reader (SAMFileReader. ^java.io.File in-file)]
                    (->> reader
                         .iterator
                         iterator-seq
                         unlikely-alleles))]
    (binding [*out* *err*]
      (println "Removing " to-remove))
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
                            out-file
                            reader
                            :compress compress)]
          (doseq [^SAMRecord read (reset-primary-record
                                   read-iterator
                                   :randomize randomize
                                   :to-remove to-remove)]
            (assert (not (nil? read)))
            (.addAlignment writer read)))))))
