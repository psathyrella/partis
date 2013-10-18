(ns ighutil.subcommands.calculate-match-probability
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency]
           [net.sf.picard.reference FastaSequenceFile]
           [io.github.cmccoy.sam SAMUtils])
  (:require [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [ighutil.fasta :refer [extract-references]]
            [ighutil.sam :refer [primary?
                                 alignment-score
                                 partition-by-name-segment
                                 read-name
                                 mapped?
                                 bam-writer
                                 TAG-EXP-MATCH]]))

(defn- cal-equal [refs reads]
  "Calculate the expectation that each position matches reference."
  (when-let [mapped (seq (filter mapped? reads))]
    (when-let [^SAMRecord primary (first (filter primary? mapped))]
      (let [sorted (sort-by (juxt primary? alignment-score)
                            #(compare %2 %1)
                            mapped)
            max-score (-> sorted
                          first
                          alignment-score)
            max-records (vec (take-while
                              (comp (partial = max-score) alignment-score)
                              sorted))
            eq-prop (SAMUtils/calculateBaseEqualProbabilities refs max-records)]
        (.setAttribute primary TAG-EXP-MATCH eq-prop))))
  reads)

(defcommand calculate-match-probability
  "Calculate the probability that each base matches the reference"
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-r" "--reference-file" "Reference sequence file"
                :required true :parse-fn io/file]
               ["-o" "--out-file" "Destination path":required true
                :parse-fn io/file]
               ["--[no-]sorted" "Input values are sorted." :default false]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [refs (->> reference-file extract-references (into {}))
          read-iterator (->> reader
                             .iterator
                             iterator-seq)
          partitioned-reads (->> read-iterator
                                 partition-by-name-segment
                                 (map vec)
                                 (mapcat #(cal-equal refs %)))]
      (with-open [writer (bam-writer out-file reader)]
        (doseq [^SAMRecordread read partitioned-reads]
          (.addAlignment writer read))))))
