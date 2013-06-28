(ns ighutil.subcommands.calculate-match-probability
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory]
           [net.sf.picard.reference
            IndexedFastaSequenceFile
            ReferenceSequenceFile]
           [io.github.cmccoy.sam SAMUtils])
  (:require [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [ighutil.sam :refer [primary? alignment-score
                                       partition-by-name read-name mapped?]]))

(defn- cal-equal [^ReferenceSequenceFile refs reads]
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
      (.setAttribute primary "bq" eq-prop))))
  reads)

(defcommand calculate-match-probability
  "Calculate the probability that each base matches the reference"
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-r" "--reference-file" "Reference sequence file"
                :required true :parse-fn io/file]
               ["-o" "--out-file" "Destination path":required true
                :parse-fn io/file]
               ["--quality" "Value to use for quality score" :default 40
                :parse-fn #(Integer/valueOf ^String %)]
               ["--[no-]sorted" "Input values are sorted." :default false]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [ref (IndexedFastaSequenceFile. reference-file)
          read-iterator (->> reader
                             .iterator
                             iterator-seq)
          partitioned-reads (->> read-iterator
                                 partition-by-name
                                 (map vec)
                                 (mapcat (partial
                                          cal-equal
                                          ref)))]
      (with-open [writer (.makeSAMOrBAMWriter (SAMFileWriterFactory.)
                                              (.getFileHeader reader)
                                              sorted
                                              ^java.io.File out-file)]
        (doseq [^SAMRecordread read partitioned-reads]
          (.addAlignment writer read))))))
