(ns ighutil.subcommands.calculate-match-probability
  "Calculate the probability of a match/mismatch at each position averaged over
  best-scoring alignments."
  (:require [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [ighutil.fasta :refer [extract-references]]
            [ighutil.sam :refer [primary?
                                 mapped?
                                 alignment-score
                                 read-name
                                 partition-by-name-segment
                                 bam-writer
                                 TAG-EXP-MATCH]])
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency]
           [net.sf.picard.reference FastaSequenceFile]
           [io.github.cmccoy.sam SAMUtils]))

(defn- cal-equal
  "Calculate the expectation that each position matches reference."
  [refs reads]
  (when-let [mapped (seq (filter mapped? reads))]
    (when-let [^SAMRecord primary (first (filter primary? mapped))]
      (let [sorted (sort-by (juxt primary? alignment-score)
                            #(compare %2 %1)
                            mapped)
            max-score (-> sorted
                          first
                          alignment-score)
            max-records (vec (take-while
                              (comp #(= max-score %) alignment-score)
                              sorted))
            eq-prop (SAMUtils/calculateBaseEqualProbabilities refs max-records)]
        (.setAttribute primary TAG-EXP-MATCH eq-prop))))
  reads)

(defn match-probability
  "Calculate the match probability for a collection of records ordered
   by name and V/D/J segment.

   Match probability (expressed as percentage) is added as
   tag 'bq', a byte array, to the primary record."
  [sam-records refs]
  (->> sam-records
       partition-by-name-segment
       (map vec)
       (mapcat #(cal-equal refs %))))

;; Command line interface to match-probability
(defcommand calculate-match-probability
  "Calculate the probability that each base matches the reference"
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-r" "--reference-file" "Reference sequence file"
                :required true :parse-fn io/file]
               ["--[no-]compress" "Compress output?" :default true]
               ["-o" "--out-file" "Destination path"
                :required true :parse-fn io/file]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [refs (->> reference-file extract-references (into {}))
          read-iterator (->> reader
                             .iterator
                             iterator-seq)]
      (with-open [writer (bam-writer out-file reader
                                     :compress compress)]
        (doseq [^SAMRecord read (match-probability read-iterator refs)]
          (.addAlignment writer read))))))
