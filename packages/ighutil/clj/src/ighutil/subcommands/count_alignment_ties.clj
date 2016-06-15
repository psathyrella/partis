(ns ighutil.subcommands.count-alignment-ties
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [plumbing.core :refer [frequencies-fast distinct-fast]]
            [ighutil.io :as zio]
            [ighutil.imgt :refer [strip-allele]]
            [ighutil.sam :as sam]))

(defn- n-ties [records]
  (let [records (vec records)
        names (mapv sam/reference-name records)
        scores (mapv sam/alignment-score records)
        max-score (apply max scores)]
    [(-> records first sam/ig-locus-segment)
     (->> (map vector scores names)
          (filter (comp (partial = max-score) first))
          (map (comp strip-allele second))
          distinct-fast
          count)]))

(defcommand count-alignment-ties
  "Count alignment ties, ignoring multiple hits to different alleles of the same
   gene"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn io/file]
               ["-o" "--out-file" "Destination path" :required true
                :parse-fn zio/writer]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [tie-histogram (->> reader
                             .iterator
                             iterator-seq
                             sam/partition-by-name-segment
                             (map n-ties)
                             frequencies-fast
                             (map #(apply conj %))
                             sort)]
      (with-open [writer ^java.io.Closeable out-file]
        (csv/write-csv writer (cons ["locus" "n_ties" "frequency"]
                                    tie-histogram))))))
