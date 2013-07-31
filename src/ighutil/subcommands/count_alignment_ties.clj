(ns ighutil.subcommands.count-alignment-ties
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [plumbing.core :refer [frequencies-fast]]
            [ighutil.io :as zio]
            [ighutil.sam :as sam]))

(defn- n-ties [records]
  (let [scores (mapv sam/alignment-score records)
        max-score (apply max scores)]
    (->> scores
         (filter (partial = max-score))
         count)))

(defcommand count-alignment-ties
  "Count alignment ties"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn io/file]
               ["-o" "--out-file" "Destination path":required true
                :parse-fn zio/writer]
               ["--quality" "Value to use for quality score" :default 40
                :parse-fn #(Integer/valueOf ^String %)]
               ["--[no-]sorted" "Input values are sorted." :default false]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)

    (let [tie-histogram (->> reader
                             .iterator
                             iterator-seq
                             sam/partition-by-name
                             (map n-ties)
                             frequencies-fast
                             sort)]
      (with-open [writer ^java.io.Closeable out-file]
        (csv/write-csv writer (concat [["n_ties" "frequency"]]
                                      tie-histogram))))))
