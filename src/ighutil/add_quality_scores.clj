(ns ighutil.add-quality-scores
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory])
  (:require [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]))

(defn- fill-read-quality [^Integer score ^SAMRecord read]
  (let [read-length (.getReadLength read)
        qarr (byte-array read-length (byte score))]
    (.setBaseQualities read qarr)
    read))

(defn seq-counter
  "calls callback after every n'th entry in sequence is evaluated."
  ([sequence n callback]
     (map #(do (if (= (rem %1 n) 0) (callback %1)) %2) (iterate inc 1) sequence)))

(defcommand add-quality-scores
  "Add quality scores"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn io/file]
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
    (with-open [writer (.makeSAMOrBAMWriter (SAMFileWriterFactory.)
                                            (.getFileHeader reader)
                                            sorted
                                            ^java.io.File out-file)]
      (doseq [read (->> reader
                        .iterator
                        iterator-seq)]
        (.addAlignment writer (fill-read-quality quality read))))))
