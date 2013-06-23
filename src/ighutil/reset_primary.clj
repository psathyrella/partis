(ns ighutil.reset-primary
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory])
  (:require [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]))

(defn primary? [^SAMRecord read]
  (not (.getNotPrimaryAlignmentFlag read)))

(defn alignment-score [^SAMRecord read]
  (.getAttribute read "AS"))

(defn- random-tiebreak [reads]
  (let [^SAMRecord primary (first (filter primary? reads))
        sorted (sort-by alignment-score #(compare %2 %1))
        max-score (-> sorted
                      first
                      alignment-score)
        max-records (vec (take-while
                          (comp (partial = max-score) (alignment-score))
                          sorted))
        selection (rand-nth max-records)]
    (when-not (= selection primary)
      (do
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
  (if (= 1 (count reads))
    reads
    (select-fn reads)))

(defcommand reset-primary
  "Reset primary alignment, breaking ties randomly"
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
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
    (let [read-iterator (->> reader
                             .iterator
                             iterator-seq)
          partitioned-reads (->> read-iterator
                                 (partition-by
                                  #(.getReadName ^SAMRecord %))
                                 (map vec)
                                 (mapcat (partial
                                          assign-primary-for-partition
                                          random-tiebreak)))]
      (with-open [writer (.makeSAMOrBAMWriter (SAMFileWriterFactory.)
                                              (.getFileHeader reader)
                                              sorted
                                              ^java.io.File out-file)]
        (doseq [^SAMRecordread read partitioned-reads]
          (.addAlignment writer read))))))
