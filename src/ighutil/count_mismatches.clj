(ns ighutil.count-mismatches
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [plumbing.core :refer [frequencies-fast]]
            [ighutil.io :as zio]))

(defn- strip-allele [^String s]
  "Remove the allele from a string"
  (let [idx (.lastIndexOf s 42)]  ; 42 == '*'
    (if (< idx 0)
      s
      (.substring s 0 idx))))

(defn- non-primary [^SAMRecord r]
  (or (.getReadUnmappedFlag r) (.getNotPrimaryAlignmentFlag r)))

(defrecord MutationSummary [^String reference n-aligned n-mismatches])

(defn- count-mutations-in-record [^SAMRecord read]
  (MutationSummary.
   ^String (.getReferenceName read)
   (- (.getAlignmentEnd read) (.getAlignmentStart read))
   (.getAttribute read "NM")))

(defn- count-mutations-in-sam [^SAMFileReader sam-file]
  (->> sam-file
       (.iterator)
       iterator-seq
       (remove non-primary)
       (map count-mutations-in-record)
       frequencies-fast))

(defcommand count-mismatches
  "Count mutations in reads"
  {:opts-spec [["-i" "--in-file" "Source file" :required true]
               ["-o" "--out-file" "Destination path"
                :parse-fn zio/writer :required true]]}
  (with-open [sam (SAMFileReader. (io/file in-file))]
    (.setValidationStringency
     sam
     SAMFileReader$ValidationStringency/SILENT)
    (let [rows (map (fn [[{:keys [reference n-aligned n-mismatches]} v]]
                      [reference n-aligned n-mismatches v])
                    (count-mutations-in-sam sam))]
      (with-open [^java.io.Closeable out-file out-file]
        (csv/write-csv out-file [["reference" "aligned_length" "mismatches" "frequency"]])
        (csv/write-csv out-file rows)))))
