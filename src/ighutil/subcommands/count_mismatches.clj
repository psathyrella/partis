(ns ighutil.subcommands.count-mismatches
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [plumbing.core :refer [frequencies-fast]]
            [ighutil.io :as zio]
            [ighutil.sam-tags :refer [TAG-COUNT
                                       TAG-N-MISMATCHES]]))

(defn- strip-allele [^String s]
  "Remove the allele from a string"
  (let [idx (.lastIndexOf s 42)]  ; 42 == '*'
    (if (< idx 0)
      s
      (.substring s 0 idx))))

(defn- non-primary [^SAMRecord r]
  (or (.getReadUnmappedFlag r) (.getNotPrimaryAlignmentFlag r)))

(defrecord ^{:private true} MutationKey [^String reference n-aligned n-mismatches])

(defn- count-mutations-in-record [^SAMRecord read]
  [(MutationKey.
      ^String (.getReferenceName read)
      (- (.getAlignmentEnd read) (.getAlignmentStart read))
      (.getAttribute read TAG-N-MISMATCHES))
   (or (.getAttribute read TAG-COUNT) 1)])

(defn- count-mutations-in-sam [^SAMFileReader sam-file]
  (letfn [(reduce-counts [xs]
            (let [m (java.util.HashMap.)]
              (doseq [[mut record-count] xs]
                (let [[n c] (or (.get m mut) [0 0])]
                  (.put m mut [(unchecked-inc n) (+ record-count c)])))
              (into {} m)))]
  (->> sam-file
       (.iterator)
       iterator-seq
       (remove non-primary)
       (map count-mutations-in-record)
       reduce-counts)))

(defcommand count-mismatches
  "Count mutations in reads"
  {:opts-spec [["-i" "--in-file" "Source file" :required true]
               ["-o" "--out-file" "Destination path"
                :parse-fn zio/writer :required true]]}
  (with-open [sam (SAMFileReader. (io/file in-file))]
    (.setValidationStringency
     sam
     SAMFileReader$ValidationStringency/SILENT)
    (let [rows (map (fn [[{:keys [reference n-aligned n-mismatches]} [n c]]]
                      [reference n-aligned n-mismatches n c])
                    (count-mutations-in-sam sam))]
      (with-open [^java.io.Closeable out-file out-file]
        (csv/write-csv out-file [["reference" "aligned_length" "mismatches" "count" "frequency"]])
        (csv/write-csv out-file rows)))))
