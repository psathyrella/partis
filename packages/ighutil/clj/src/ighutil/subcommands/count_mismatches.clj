(ns ighutil.subcommands.count-mismatches
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency]
           [io.github.cmccoy.sam SAMUtils])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [plumbing.core :refer [frequencies-fast safe-get]]
            [ighutil.fasta :refer [extract-references]]
            [ighutil.imgt :refer [strip-allele]]
            [ighutil.io :as zio]
            [ighutil.sam :refer [TAG-COUNT
                                 TAG-N-MISMATCHES]]))

(defn- non-primary [^SAMRecord r]
  (or (.getReadUnmappedFlag r) (.getNotPrimaryAlignmentFlag r)))

(defrecord ^{:private true} MutationKey
  [^String reference n-aligned n-mismatches])

(defn- count-mutations-in-record [^SAMRecord read ref-map]
  [(->MutationKey
    ^String (.getReferenceName read)
    (SAMUtils/countAlignedBases read)
    (SAMUtils/countMutations
     read
     ^bytes (safe-get ref-map (.getReferenceName read))))
   (or (.getAttribute read TAG-COUNT) 1)])

(defn- count-mutations-in-sam [^SAMFileReader sam-file ref-map]
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
         (map #(count-mutations-in-record % ref-map))
         reduce-counts)))

(defcommand count-mismatches
  "Count mutations in reads"
  {:opts-spec [["-i" "--in-file" "Source file" :required true]
               ["-r" "--reference-file" "Reference sequence file"
                :required true]
               ["-o" "--out-file" "Destination path"
                :parse-fn zio/writer :required true]]}

  (let [ref-map (into {} (extract-references reference-file))]
    (with-open [sam (SAMFileReader. (io/file in-file))]
      (.setValidationStringency
       sam
       SAMFileReader$ValidationStringency/SILENT)
      (let [rows (map (fn [[{:keys [reference n-aligned n-mismatches]} [n c]]]
                        [reference n-aligned n-mismatches n c])
                      (count-mutations-in-sam sam ref-map))]
        (with-open [^java.io.Closeable out-file out-file]
          (csv/write-csv out-file [["reference" "aligned_length"
                                    "mismatches" "count" "frequency"]])
          (csv/write-csv out-file rows))))))
