(ns ighutil.subcommands.distinct-by-vdjcdr3
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [ighutil.imgt :refer [strip-allele]]
            [ighutil.io :as zio]
            [ighutil.sam :as sam]
            [plumbing.core :refer [?>> distinct-by]])
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileHeader$SortOrder]))


(defn- vdjcdr3 [sam-records]
  (let [^SAMRecord record (first sam-records)
        vdj (into {} (for [^String r (mapv sam/reference-name sam-records)]
                       [(.substring r 0 4) r]))
        cdr3-length (.getIntegerAttribute record "XL")]
    [(get vdj "IGHV") (get vdj "IGHD") (get vdj "IGHJ") cdr3-length]))

(defn- frequency-atom [f a]
  (fn [x] (let [r (f x)] (swap! a update-in [r] (fnil inc 0)) r)))

(defcommand distinct-by-vdjcdr3
  "Distinct records by V/D/J/cdr3"
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-o" "--out-file" "Destination path" :required true
                :parse-fn io/file]
               ["--[no-]all-vdj" "require a V, D, and J alignment" :default true]
               ["-c" "--count-file" "Store input counts of each V/D/J/CDR3 combination to this path"
                :parse-fn zio/writer]
               ["--[no-]compress" "Compress BAM output?" :default true]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [header (.getFileHeader reader)
          read-iterator (->> reader
                             .iterator
                             iterator-seq)
          vdj-freqs (atom {})
          to-write (->> read-iterator
                        (partition-by sam/read-name)
                        (?>> all-vdj remove (fn [x] (->> x vdjcdr3 (some nil?))))
                        (distinct-by
                         (if count-file
                           (frequency-atom vdjcdr3 vdj-freqs)
                           vdjcdr3))
                        (apply concat))]
      (.setSortOrder header SAMFileHeader$SortOrder/unsorted)
      (with-open [writer (sam/bam-writer
                          out-file
                          reader
                          :compress compress)]
        (doseq [^SAMRecord read to-write]
          (assert (not (nil? read)))
          (.addAlignment writer read)))
      (when count-file
        (with-open [^java.io.Writer f count-file]
          (->> @vdj-freqs
               seq
               sort
               (map #(apply conj %))
               (cons ["v_gene" "d_gene" "j_gene" "cdr3_length" "count"])
               (csv/write-csv count-file)))))))
