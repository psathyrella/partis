(ns ighutil.subcommands.distinct-by-vdjcdr3
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory
            SAMFileHeader$SortOrder])
  (:require [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [ighutil.sam :as sam]
            [ighutil.imgt :refer [strip-allele]]
            [plumbing.core :refer [frequencies-fast distinct-by]]))

(defn- vdjcdr3 [sam-records]
  (let [^SAMRecord record (first sam-records)
        vdj (into #{} (map sam/reference-name sam-records))
        cdr3-length (.getIntegerAttribute record "XL")]
    [vdj cdr3-length]))

(defcommand distinct-by-vdjcdr3
  "Distinct records by V/D/J/cdr3"
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-o" "--out-file" "Destination path":required true
                :parse-fn io/file]
               ["--[no-]compress" "Compress output?" :default true]]}
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
          to-write (->> read-iterator
                          (partition-by sam/read-name)
                          (distinct-by vdjcdr3)
                          (apply concat))]
      (.setSortOrder header SAMFileHeader$SortOrder/unsorted)
      (with-open [writer (sam/bam-writer
                          out-file
                          reader
                          :compress compress)]
        (doseq [^SAMRecord read to-write]
          (assert (not (nil? read)))
          (.addAlignment writer read))))))
