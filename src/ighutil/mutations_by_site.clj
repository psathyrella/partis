(ns ighutil.mutations-by-site
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency]
           [net.sf.picard.util
            SamLocusIterator
            SamLocusIterator$LocusInfo
            SamLocusIterator$RecordAndOffset]
           [net.sf.picard.reference
            ReferenceSequenceFileFactory
            ReferenceSequenceFile
            ReferenceSequence])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [ighutil.enumerate-mutations :refer [extract-refs]]
            [ighutil.io :as zio]))

(defn count-mutations-by-position [^SAMFileReader reader ref-map]
  (let [locus-iterator (SamLocusIterator. reader)
        summarize-position
        (fn [^SamLocusIterator$LocusInfo locus]
          (let [ref-name (.getSequenceName locus)
                ref-bases (get ref-map ref-name)
                pos (dec (.getPosition locus))
                ref-base (char (aget ^bytes ref-bases pos))
                record-pos (.getRecordAndPositions locus)
                bases (map #(-> (.getReadBase ^SamLocusIterator$RecordAndOffset %) char str keyword)
                           record-pos)]
            (assoc (frequencies bases)
              :position pos
              :reference ref-name
              :n-reads (.size record-pos)
              :ref-base ref-base)))]
    (->> locus-iterator
         (.iterator)
         iterator-seq
         (map summarize-position)
         (remove (comp zero? :n-reads)))))

(defcommand mutations-by-site
  "Count mutations in reads"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn io/file]
               ["-o" "--out-file" "Destination path"
                :parse-fn zio/writer :required true]
               ["-r" "--reference-file" "Reference file" :required true]]}
  (let [ref (-> reference-file
                 io/file
                 ReferenceSequenceFileFactory/getReferenceSequenceFile)
         ref-map (into {}  (extract-refs ref))]
     (with-open [sam (SAMFileReader. in-file)]
       (.setValidationStringency
        sam
        SAMFileReader$ValidationStringency/SILENT)
       (with-open [out-file out-file]
         (csv/write-csv out-file [["reference" "position" "ref-base" "n-reads" "A" "C" "G" "T" "N"]])
         (let [base-freqs (count-mutations-by-position sam ref-map)
               rows (map (juxt :reference :position :ref-base :n-reads :A :C :G :T :N) base-freqs)]
           (csv/write-csv out-file rows))))))
