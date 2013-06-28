(ns ighutil.subcommands.mutations-by-site
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
            [clojure.core.reducers :as r]
            [cliopatra.command :refer [defcommand]]
            [ighutil.fasta :refer [extract-references]]
            [ighutil.io :as zio]
            [ighutil.imgt :as imgt]
            [plumbing.core :refer [map-vals frequencies-fast]]))

(defn count-mutations-by-position [^SAMFileReader reader ref-map]
  (let [locus-iterator (SamLocusIterator. reader)
        position-translation (map-vals (comp (partial into {})
                                             :translation)
                                       imgt/v-gene-meta)
        summarize-position
        (fn [^SamLocusIterator$LocusInfo locus]
          (let [ref-name (.getSequenceName locus)
                ref-bases (get ref-map ref-name)
                pos (dec (.getPosition locus))
                ref-base (char (aget ^bytes ref-bases pos))
                record-pos (.getRecordAndPositions locus)
                extract-bq (fn [^SamLocusIterator$RecordAndOffset x]
                             (let [offset (.getOffset x)
                                   ^SAMRecord record (.getRecord x)
                                   ^bytes bq (.getAttribute record "bq")]
                               (aget bq offset)))
                bases (map #(->
                             (.getReadBase ^SamLocusIterator$RecordAndOffset %)
                             char str keyword)
                           record-pos)
                exp-match (->> record-pos
                               (r/map extract-bq)
                               (r/filter (partial <= 0)) ;; Missing
                               (r/map (partial * 0.01))
                               (r/reduce +))]
            (assoc (frequencies-fast bases)
              :position pos
              :alignment-position (get-in position-translation [ref-name pos])
              :reference ref-name
              :n-reads (.size record-pos)
              :ref-base ref-base
              :exp-match exp-match)))]
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
        ref-map (into {}  (extract-references ref))]
    (with-open [sam (SAMFileReader. ^java.io.File in-file)]
      (.setValidationStringency
       sam
       SAMFileReader$ValidationStringency/SILENT)
      (with-open [^java.io.Closeable out-file out-file]
        (csv/write-csv out-file [["reference" "position"
                                  "alignment_position"
                                  "ref_base" "n_reads" "exp_matching"
                                  "A" "C" "G" "T" "N"]])
        (let [base-freqs (count-mutations-by-position sam ref-map)
              rows (map (juxt :reference :position :alignment-position
                              :ref-base :n-reads :exp-match
                              :A :C :G :T :N) base-freqs)]
          (csv/write-csv out-file rows))))))
