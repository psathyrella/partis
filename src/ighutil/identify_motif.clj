(ns ighutil.identify-motif
  (:import [net.sf.picard.reference
            ReferenceSequenceFileFactory]
           [io.github.cmccoy.dna IUPACUtils])
  (:require [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [clojure.data.csv :as csv]
            [ighutil.io :as zio]
            [ighutil.fasta :refer [extract-references]]))

(defn- motif-occurrences [^shorts motif ^bytes ref]
  "Finds occurrences of motif in ref, allowing for ambiguity"
  (let [packed-ref (IUPACUtils/packBytes ref)
        motif-length (alength motif)
        ref-length (alength ref)]
    (when (>= ref-length motif-length)
      (for [i (range (- ref-length motif-length))
            :when (IUPACUtils/areCompatible packed-ref motif i motif-length)]
        [i (String. ref i motif-length)]))))

(defcommand identify-motif
  "Find occurrences of a motif in a FASTA file"
  {:opts-spec [["-i" "--in-file" "Source FASTA file" :required true
                :parse-fn io/file]
               ["-o" "--out-file" "Destination path" :required true
                :parse-fn io/file]
               ["-m" "--motif" "Motif (IUPAC Ambiguous)"
                :required true]]}
  (let [reference-seqs (-> in-file
                           ReferenceSequenceFileFactory/getReferenceSequenceFile
                           extract-references)
        packed-motif (IUPACUtils/packBytes (.getBytes ^String motif))]
    (with-open [out (zio/writer out-file)]
      (csv/write-csv out [["reference" "position" "motif" "sequence"]])
      (doseq [[ref-name ^bytes ref-bytes] reference-seqs]
        (println ref-name)
        (let [occ (->> ref-bytes
                       (motif-occurrences packed-motif))]
          (csv/write-csv out (for [[l s] occ] [ref-name l motif s])))))))
