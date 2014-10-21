(ns ighutil.subcommands.identify-motif
  (:import  [io.github.cmccoy.dna IUPACUtils])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.string :as string]
            [cliopatra.command :refer [defcommand]]
            [ighutil.io :as zio]
            [ighutil.fasta :refer [extract-references]]))

(defn- motif-occurrences
  "Finds occurrences of motif in ref, allowing for ambiguity
  returns a lazy-seq of (location, sequence)"
  [^bytes motif ^bytes s]
  (let [packed (IUPACUtils/packBytes s)
        motif-length (alength motif)
        s-length (alength s)]
    (when (>= s-length motif-length)
      (for [i (range (- s-length motif-length))
            :when (IUPACUtils/areCompatible packed motif i motif-length)]
        [i (String. s (int i) (int motif-length))]))))

(defcommand identify-motif
  "Find occurrences of a motif in a FASTA file"
  {:opts-spec [["-i" "--in-file" "Source FASTA file" :required true
                :parse-fn io/file]
               ["-o" "--out-file" "Destination path" :required true]
               ["-m" "--motifs" "Motif (IUPAC Ambiguous), comma separated"
                :required true]]}
  (with-open [out (zio/writer out-file)]
    (csv/write-csv out [["reference" "position" "motif" "sequence"]])
    (let [reference-seqs (extract-references in-file)]
      (doseq [^String motif (string/split motifs #",")]
        (let [packed-motif (IUPACUtils/packBytes (.getBytes motif))]
          (doseq [[ref-name ^bytes ref-bytes] reference-seqs]
            (let [occ (->> ref-bytes
                           (motif-occurrences packed-motif))]
              (csv/write-csv out (for [[l s] occ] [ref-name l motif s])))))))))
