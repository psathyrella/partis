(ns ighutil.subcommands.kmer-matrix
  (:import [net.sf.samtools
            AlignmentBlock
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory]
           [net.sf.picard.reference FastaSequenceFile]
           [io.github.cmccoy sam.SAMUtils dna.IUPACUtils])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.string :as string]
            [clojure.math.combinatorics :refer [cartesian-product]]
            [cliopatra.command :refer [defcommand]]
            [ighutil.fasta :refer [extract-references]]
            [ighutil.io :as zio]
            [plumbing.core :refer [frequencies-fast]]))

(defn- kmer-mutations [^Integer k ^SAMRecord read ^bytes ref]
  "Yields [ref-kmer query-kmer] tuples for all kmers in aligned blocks of
length >= k"
  (let [^bytes query (.getReadBases read)
        too-short? (fn [^AlignmentBlock b] (< (.getLength b) k))
        mutations-in-block (fn [^AlignmentBlock b]
                             (let [qstart (dec (.getReadStart b))
                                   rstart (dec (.getReferenceStart b))
                                   l (.getLength b)]
                               (for [i (range (- l (dec k)))]
                                 [(String. ref ^Integer (+ rstart i) k)
                                  (String. query ^Integer (+ qstart i) k)])))]
    (->> read
         .getAlignmentBlocks
         (remove too-short?)
         (mapcat mutations-in-block))))

(defn- resolve-ambiguous-in-ref [[[r q] c]]
  (let [packed-ref (-> ^String r .getBytes
                       IUPACUtils/packBytes)
        ref-card (IUPACUtils/cardinality packed-ref)]
    (if (= 1 ref-card)
      [[[r q] c]]
      (let [unambig (IUPACUtils/disambiguate packed-ref)
            c (/ (float c) (alength unambig))]
        (vec (for [r (seq unambig)] [[(String. ^bytes r) q] c]))))))

(defn- print-dna-matrix [m k]
  (let [bases [\A \C \G \T]
        kmers (->> bases
                   (repeat k)
                   (apply cartesian-product)
                   (mapv string/join))
        header (cons "" kmers)
        rows (for [i kmers]
               (vec (cons i (for [j kmers] (get m [i j] 0)))))]
    (cons header rows)))

(defcommand kmer-matrix
  ""
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-r" "--reference-file" "Reference sequence file"
                :required true :parse-fn io/file]
               ["-k" "Kmer size" :default 4 :parse-fn #(Integer/parseInt %)]
               ["-o" "--out-file" "Destination path":required true
                :parse-fn io/file]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [refs (with-open [f (FastaSequenceFile. reference-file true)]
                 (->> f extract-references (into {})))
          mutations-for-read (fn [^SAMRecord r]
                               (kmer-mutations
                                k
                                r
                                (get refs (.getReferenceName r))))
          mutation-counts (->> reader
                               .iterator
                               iterator-seq
                               (mapcat mutations-for-read)
                               frequencies-fast
                               seq
                               (mapcat resolve-ambiguous-in-ref)
                               (into {}))]
      (with-open [out (zio/writer out-file)]
        (csv/write-csv out (print-dna-matrix mutation-counts k))))))
