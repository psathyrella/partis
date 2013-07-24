(ns ighutil.subcommands.kmer-matrix
  (:import [net.sf.samtools
            AlignmentBlock
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory]
           [net.sf.picard.reference FastaSequenceFile]
           [net.sf.picard.util IntervalTree]
           [io.github.cmccoy sam.SAMUtils dna.IUPACUtils])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.string :as string]
            [clojure.math.combinatorics :refer [cartesian-product]]
            [cliopatra.command :refer [defcommand]]
            [ighutil.fasta :refer [extract-references]]
            [ighutil.io :as zio]
            [ighutil.csv :refer [read-typed-csv int-of-string]]
            [plumbing.core :refer [frequencies-fast ?>> safe-get]]
            [primitive-math :as p]))

(set! *unchecked-math* true)

(defn- kmer-mutations [k ^SAMRecord read ^bytes ref & {:keys [exclude]
                                                       :or {exclude (IntervalTree.)}}]
  "Yields [ref-kmer query-kmer] tuples for all kmers in aligned blocks of
length >= k"
  (let [k (int k)
        ^bytes query (.getReadBases read)
        too-short? (fn [^AlignmentBlock b] (p/< (.getLength b) k))
        mutations-in-block (fn [^AlignmentBlock b]
                             (let [qstart (dec (.getReadStart b))
                                   rstart (dec (.getReferenceStart b))
                                   l (.getLength b)]
                               (for [i (range (p/- l (dec k)))]
                                 (let [r (p/+ rstart (int i))
                                       q (p/+ qstart (int i))]
                                   ;; Skip any kmers which overlap with sites in
                                   ;; exclude
                                   (when-not
                                       (.minOverlapper ^IntervalTree exclude
                                                       r (p/+ r (int k)))
                                     [(String. ref r k)
                                      (String. query q k)])))))]
    (->> read
         .getAlignmentBlocks
         (remove too-short?)
         (mapcat mutations-in-block)
         (filter identity))))

(defn- resolve-ambiguous-in-ref [[[^String r q] c]]
  "Distributes c among all unambiguous kmers encoded by r"
  (let [packed-ref (-> r .getBytes
                       IUPACUtils/packBytes)
        ref-card (IUPACUtils/cardinality packed-ref)]
    (if (p/== 1 ref-card)
      [[[r q] c]]
      (let [unambig (IUPACUtils/disambiguate packed-ref)
            c (p/div (float c) (float (alength unambig)))]
        (vec (for [r (seq unambig)] [[(String. ^bytes r) q] c]))))))

(defn- kmer-mutation-matrix [m k]
  "Given a map from [ref-kmer qry-kmer] -> count,
   generates a square matrix suitable for printing."
  (let [bases [\A \C \G \T]
        kmers (->> bases
                   (repeat k)
                   (apply cartesian-product)
                   (mapv string/join))
        header (cons "reference" kmers)
        rows (for [i kmers]
               (vec (cons i (for [j kmers] (get m [i j] 0)))))]
    (cons header rows)))

(defn- ^IntervalTree interval-list-of-positions [positions]
  (when (seq positions)
    (let [f (first positions)
          positions (rest positions)]
      (loop [result (IntervalTree.)
             [start end] [f f]
             positions positions]
        (if-let [[i & r] positions]
          (if (or (= end i) (= (inc end) i))
            (recur result [start i] r)
            (do
              (.put result start end 1)
              (recur result [i i] r)))
          (do
            (.put result start end 1)
            result))))))

(defn- parse-exclude-file [^java.io.Reader r]
  (->> (read-typed-csv r {:position int-of-string})
       (sort-by (juxt :reference :position))
       (partition-by :reference)
       (map (fn [xs] [(-> xs first :reference)
                      (->> xs (map :position)
                           interval-list-of-positions)]))
       (into {})))

(defcommand kmer-matrix
  ""
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-r" "--reference-file" "Reference sequence file"
                :required true :parse-fn io/file]
               ["-k" "Kmer size" :default 4 :parse-fn #(Integer/parseInt %)]
               ["-o" "--out-file" "Destination path":required true]
               ["--exclude-positions" "Positions to exclude"
                :parse-fn zio/reader]
               ["--[no-]ambiguous" "Assign ambiguous references equally to each?"
                :default false]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [exclude (when exclude-positions
                    (with-open [^java.io.Reader r exclude-positions]
                      (parse-exclude-file r)))
          refs (with-open [f (FastaSequenceFile. reference-file true)]
                 (->> f extract-references (into {})))
          mutations-for-read (fn [^SAMRecord r]
                               (kmer-mutations
                                k
                                r
                                (safe-get refs (.getReferenceName r))
                                :exclude (get exclude
                                              (.getReferenceName r)
                                              (IntervalTree.))))
          mutation-counts (->> reader
                               .iterator
                               iterator-seq
                               (mapcat mutations-for-read)
                               frequencies-fast
                               (?>> ambiguous
                                    mapcat resolve-ambiguous-in-ref)
                               (?>> ambiguous into {}))]
      (with-open [out (zio/writer out-file)]
        (csv/write-csv out (kmer-mutation-matrix mutation-counts k))))))
