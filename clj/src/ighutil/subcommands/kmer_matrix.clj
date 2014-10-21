(ns ighutil.subcommands.kmer-matrix
  (:require [cliopatra.command :refer [defcommand]]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [clojure.math.combinatorics :refer [cartesian-product]]
            [clojure.string :as string]
            [ighutil.csv :refer [read-typed-csv]]
            [ighutil.fasta :refer [extract-references]]
            [ighutil.io :as zio]
            [ighutil.csv :refer [read-typed-csv int-of-string]]
            [ighutil.sam :as sam]
            [plumbing.core :refer [frequencies-fast ?>> safe-get]]
            [primitive-math :as p])
  (:import [net.sf.samtools
            AlignmentBlock
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency]
           [net.sf.picard.util IntervalTree]
           [io.github.cmccoy sam.SAMUtils dna.IUPACUtils]))

(set! *unchecked-math* true)

(defn kmer-mutations
  "Yields [ref-kmer query-kmer] tuples for all kmers in aligned blocks
of length >= k"
  [^Integer k
   ^SAMRecord read
   ^bytes ref
   & {:keys [exclude drop-uncertain?  frame]
      :or {exclude (IntervalTree.) drop-uncertain false}}]
  (let [k (int k)
        ^bytes query (.getReadBases read)
        uncertain (sam/uncertain-sites read)
        certain? (fn [start]
                   (= 0 (.. uncertain
                            (get (int start) (p/+ (int start) k))
                            (cardinality))))
        too-short? (fn [^AlignmentBlock b] (p/< (.getLength b) k))
        mutations-in-block (fn [^AlignmentBlock b]
                             (let [qstart (dec (.getReadStart b))
                                   rstart (dec (.getReferenceStart b))
                                   l (.getLength b)]
                               (for [i (range (p/- l (dec k)))]
                                 (let [r (p/+ rstart (int i))
                                       q (p/+ qstart (int i))]
                                   ;; Skip any k-mers which overlap with sites in
                                   ;; exclude
                                   (when-not
                                       (or
                                        ;; Frame specified and not in frame
                                        (and (not (nil? frame)) (not= frame (mod r k)))
                                        ;; Uncertain about ref base
                                        (and drop-uncertain? (not (certain? q)))
                                        ;; In list of positions to exclude
                                        (.minOverlapper ^IntervalTree exclude
                                                        r (p/+ r (int k))))
                                     [(String. ref r k)
                                      (String. query q k)])))))]
    (->> read
         .getAlignmentBlocks
         (remove too-short?)
         (mapcat mutations-in-block)
         (filter identity))))

(defn- resolve-ambiguous-in-ref
  "Distributes c among all unambiguous kmers encoded by r"
  [[[^String r q] c]]
  (let [packed-ref (-> r .getBytes
                       IUPACUtils/packBytes)
        ref-card (IUPACUtils/cardinality packed-ref)]
    (if (p/== 1 ref-card)
      [[[r q] c]]
      (let [unambig (IUPACUtils/disambiguate packed-ref)
            c (p/div (float c) (float (alength unambig)))]
        (vec (for [r (seq unambig)] [[(String. ^bytes r) q] c]))))))

(defn- kmer-mutation-matrix
  "Given a map from [ref-kmer qry-kmer] -> count,
   generates a square matrix suitable for printing."
  [m k]
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
      (if (seq positions)
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
              result)))
        (doto (IntervalTree.)
          (.put f f 1))))))

(defn- parse-exclude-file [^java.io.Reader r]
  (->> (read-typed-csv r {:position int-of-string})
       (sort-by (juxt :reference :position))
       (partition-by :reference)
       (map (fn [xs] [(-> xs first :reference)
                      (->> xs (map :position)
                           interval-list-of-positions)]))
       (into {})))

(defn- atoi [^String s]
  (Integer/parseInt s))

(defcommand kmer-matrix
  "k-mer mutation matrix"
  {:opts-spec [["-i" "--in-file" "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-r" "--reference-file" "Reference sequence file"
                :required true :parse-fn io/file]
               ["-k" "Kmer size" :default 4 :parse-fn atoi]
               ["-o" "--out-file" "Destination path":required true]
               ["--exclude-positions" "Positions to exclude"
                :parse-fn zio/reader]
               ["-f" "--frame" "Frame to use" :parse-fn atoi]
               ["--[no-]uncertain" "Allow uncertain positions?" :default true]
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
          refs (->> reference-file extract-references (into {}))
          mutations-for-read (fn [^SAMRecord r]
                               (kmer-mutations
                                k
                                r
                                (safe-get refs (.getReferenceName r))
                                :exclude (get exclude
                                              (.getReferenceName r)
                                              (IntervalTree.))
                                :drop-uncertain? (not uncertain)
                                :frame frame))
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
