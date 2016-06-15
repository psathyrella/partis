(ns ighutil.subcommands.mutations-by-site
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            AlignmentBlock])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.core.reducers :as r]
            [cliopatra.command :refer [defcommand]]
            [ighutil.fasta :refer [extract-references]]
            [ighutil.io :as zio]
            [ighutil.imgt :as imgt]
            [ighutil.sam :as sam]
            [primitive-math :as p]
            [hiphip.double :as dbl]
            [hiphip.long :as long]
            [plumbing.core :refer [map-vals frequencies-fast safe-get]]))

(defn- percent-to-proportion [^bytes b]
  (let [^doubles xs (io.github.cmccoy.Primitives/bytesToDoubles b)]
    (dbl/afill! [x xs] (p/div x 100.0))))

(defn- count-bases-in-read [^SAMRecord record ref-length &
                            {:keys [drop-uncertain]
                             :or {drop-uncertain false}}]
  (let [ref-length (int ref-length)
        ^bytes qbases (.getReadBases record)
        ^bytes bq (sam/exp-match record)
        uncertain (sam/uncertain-sites record)
        ^doubles bqd (percent-to-proportion bq)
        ^longs abases (long-array (p/* 4 ref-length) 0)
        ^doubles ematch (double-array ref-length 0.0)
        ^longs coverage (long-array ref-length 0.0)]
    (doseq [^AlignmentBlock ab (.getAlignmentBlocks record)]
      (let [qstart (-> ab .getReadStart dec int)
            rstart (-> ab .getReferenceStart dec int)
            length (-> ab .getLength int)]
        (doseq [i (range length)]
          (let [ridx (p/+ rstart (int i))
                qidx (p/+ qstart (int i))
                is-certain (not (.get uncertain qidx))
                b (char (aget qbases (p/+ (int i) qstart)))
                offset (case b
                         \A 0
                         \C 1
                         \G 2
                         \T 3
                         (throw (IllegalArgumentException. "Unknown base!")))]
            (when (or (not drop-uncertain) is-certain)
              (aset-long abases (p/+ offset (p/* 4 ridx)) 1)
              (aset-double ematch ridx (aget bqd qidx))
              (aset-long coverage ridx 1))))))
    [abases ematch coverage]))

(defn count-bases-by-position [records ref-map
                               & {:keys [drop-uncertain]
                                  :or {drop-uncertain false}}]
  (let [position-translation (map-vals (comp (partial into {}) :translation)
                                       @imgt/v-gene-meta)
        result (map-vals (fn [^bytes x] (let [l (alength x)]
                                          [x
                                           (long-array (* 4 l))
                                           (double-array l)
                                           (long-array l)])) ref-map)
        process-read (fn [^SAMRecord read]
                       (let [[^bytes rbases ra re rc] (safe-get
                                                       result
                                                       (.getReferenceName read))
                             [qa qe qc] (count-bases-in-read
                                         read
                                         (alength rbases)
                                         :drop-uncertain drop-uncertain)]
                         (long/afill! [r ra q qa] (+ r q))
                         (dbl/afill! [r re q qe] (+ r q))
                         (long/afill! [r rc q qc] (+ r q))))]
    (doall (map process-read records))
    (for [[name [^bytes bases ^longs counts ^doubles ematch ^longs cov]] result
          i (range (alength bases))]
      (let [i (int i)
            ref-base (char (aget bases i))
            [a c g t] (java.util.Arrays/copyOfRange
                       counts
                       (p/* i 4)
                       (p/* (inc i) 4))]
        {:position i
         :alignment-position (get-in position-translation [name i])
         :reference name
         :n-reads (aget cov i)
         :ref-base ref-base
         :A a :C c :G g :T t
         :exp-match (aget ematch i)}))))

(defcommand mutations-by-site
  "Count mutations in reads"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn io/file]
               ["-o" "--out-file" "Destination path"
                :parse-fn zio/writer :required true]
               ["-r" "--reference-file" "Reference file" :required true]
               ["--[no-]uncertain" "Allow uncertain positions?" :default true]]}
  (let [ref-map (into {}  (extract-references reference-file))]
    (with-open [sam (sam/bam-reader in-file)]
      (.setValidationStringency
       sam
       SAMFileReader$ValidationStringency/SILENT)
      (with-open [^java.io.Closeable out-file out-file]
        (let [header ["reference" "position"
                      "alignment_position"
                      "ref_base" "n_reads" "exp_matching"
                      "A" "C" "G" "T" "N"]
              base-freqs (count-bases-by-position
                          (-> sam .iterator iterator-seq)
                          ref-map
                          :drop-uncertain (not uncertain))
              rows (->> base-freqs
                        (remove (comp zero? :n-reads))
                        (map (juxt :reference :position :alignment-position
                                   :ref-base :n-reads :exp-match
                                   :A :C :G :T :N)))]
          (csv/write-csv out-file (cons header rows)))))))
