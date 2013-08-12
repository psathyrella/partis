(ns ighutil.subcommands.mutation-correlation
  (:import [net.sf.samtools
            SAMRecord
            SAMFileHeader
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMSequenceRecord])
  (:require [cheshire.core :refer [generate-stream]]
            [cheshire.generate :refer [add-encoder encode-seq]]
            [clojure.core.reducers :as r]
            [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [hiphip.long :as long]
            [ighutil.io :as zio]
            [ighutil.sam :as sam]
            [ighutil.sam-tags :refer [TAG-EXP-MATCH]]
            [plumbing.core :refer [safe-get]]))

(add-encoder (resolve (symbol "[J")) encode-seq)

(defn- unmask-base-exp-match [^SAMRecord read ^bytes bq ref-len]
  "Generates two arrays: [counts matches]
  counts contains 1/0 indicating whether a base aligns to each
  position.

  matches contains 1/0 indicating whether there is a match at
  each position."
  (let [counts (long-array ref-len)
        matches (long/aclone counts)]
    (doseq [i (range (alength bq))]
      (let [ref-idx (.getReferencePositionAtReadPosition read i)
            b (long (aget bq i))]
        (when (and (not= 0 ref-idx) (>= b 0) (= 0 (mod b 100)))
          (long/aset counts ref-idx 1)
          (long/aset matches ref-idx (/ b 100)))))
    [counts matches]))

(defn- match-by-site-of-read [^SAMRecord read ref-length]
  "Given a SAM record with the TAG-EXP-MATCH tag,
   generates a vector of
   [reference-name {site-index {matches-at-site }}]"
  (let [^bytes bq (.getAttribute read TAG-EXP-MATCH)
        [^longs counts ^longs arr] (unmask-base-exp-match read bq ref-length)
        reference-name (sam/reference-name read)
        f (fn [i] (if (not= 0 (long/aget counts i))
                    {:index i
                     :count counts
                     (safe-get {true :unmutated false :mutated}
                               (not= (long/aget arr i) 0)) arr}
                    {:index i}))]
    (assert (not (nil? bq)))
    {reference-name (->> (vec (range ref-length))
                         (r/map f)
                         (into []))}))

(defn- ^longs sum-longs [^longs xs ^longs ys]
  (long/amap [x xs y ys] (+ x y)))

(defn- match-by-site-of-records [ref-lengths records]
  (letfn [(f [^SAMRecord read]
            (match-by-site-of-read
             read
             (safe-get ref-lengths (sam/reference-name read))))
          (merge-vector-item [x y]
            (assert (= (safe-get x :index) (safe-get y :index)))
            (let [sk #(select-keys % [:count :mutated :unmutated])]
              (assoc (merge-with sum-longs (sk x) (sk y))
                :index (safe-get x :index))))
          (merge-fn
            ([x y] (mapv merge-vector-item x y)))]
    (->> records
         (map f)
         (apply merge-with merge-fn))))


(defn- reference-lengths [^SAMFileHeader header]
  "Generate a map of reference name -> length"
  (->> header
       .getSequenceDictionary
       .getSequences
       (r/map (fn [^SAMSequenceRecord r]
                [(.getSequenceName r) (.getSequenceLength r)]))
       (into {})))


(defcommand mutation-correlation
  "Set up data structures for summarizing mutation correlation"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn io/file]
               ["-o" "--out-file" "Destination path":required true
                :parse-fn zio/writer]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)
              writer ^java.io.Closeable out-file]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [ref-lengths (-> reader .getFileHeader reference-lengths)
          m (->> reader
                 .iterator
                 iterator-seq
                 (match-by-site-of-records ref-lengths))]
      (generate-stream m writer {:pretty true}))))
