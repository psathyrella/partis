(ns ighutil.subcommands.mutation-correlation
  (:import [net.sf.samtools
            SAMRecord
            SAMFileHeader
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMSequenceRecord])
  (:require [cheshire.core :refer [generate-stream]]
            [cheshire.generate :refer [add-encoder remove-encoder encode-seq]]
            [clojure.core.reducers :as r]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [hiphip.long :as long]
            [ighutil.io :as zio]
            [ighutil.sam :as sam]
            [ighutil.sam-tags :refer [TAG-EXP-MATCH]]
            [plumbing.core :refer [safe-get map-vals]]
            [primitive-math :as p]))


(defn- result-skeleton [ref-lengths]
  (into
   {}
   (for [[name length] ref-lengths]
     (let [msize (* length length)]
       [name {:length length
              :mutated (long-array msize)
              :unmutated (long-array msize)
              :count (long-array length)}] ))))

(defn- finalize-result [result]
  (map-vals
   (fn [{:keys [mutated unmutated count length] :as m}]
     (assoc m
       :mutated (vec mutated)
       :unmutated (vec unmutated)
       :count (vec count)))
   result))

(defn- unmask-base-exp-match [^SAMRecord read ^bytes bq ref-len]
  "Generates two arrays: [counts matches]
  counts contains 1/0 indicating whether a base aligns to each
  position.

  matches contains 1/0 indicating whether there is a match at
  each position."
  (let [counts (long-array ref-len 0)
        matches (long/aclone counts)]
    (doseq [i (range (alength bq))]
      (let [ref-idx (.getReferencePositionAtReadPosition read i)
            b (long (aget bq i))]
        (when (and (not= 0 ref-idx) (>= b 0) (= 0 (mod b 100)))
          (let [ref-idx (dec (int ref-idx))]
            (long/aset counts ref-idx 1)
            (long/aset matches ref-idx (/ b 100))))))
    [counts matches]))

(defn- match-by-site-of-read [^SAMRecord read
                              {:keys [length mutated unmutated count]}]
  "Given a SAM record with the TAG-EXP-MATCH tag,
   generates a vector of
   [reference-name {site-index {matches-at-site }}]"
  (let [^bytes bq (.getAttribute read TAG-EXP-MATCH)
        [^longs read-counts ^longs read-matches] (unmask-base-exp-match read bq length)]
    (assert (not (nil? bq)))
    (long/afill! [c count i read-counts] (p/+ i c))
    (doseq [i (range (alength read-matches))]
      (let [^longs tgt (if (p/not== 0 (long/aget read-matches i))
                         unmutated
                         mutated)]
        (long/doarr [[j c] read-matches]
                    (long/ainc tgt (+ (* i length) j) c))))))

(defn- match-by-site-of-records [ref-lengths records]
  "Get the number of records that match / mismatch at each other site
   conditioning on a mutation in each position"
  (let [result (result-skeleton ref-lengths)]
    (doseq [^SAMRecord read records]
      (match-by-site-of-read
       read
       (safe-get result  (sam/reference-name read))))
    result))

(defn- reference-lengths [^SAMFileHeader header]
  "Generate a map of reference name -> length"
  (->> header
       .getSequenceDictionary
       .getSequences
       (r/map (fn [^SAMSequenceRecord r]
                [(.getSequenceName r) (.getSequenceLength r)]))
       (into {})))

(defcommand mutation-correlation
  "Set up for summarizing mutation correlation"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn io/file]
               ["-o" "--out-file" "Destination path" :required true
                :parse-fn zio/writer]]}
  (assert (not= in-file out-file))
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)
              writer ^java.io.Closeable out-file]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [long-array-cls (resolve (symbol "[J"))
          ref-lengths (-> reader .getFileHeader reference-lengths)
          m (->> reader
                 .iterator
                 iterator-seq
                 (match-by-site-of-records ref-lengths))]
      (add-encoder long-array-cls encode-seq)
      (try
        (generate-stream m writer {:pretty true})
        (finally (remove-encoder long-array-cls))))))
