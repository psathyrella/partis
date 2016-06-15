(ns ighutil.subcommands.mutation-correlation
  (:require [cheshire.core :refer [generate-stream]]
            [cheshire.generate :refer [add-encoder remove-encoder encode-seq]]
            [clojure.core.reducers :as r]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]]
            [hiphip.long :as long]
            [ighutil.fasta :refer [extract-references]]
            [ighutil.imgt :as imgt]
            [ighutil.io :as zio]
            [ighutil.sam :as sam]
            [plumbing.core :refer [safe-get map-vals]]
            [primitive-math :as p]
            [taoensso.timbre :as timbre])
  (:import [net.sf.samtools
            SAMRecord
            SAMFileHeader
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMSequenceRecord]))


(defn- result-skeleton [refs]
  (into
   {}
   (for [[name ^bytes bases] refs]
     (let [length (or (get-in @imgt/v-gene-meta
                              [name :aligned-length])
                      (alength bases))
           msize (* length length)]
       [name {:length length
              :bases (or (get-in @imgt/v-gene-meta
                                 [name :aligned])
                         (String. bases))
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

(defn- unmask-base-exp-match
  "Generates two arrays: [counts mutations]

  counts contains 1/0 indicating whether a base aligns to each
  position.

  mutations contains 1/0 indicating whether there is a mutation at
  each position."
  [^SAMRecord read ^bytes bq ref-len]
  (let [counts (long-array ref-len 0)
        muts (long/aclone counts)
        ref-name (.getReferenceName read)]
    (doseq [i (range (alength bq))]
      (let [i (int i)
            ref-idx (.getReferencePositionAtReadPosition read i)
            idx (or (get-in @imgt/v-gene-meta
                            [ref-name :translation (dec ref-idx)])
                    (dec ref-idx))
            b (long (aget bq i))]
        (when (and (not= 0 ref-idx)            ; Is aligned
                   (>= b 0) (= 0 (mod b 100))) ; Certain match/mismatch
          (let [idx (int idx)]
            (long/aset counts idx 1)
            (long/aset muts idx (- 1 (/ b 100)))))))
    [counts muts]))

(defn- conditional-mutations-of-read
  "Given a SAM record with the TAG-EXP-MATCH tag,
   generates a vector of
   [reference-name {site-index {matches-at-site }}]"
  [^SAMRecord read
   {:keys [length mutated unmutated count]}]
  (let [length (int length)
        ^bytes bq (sam/exp-match read)
        [^longs read-counts
         ^longs read-muts] (unmask-base-exp-match read bq length)]
    (assert (not (nil? bq)))
    (long/afill! [c count i read-counts] (p/+ i c))
    (doseq [i (range (alength read-muts))]
      (let [i (int i)
            ^longs tgt (if (p/== 0 (long/aget read-muts i))
                         unmutated
                         mutated)]
        (long/doarr [[j c] read-muts]
                    (long/ainc tgt (p/+ (p/* i length) j) c))))))

(defn- conditional-mutations-of-records
  "Get the number of records that match / mismatch at each other site
   conditioning on a mutation in each position"
  [refs records]
  (let [result (result-skeleton refs)]
    (doseq [^SAMRecord read records]
      (conditional-mutations-of-read
       read
       (safe-get result (sam/reference-name read))))
    result))

(defcommand mutation-correlation
  "Set up for summarizing mutation correlation

   Resulting matrices have a row for each position, with the number of mutations
   at the column j conditioned on the state in row i"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn io/file]
               ["-r" "--reference-file" "Reference file" :required true
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
          ref-map (->> reference-file
                       extract-references
                       (into {}))
          m (->> reader
                 .iterator
                 iterator-seq
                 (conditional-mutations-of-records ref-map)
                 (filter (fn [x] (let [^longs l (-> x second :count)]
                                   (-> l long/asum (> 0)))))
                 (into {}))]
      (add-encoder long-array-cls encode-seq)
      (timbre/info "Finished: writing results")
      (try
        (generate-stream m writer {:pretty true})
        (finally (remove-encoder long-array-cls)))))
  nil)
