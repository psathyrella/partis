(ns ighutil.enumerate-mutations
  (:import [net.sf.picard.reference
            ReferenceSequenceFileFactory
            ReferenceSequenceFile]
           [net.sf.samtools SAMRecord AlignmentBlock SAMFileReader])
  (:require [clojure.java.io :as io]
            [clojure.core.memoize :refer [memo-lu]]
            [cliopatra.command :refer [defcommand]]
            [ighutil.ubtree :as ub]))

;;; TODO: Better as for?
(defn- mutations-in-block [^bytes query-bases
                           ^bytes reference-bases
                           ^AlignmentBlock block]
  (let [l (.getLength block)]
    (loop [res []
           qi (dec (.getReadStart block))
           ri (dec (.getReferenceStart block))]
      (if (< qi l)
        (let [q (aget query-bases qi)
              r (aget reference-bases ri)
              eq (= q r)]
          (recur (if (not= q r)
                   (conj res {:ref-idx ri :ref (char r) :qry (char q)})
                   res)
                 (inc qi)
                 (inc ri)))
        res))))

(defn- strip-allele [^String s]
  (let [idx (.lastIndexOf s 42)]  ; 42 == '*'
    (if (< idx 0)
      s
      (.substring s 0 idx))))

(defn- identify-mutations-in-read [^SAMRecord read reference-bases]
  (let [aligned-blocks (.getAlignmentBlocks read)
        s (.getReadBases read)]
    (into [] (mapcat
              (partial mutations-in-block s reference-bases)
              aligned-blocks))))

(defn- extract-tags
  [^SAMRecord read & tag-names]
  (into {} (map
            (fn [^String n] [(keyword n) (.getAttribute read n)])
            tag-names)))

(defn- non-primary [^SAMRecord r]
  (or (.getReadUnmappedFlag r) (.getNotPrimaryAlignmentFlag r)))

(defn- identify-mutations-in-sam [^SAMFileReader sam ref-get]
  (let [sam-records (-> sam (.iterator) iterator-seq)]
    (for [^SAMRecord r (remove non-primary sam-records)]
      (let [ref-bases (ref-get (.getReferenceName r))
            muts (identify-mutations-in-read r ref-bases)]
        (merge
         {:name (.getReadName r)
          :reference (strip-allele (.getReferenceName r))
          :mutations muts
          :n-mutations (count muts)}
         (extract-tags r "XL" "XJ" "XC"))))))

(defn- reference-per-read [^SAMFileReader sam]
  (let [sam-records (-> sam (.iterator) iterator-seq)]
    (into
     {}
     (for [^SAMRecord r (remove non-primary sam-records)]
       [(.getReadName r) (strip-allele (.getReferenceName r))]))))

(defn- ref-getter [^ReferenceSequenceFile f]
  (memo-lu
   (fn [^String contig]
     (.getBases (.getSequence f contig)))))

(defn summarize-mutation-partition [coll]
  (let [sort-mutations (partial sort-by (juxt :ref-idx :ref :qry))
        {j-gene :XJ cdr3-length :XL v-gene :reference} (first coll)]
    (loop [[r & rest] coll t ub/ubtree smap {} edges {}]
      (if r
        (let [{m' :mutations name :name} r
              mutations (sort-mutations m')]
          (if-let [hit (first (ub/lookup-subs t mutations))]
            (recur
             rest
             t
             smap
             (assoc edges name (get smap hit)))
            (recur
             rest
             (ub/insert t mutations)
             (assoc smap mutations name)
             edges)))
        {:v-gene v-gene :j-gene j-gene :cdr3-length cdr3-length :edges edges}))))

(defcommand enumerate-mutations
  "Enumerate mutations by V / J"
  {:opts-spec [["-j" "--jobs" "Number of processors" :default 1
                :parse-fn #(Integer. ^String %)]]
   :bind-args-to [reference v-sam-path]}
  (let [ref (-> reference
                io/file
                ReferenceSequenceFileFactory/getReferenceSequenceFile)
        ref-get (ref-getter ref)]
    (with-open [sam (SAMFileReader. (io/file v-sam-path))]
      (let [muts (identify-mutations-in-sam sam ref-get)]
        (->> muts
             (remove (comp zero? :n-mutations)) ;  Drop sequences without mutations from
                                        ;  germline
             (sort-by (juxt :reference :XJ :XL :n-mutations))
             (partition-by (juxt :reference :XJ :XL))
             (map summarize-mutation-partition)
             (take 2)
             doall)))))

(defn test-enumeration []
  (enumerate-mutations ["ighv.fasta" "AJ_memory_001_v.bam"]))
