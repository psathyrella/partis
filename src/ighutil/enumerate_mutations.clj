(ns ighutil.enumerate-mutations
  (:import [net.sf.picard.reference
            ReferenceSequenceFileFactory
            ReferenceSequenceFile
            ReferenceSequence]
           [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency]
           [io.github.cmccoy.sam SAMUtils])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [ighutil.io :as zio]
            [ighutil.ubtree :as ub]))

(def ^{:private true} n-records (atom 0))

(defn- report! []
  (let [n (swap! n-records inc)]
    (when (= 0 (mod n 50000))
      (.print System/err (format "Read record %10d\r" n)))))

(defn- strip-allele [^String s]
  (let [idx (.lastIndexOf s 42)]  ; 42 == '*'
    (if (< idx 0)
      s
      (.substring s 0 idx))))

(defn- non-primary [^SAMRecord r]
  (or (.getReadUnmappedFlag r) (.getNotPrimaryAlignmentFlag r)))

(defn- identify-mutations-on-ref [records ^bytes ref-bases]
  (for [^SAMRecord r (remove non-primary records)]
    (let [muts (SAMUtils/enumerateMutations r ref-bases)
          tag #(.getAttribute r ^String %)]
      (report!)
      {:name (.getReadName r)
       :reference (strip-allele (.getReferenceName r))
       :mutations muts
       :n-mutations (count muts)
       :cdr3-length (tag "XL")
       :j-gene (tag "XJ")
       :count (tag "XC")})))

(defn- identify-mutations-in-sam [^SAMFileReader reader ref-map]
  (when-not (.hasIndex reader)
    (throw (IllegalArgumentException. (format "SAM file lacks index!"))))
  (for [v-group (->> ref-map
                     (sort-by first)
                     (partition-by (comp strip-allele first)))]
    (apply concat
     (for [[^String ref-name ^bytes ref-bases] v-group]
       (with-open [reads (.query reader ref-name 0 0 false)]
         (let [mutations (-> reads
                             iterator-seq
                             (identify-mutations-on-ref ref-bases))]
           (doall mutations)))))))

(defn- extract-refs [^ReferenceSequenceFile f]
  (.reset f)

  (loop [sequences (transient [])]
    (if-let [^ReferenceSequence s (.nextSequence f)]
      (recur (conj! sequences [(.getName s) (.getBases s)]))
      (persistent! sequences))))

(defn summarize-mutation-partition [coll]
  (let [{:keys [j-gene cdr3-length reference]} (first coll)
        coll (vec coll)
        n-seqs (count coll)]
    (println reference j-gene cdr3-length n-seqs)
    (loop [[r & rest] coll t ub/ubtree smap {} edges {}]
      (if r
        (let [{m :mutations name :name} r
              mutations (vec (sort m))]
          (if-let [hit (first (ub/lookup-subs t mutations))]
            (recur
             rest
             t
             smap
             (assoc edges name [(get smap hit) (- (count mutations) (count hit))]))
            (recur
             rest
             (ub/insert t mutations)
             (assoc smap mutations name)
             edges)))
        {:v-gene reference :j-gene j-gene :cdr3-length cdr3-length :edges edges
         :n-seqs n-seqs}))))

(defcommand enumerate-mutations
  "Enumerate mutations by V / J"
  {:opts-spec [["-j" "--jobs" "Number of processors" :default 1
                :parse-fn #(Integer. ^String %)]
               ["-o" "--out-file" "Destination path"
                :parse-fn zio/writer :required true]]
   :bind-args-to [reference v-sam-path]}
  (let [ref (-> reference
                io/file
                ReferenceSequenceFileFactory/getReferenceSequenceFile)
        ref-map (extract-refs ref)]
    (with-open [sam (SAMFileReader. (io/file v-sam-path))]
      (.setValidationStringency
       sam
       SAMFileReader$ValidationStringency/SILENT)
      (let [muts (identify-mutations-in-sam sam ref-map)]
        (with-open [out-file out-file]
          (csv/write-csv out-file [["v_gene" "j_gene" "n_seqs" "a" "b" "i"]])
          (doseq [mut muts]
            (let [summaries (->> mut
                                 (remove (comp zero? :n-mutations))
                                 (sort-by (juxt :reference
                                                :j-gene
                                                :cdr3-lenth
                                                :n-mutations))
                                 (partition-by (juxt :reference
                                                     :j-gene
                                                     :cdr3-length))
                                 (map summarize-mutation-partition))]
              (csv/write-csv
               out-file
               (for [{:keys [reference j-gene edges n-seqs]} summaries
                     [from [to i]] edges]
                 [reference j-gene n-seqs from to i])))))))))

(defn test-enumerate []
  (enumerate-mutations ["-o" "test.csv" "ighv.fasta" "AJ_memory_001_v.bam"]))
