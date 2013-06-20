(ns ighutil.enumerate-mutations
  "Enumerate all mutations from germline.
   Output is a graph containing edges between nodes with subsets of mutations."
  (:import [net.sf.picard.reference
            ReferenceSequenceFileFactory
            ReferenceSequenceFile
            ReferenceSequence]
           [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency]
           [io.github.cmccoy.sam SAMUtils Mutation])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [ighutil.io :as zio]
            [ighutil.ubtree :as ub]
            [ighutil.fasta :refer [extract-references]]
            [clojure.core.reducers :as r]))

(def ^{:private true} n-records (atom 0))

(defn- report! []
  (let [n (swap! n-records inc)]
    (when (= 0 (mod n 50000))
      (.print System/err (format "Read record %10d\r" n)))))

(defn- strip-allele [^String s]
  "Remove the allele from a string"
  (let [idx (.lastIndexOf s 42)]  ; 42 == '*'
    (if (< idx 0)
      s
      (.substring s 0 idx))))

(defn- non-primary [^SAMRecord r]
  (or (.getReadUnmappedFlag r) (.getNotPrimaryAlignmentFlag r)))

(defn- identify-mutations-for-ref [records ^bytes ref-bases ref-name]
  "Identify mutations for reads mapped to a single reference."
  (->>
   records
   (remove non-primary)
   vec
   (r/map (fn [^SAMRecord r]
            (let [muts (SAMUtils/enumerateMutations r ref-bases)
                  tag #(.getAttribute r ^String %)]
              {:name (.getReadName r)
               :reference ref-name
               :mutations muts
               :n-mutations (count muts)
               :cdr3-length (tag "XL")
               :j-gene (tag "XJ")
               :count (tag "XC")})))
   (into [])))

(defn identify-mutations-in-sam [^SAMFileReader reader ref-map]
  "Identify mutations from reference sequence in a SAM file."
  (when-not (.hasIndex reader)
    (throw (IllegalArgumentException. (format "SAM file lacks index!"))))
  (apply
   concat
   (for [v-group (->> ref-map
                      (sort-by first)
                      (partition-by (comp strip-allele first)))]
     (apply
      concat
      (for [[^String ref-name ^bytes ref-bases] v-group]
        (with-open [reads (.query reader ref-name 0 0 false)]
          (let [mutations (-> reads
                              iterator-seq
                              (identify-mutations-for-ref
                               ref-bases
                               (strip-allele ref-name)))]
            ;; Have to consume the seq here:
            ;; Only one SAM iterator may be open at a time.
            (doall mutations))))))))

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
        ref-map (extract-references ref)]
    (with-open [sam (SAMFileReader. (io/file v-sam-path))]
      (.setValidationStringency
       sam
       SAMFileReader$ValidationStringency/SILENT)
      (let [muts (identify-mutations-in-sam sam ref-map)]
        (with-open [^java.io.Closeable out-file out-file]
          (csv/write-csv out-file [["sequence" "reference" "location" "type" "wt" "mut"]])
          (let [rows (for [{:keys [mutations name reference] :as read} muts
                            ^Mutation m mutations]
                        [name reference (.. m (getType) (getCode))
                         (.getPosition m)
                         (.getWildType m)
                         (.getMutant m)])]
            (csv/write-csv out-file rows)))))))
