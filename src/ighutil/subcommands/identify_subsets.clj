(ns ighutil.subcommands.identify-subsets
  "Enumerate all mutations from germline.
   Output is a graph containing edges between nodes with subsets of mutations."
  (:import [net.sf.picard.reference
            ReferenceSequenceFileFactory]
           [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency]
           [io.github.cmccoy.sam SAMUtils])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [ighutil.io :as zio]
            [ighutil.ubtree :as ub]
            [clojure.core.reducers :as r]
            [ighutil.subcommands.enumerate-mutations :refer [identify-mutations-in-sam]]
            [ighutil.fasta :refer [extract-references]]))

(defn summarize-mutation-partition [coll]
  "Look for mutations which are supersets of the mutations in other reads."
  (let [{:keys [j-gene cdr3-length reference]} (first coll)
        coll (vec coll)
        n-seqs (count coll)]
    (loop [[r & rest] coll t ub/ubtree smap {} edges {}]
      (if r
        (let [{m :mutations name :name} r
              mutations (vec (sort m))]
          (if-let [hits (seq (ub/lookup-subs t mutations))]
            (let [subsets (map
                           (fn [h] [(get smap h)
                                    (- (count mutations) (count h))])
                           hits)]
              (assert (> (count hits) 0))
              (recur
               rest
               t
               smap
               (assoc edges name subsets)))
            (recur
             rest
             (ub/insert t mutations)
             (assoc smap mutations name)
             edges)))
        {:v-gene reference :j-gene j-gene :cdr3-length cdr3-length :edges edges
         :n-seqs n-seqs}))))

(defcommand identify-subsets
  "Identify reads containing supersets of mutations on other reads."
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
          (csv/write-csv out-file [["v_gene" "j_gene" "n_seqs" "a" "b" "dist"]])
          (doseq [mut (partition-by :reference muts)]
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
               (for [{:keys [v-gene j-gene edges n-seqs]} summaries
                     [from kvs] edges
                     [to i] kvs]
                 [v-gene j-gene n-seqs from to i])))))))))
