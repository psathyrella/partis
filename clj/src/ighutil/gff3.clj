(ns ighutil.gff3
  (:require [clojure.string :as string]
            [clojure.java.io :as io]
            [plumbing.core :refer [keywordize-map ?>>]])
  (:import [net.sf.picard.util IntervalTree
            IntervalTree$Node]))

(defn- parse-gff3-attributes [^String attributes &
                              {:keys [keywordize?]}]
  (let [parts (string/split attributes #";")]
    (->> parts
         (map string/trim)
         (map #(string/split % #"=" 2))
         (into {})
         (?>> keywordize? keywordize-map))))

(defn- mask-missing [^String s]
  (when (not= s ".") s))

(defn parse-gff3-record [^String line]
  (let [[seqid source type start end score strand phase attr]
        (map mask-missing (string/split line #"\t"))
        start (Integer/parseInt start)
        end (Integer/parseInt end)]
    {:seqid seqid
     :source source
     :start start
     :end end
     :attributes (parse-gff3-attributes attr :keywordize? true)}))

(defn parse-gff3 [line-iter]
  (->> line-iter
       (remove #(or (string/blank? %) (.startsWith ^String % "#")))
       (map parse-gff3-record)))

(defn gff3-to-interval-map [gff-records &
                            {:keys [key] :or {key :Name}}]
  "Creates a map of seqid -> IntervalTree, with `key` as the value in the
   tree."
  (letfn [(r [m {:keys [seqid start end attributes]}]
            (let [^IntervalTree t (get m seqid (IntervalTree.))]
              (.put t start end (get attributes key))
              (assoc! m seqid t)))]
    (persistent! (reduce r (transient {}) gff-records))))

(defn overlapping-vals [^IntervalTree tree pos]
  (let [pos (int pos)
        overlappers (.overlappers tree pos pos)]
    (->> overlappers
         iterator-seq
         (mapv #(.getValue ^IntervalTree$Node %)))))
