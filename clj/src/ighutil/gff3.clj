(ns ighutil.gff3
  (:require [clojure.string :as string]
            [clojure.java.io :as io]
            [ighutil.io :as zio]
            [plumbing.core :refer [keywordize-map ?>>]])
  (:import [net.sf.picard.util Interval IntervalTreeMap]))

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

(defn slurp-gff3 [file-name]
  (-> file-name
      zio/reader
      file-seq
      parse-gff3
      (into [])))

(defn ^IntervalTreeMap gff3-to-interval-map [gff-records &
                                             {:keys [key]
                                              :or {key identity}}]
  "Creates an IntervalTreeMap, with reference intervals as keys,
   and gff-records as values"
  (let [result (IntervalTreeMap.)]
    (doseq [{:keys [seqid start end attributes] :as record} gff-records]
      (.put result (Interval. seqid start end) (key record)))
    result))

(defn overlapping [^IntervalTreeMap tree ^String seqid pos]
  "Get all intervals overlapping position `pos`"
  (let [pos (int pos)]
    (vec (.getOverlapping tree (Interval. seqid pos pos)))))

(defn all-feature-names [^IntervalTreeMap tree &
                         {:keys [key]
                          :or {key (comp :Name :attributes)}}]
  "Get all feature names from an interval tree map"
  (vec
   (for [^java.util.Map$Entry item tree]
     [(.getSequence ^Interval (.getKey item))
      (-> item .getValue key)])))
