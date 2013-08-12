(ns ighutil.sam
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory]
           [java.util BitSet])
  (:require [clojure.java.io :as io]
            [ighutil.sam-tags :refer [TAG-EXP-MATCH]]))

(defn bam-reader [path]
  (-> path io/file SAMFileReader.))

(defn bam-writer [path ^SAMFileReader reader &
                  {:keys [sorted] :or {sorted false}}]
  (.makeSAMOrBAMWriter (SAMFileWriterFactory.)
                       (.getFileHeader reader)
                       sorted
                       (io/file path)))

;; Accessors
(defn read-name [^SAMRecord read]
  (.getReadName read))

(defn position [^SAMRecord read]
  (-> read
      .getAlignmentStart
      int
      dec))

(defn reference-name [^SAMRecord r]
  (.getReferenceName r))

(defn alignment-score [^SAMRecord read]
  (.getAttribute read "AS"))

(defn nm [^SAMRecord read]
  (.getAttribute read "NM"))

(defn primary? [^SAMRecord read]
  (not (.getNotPrimaryAlignmentFlag read)))

(defn mapped? [^SAMRecord read]
  (not (.getReadUnmappedFlag read)))

(defn- ^BitSet byte-array->uncertain-sites [^bytes xs]
  (let [l (alength xs)
        ^BitSet bs (BitSet. l)]
    (doseq [i (range l)]
      (let [b (aget xs i)
            uncertain? (and (not= b 100) (not= b 0))]
        (.set bs (int i) uncertain?)))
    bs))

;; Handle record base expected match
(defn ^BitSet uncertain-sites [^SAMRecord r]
  (-> r (.getAttribute TAG-EXP-MATCH) byte-array->uncertain-sites))

(def partition-by-name (partial partition-by read-name))
