(ns ighutil.sam
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory
            SAMFileWriter]
           [java.util BitSet])
  (:require [clojure.java.io :as io]
            [ighutil.sam-tags :refer [TAG-EXP-MATCH TAG-STATUS]]))

(defn bam-reader [path]
  (-> path io/file SAMFileReader.))

(defn ^SAMFileWriter bam-writer [path ^SAMFileReader reader &
                  {:keys [sorted]
                   :or {sorted true}}]
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

(defn ^String reference-name [^SAMRecord r]
  (.getReferenceName r))

(defn ^Integer alignment-score [^SAMRecord read]
  (.getAttribute read "AS"))

(defn ^Integer nm [^SAMRecord read]
  (.getAttribute read "NM"))

(defn ^Integer sequence-status [^SAMRecord read]
  (.getAttribute read TAG-STATUS))

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

(defn gene-type [^SAMRecord r]
  (-> r reference-name (.charAt 3)))

(def partition-by-name-type (partial partition-by (juxt read-name gene-type)))
(def partition-by-name (partial partition-by read-name))
