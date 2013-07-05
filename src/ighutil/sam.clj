(ns ighutil.sam
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory])
  (:require [clojure.java.io :as io]))

(defn bam-reader [path]
  (-> path io/file SAMFileReader.))

(defn bam-writer [path ^SAMFileReader reader &
                  {:keys [sorted] :or {sorted false}}]
  (.makeSAMOrBAMWriter (SAMFileWriterFactory.)
                       (.getFileHeader reader)
                       sorted
                       (io/file path)))

(defn read-name [^SAMRecord read]
  (.getReadName read))

(def partition-by-name (partial
                        partition-by
                        read-name))

;; Accessors
(defn read-name [^SAMRecord read]
  (.getReadName read))

(defn alignment-score [^SAMRecord read]
  (.getAttribute read "AS"))

(defn nm [^SAMRecord read]
  (.getAttribute read "NM"))

(defn primary? [^SAMRecord read]
  (not (.getNotPrimaryAlignmentFlag read)))

(defn mapped? [^SAMRecord read]
  (not (.getReadUnmappedFlag read)))
