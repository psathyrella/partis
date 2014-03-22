(ns ighutil.sam
  "Functions for working with SAM/BAM records"
  (:import [net.sf.samtools
            Defaults
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory
            SAMFileWriter]
           [java.util BitSet])
  (:require [clojure.java.io :as io]))

;;;;;;;;;;;;;;;;
;; SAM tag names
(def ^String TAG-N-MISMATCHES "NM")
(def ^String TAG-COUNT "XC")
(def ^String TAG-EXP-MATCH "bq")
(def ^String TAG-STATUS "XS")
;;;;;;;;;;;;;;;;

;;;;;;
;; I/O
;;;;;;
(defn ^SAMFileReader bam-reader [path]
  (-> path io/file SAMFileReader.))

(defn ^SAMFileWriter bam-writer [path ^SAMFileReader reader &
                                 {:keys [sorted compress]
                                  :or {sorted true
                                       compress true}}]
  (.makeBAMWriter (SAMFileWriterFactory.)
                  (.getFileHeader reader)
                  sorted
                  (io/file path)
                  (if compress Defaults/COMPRESSION_LEVEL 0)))

(defn reference-names [^SAMFileReader reader]
  "Get the names of the reference sequences in this file."
  (letfn [(sequence-name [^net.sf.samtools.SAMSequenceRecord x]
            (.getSequenceName x))]
    (->> reader
         .getFileHeader
         .getSequenceDictionary
         .getSequences
         (mapv sequence-name))))

;;;;;;;;;;;;;
;; Accessors
;;;;;;;;;;;;;
(defn read-name [^SAMRecord read]
  (.getReadName read))

(defn ^Integer position [^SAMRecord read]
  "*0-based* position of read along reference sequence"
  (-> read
      .getAlignmentStart
      int
      dec))

(defn ^String reference-name [^SAMRecord r]
  (.getReferenceName r))

(defn ^Integer alignment-score [^SAMRecord read]
  (.getIntegerAttribute read "AS"))

(defn ^Integer nm [^SAMRecord read]
  "Number of mismatches"
  (.getIntegerAttribute read "NM"))

(defn ^Integer sequence-status [^SAMRecord read]
  (.getIntegerAttribute read TAG-STATUS))

(defn ^Boolean primary? [^SAMRecord read]
  (not (.getNotPrimaryAlignmentFlag read)))

(defn ^Boolean mapped? [^SAMRecord read]
  (not (.getReadUnmappedFlag read)))

(defn ^Boolean supplementary? [^SAMRecord read]
  (.getSupplementaryAlignmentFlag read))

(defn exp-match [^SAMRecord read]
  "Matching probabilities, expressed as percentage"
  (.getByteArrayAttribute read TAG-EXP-MATCH))

;;;;;;;;;;;;;;;;;;
;; Site-certainty
;;;;;;;;;;;;;;;;;;
(defn- ^BitSet byte-array->uncertain-sites [^bytes xs]
  "Convert a byte array [e.g., from the bq tag] to a bitset with uncertain
   sites."
  (let [l (alength xs)
        ^BitSet bs (BitSet. l)]
    (doseq [i (range l)]
      (let [b (aget xs i)
            uncertain? (and (not= b 100) (not= b 0))]
        (.set bs (int i) uncertain?)))
    bs))

;; Handle record base expected match
(defn ^BitSet uncertain-sites [^SAMRecord r]
  "Returns a BitSet where set bits indicate that a site is certain"
  (-> r exp-match byte-array->uncertain-sites))

(defn ^Character ig-segment [^SAMRecord r]
  "Gene segment type (e.g. V / D / J)"
  (-> r reference-name (.charAt 3)))

(defn ^String ig-locus-segment [^SAMRecord r]
  "Gene locus and segment type (e.g. IGHV, IGHD, IGLJ)"
  (-> r reference-name (.substring 0 4)))

;;;;;;;;;;;;;;;;;;;;;;
;; Record partitioning
;;;;;;;;;;;;;;;;;;;;;;
(def partition-by-name-segment (partial partition-by (juxt read-name ig-segment)))
(def partition-by-name (partial partition-by read-name))
