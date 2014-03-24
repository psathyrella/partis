(ns ighutil.sam
  "Functions for working with SAM/BAM records"
  (:require [clojure.java.io :as io]
            [schema.core :as s])
  (:import [net.sf.samtools
            Defaults
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileWriterFactory
            SAMFileWriter]
           [java.util BitSet]))

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
(s/defn bam-reader :- SAMFileReader
  [path]
  (-> path io/file SAMFileReader.))

(defn ^SAMFileWriter bam-writer
  [path
   ^SAMFileReader reader &
   {:keys [sorted compress]
    :or {sorted true
         compress true}}]
  (.makeBAMWriter (SAMFileWriterFactory.)
                  (.getFileHeader reader)
                  sorted
                  (io/file path)
                  (if compress Defaults/COMPRESSION_LEVEL 0)))

(s/defn reference-names :- [String]
  "Get the names of the reference sequences in this file."
  [reader :- SAMFileReader]
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
(s/defn read-name :- String [r :- SAMRecord]
  (.getReadName r))

(s/defn position :- s/Int
  "*0-based* position of read along reference sequence"
  [r :- SAMRecord]
  (-> r
      .getAlignmentStart
      int
      dec))

(s/defn reference-name :- String [r :- SAMRecord]
  (.getReferenceName r))

(s/defn alignment-score :- s/Int [r :- SAMRecord]
  (.getIntegerAttribute r "AS"))

(s/defn nm :- s/Int
  "Number of mismatches"
  [r :- SAMRecord]
  (.getIntegerAttribute r "NM"))

(s/defn sequence-status :- s/Int [r :- SAMRecord]
  (.getIntegerAttribute r TAG-STATUS))

(s/defn primary? :- s/Bool [r :- SAMRecord]
  (not (.getNotPrimaryAlignmentFlag r)))

(s/defn mapped? :- s/Bool [r :- SAMRecord]
  (not (.getReadUnmappedFlag r)))

(s/defn supplementary? :- s/Bool [r :- SAMRecord]
  (.getSupplementaryAlignmentFlag r))

(s/defn exp-match
  "Matching probabilities, expressed as percentage"
  [r :- SAMRecord]
  (.getByteArrayAttribute r TAG-EXP-MATCH))

;;;;;;;;;;;;;;;;;;
;; Site-certainty
;;;;;;;;;;;;;;;;;;
(defn- ^java.util.BitSet byte-array->uncertain-sites
  "Convert a byte array [e.g., from the bq tag] to a bitset with uncertain
   sites."
  [^bytes xs]
  (let [l (alength xs)
        ^java.util.BitSet bs (java.util.BitSet. l)]
    (doseq [i (range l)]
      (let [b (aget xs i)
            uncertain? (and (not= b 100) (not= b 0))]
        (.set bs (int i) uncertain?)))
    bs))

;; Handle record base expected match
(s/defn uncertain-sites :- java.util.BitSet
  "Returns a BitSet where set bits indicate that a site is certain"
  [r :- SAMRecord]
  (-> r exp-match byte-array->uncertain-sites))

(s/defn ig-segment :- Character
  "Gene segment type (e.g. V / D / J)"
  [r :- SAMRecord]
  (-> r reference-name (.charAt 3)))

(s/defn ig-locus-segment :- String
  "Gene locus and segment type (e.g. IGHV, IGHD, IGLJ)"
  [r :- SAMRecord]
  (-> r reference-name (.substring 0 4)))

;;;;;;;;;;;;;;;;;;;;;;
;; Record partitioning
;;;;;;;;;;;;;;;;;;;;;;
(def partition-by-name-segment (partial partition-by (juxt read-name ig-segment)))
(def partition-by-name (partial partition-by read-name))
