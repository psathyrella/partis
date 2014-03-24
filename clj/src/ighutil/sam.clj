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
  [reader :- SAMFileReader]
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
(s/defn read-name :- String [r :- SAMRecord]
  (.getReadName r))

(s/defn position :- s/Int [r :- SAMRecord]
  "*0-based* position of read along reference sequence"
  (-> r
      .getAlignmentStart
      int
      dec))

(s/defn reference-name :- String [r :- SAMRecord]
  (.getReferenceName r))

(s/defn alignment-score :- s/Int [r :- SAMRecord]
  (.getIntegerAttribute r "AS"))

(s/defn nm :- s/Int [r :- SAMRecord]
  "Number of mismatches"
  (.getIntegerAttribute r "NM"))

(s/defn sequence-status :- s/Int [r :- SAMRecord]
  (.getIntegerAttribute r TAG-STATUS))

(s/defn primary? :- s/Bool [r :- SAMRecord]
  (not (.getNotPrimaryAlignmentFlag r)))

(s/defn mapped? :- s/Bool [r :- SAMRecord]
  (not (.getReadUnmappedFlag r)))

(s/defn supplementary? :- s/Bool [r :- SAMRecord]
  (.getSupplementaryAlignmentFlag r))

(s/defn exp-match :- bytes [r :- SAMRecord]
  "Matching probabilities, expressed as percentage"
  (.getByteArrayAttribute r TAG-EXP-MATCH))

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
(s/defn uncertain-sites :- BitSet [r :- SAMRecord]
  "Returns a BitSet where set bits indicate that a site is certain"
  (-> r exp-match byte-array->uncertain-sites))

(s/defn ig-segment :- Character [r :- SAMRecord]
  "Gene segment type (e.g. V / D / J)"
  (-> r reference-name (.charAt 3)))

(s/defn ig-locus-segment :- String [r :- SAMRecord]
  "Gene locus and segment type (e.g. IGHV, IGHD, IGLJ)"
  (-> r reference-name (.substring 0 4)))

;;;;;;;;;;;;;;;;;;;;;;
;; Record partitioning
;;;;;;;;;;;;;;;;;;;;;;
(def partition-by-name-segment (partial partition-by (juxt read-name ig-segment)))
(def partition-by-name (partial partition-by read-name))
