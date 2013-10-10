(ns ighutil.sam-utils-test
  (:import [net.sf.samtools SAMRecord SAMFileHeader SAMSequenceRecord]
           [io.github.cmccoy.sam SAMUtils AlignedPair$MatchesReference])
  (:require [clojure.test :refer :all]))


(def REF-MAP {"sequence1" (.getBytes "ACGTA")})

(def HEADER (doto (SAMFileHeader.)
              (.addSequence (SAMSequenceRecord. "sequence1" 5))))

(defn- create-read [^String ref-name ^String bases ^String cigar pos & r]
  (doto (SAMRecord. HEADER)
    (.setReadBases (.getBytes bases))
    (.setCigarString cigar)
    (.setAlignmentStart (int pos))
    (.setAttribute "bq" (or (first r) (byte-array (count bases) (byte 100))))
    (.setReferenceName ref-name)))

(deftest test-aligned-pairs
  (let [bq (->> [100 50 0 100 0 0] (map byte) byte-array)
        read (create-read "sequence1" "ACGGCT" "3M1D1M2I" 2 bq)]
    (let [r (SAMUtils/getAlignedPairs read)
          F AlignedPair$MatchesReference/FALSE
          T AlignedPair$MatchesReference/TRUE
          U AlignedPair$MatchesReference/UNKNOWN]
      (is (= 7 (count r)))
      (is (= [T U F F T F F] (map #(.getMatchesReference %) r)))
      (is (= [0 1 2 -1 3 4 5] (map #(.getQueryPosition %) r)))
      (is (= [1 2 3 4 5 -1 -1] (map #(.getReferencePosition %) r))))))
