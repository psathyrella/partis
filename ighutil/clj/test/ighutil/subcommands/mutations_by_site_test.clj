(ns ighutil.subcommands.mutations-by-site-test
  (:import [net.sf.samtools SAMRecord SAMFileHeader SAMSequenceRecord])
  (:require [ighutil.subcommands.mutations-by-site :refer :all]
            [clojure.test :refer :all]))


(def REF-MAP {"sequence1" (.getBytes "ACGTA")})

(def HEADER (doto (SAMFileHeader.)
              (.addSequence (SAMSequenceRecord. "sequence1" 5))))

(defn- create-read [^String ref-name ^String bases ^String cigar pos]
  (doto (SAMRecord. HEADER)
    (.setReadBases (.getBytes bases))
    (.setCigarString cigar)
    (.setAlignmentStart (int pos))
    (.setAttribute "bq" (byte-array (count bases) (byte 100)))
    (.setReferenceName ref-name)))

(deftest test-count-bases-by-position
  (let [reads [(create-read "sequence1" "ACG" "3M" 1)
               (create-read "sequence1" "CGCT" "3M1D1M" 1)]
        result (vec (count-bases-by-position reads REF-MAP))
        exp {:position 0
             :exp-match 2.0
             :A 1
             :C 1
             :G 0
             :T 0
             :n-reads 2
             :reference "sequence1"
             :ref-base \A
             :alignment-position nil}]
    (is (= 5 (count result)))
    (is (= [2 2 2 0 1] (map :n-reads result)))
    (is (= (seq "ACGTA") (map :ref-base result)))
    (is (= exp (first result)))))
