(ns ighutil.subcommands.kmer-matrix-test
  (:import [net.sf.samtools SAMRecord SAMSequenceRecord SAMFileHeader])
  (:require [clojure.test :refer :all]
            [ighutil.subcommands.kmer-matrix :refer [kmer-mutations]]))

(def SEQUENCE1 (.getBytes "ACGT"))

(def HEADER (doto (SAMFileHeader.)
              (.addSequence (SAMSequenceRecord. "sequence1" (alength SEQUENCE1)))))

(defn- create-read [^String ref-name ^String bases ^String cigar pos]
  (doto (SAMRecord. HEADER)
    (.setReadBases (.getBytes bases))
    (.setCigarString cigar)
    (.setAttribute "bq" (byte-array (.length bases) (byte 100)))
    (.setAlignmentStart (int pos))
    (.setReferenceName ref-name)))

(defn- interval-tree [coll]
  (let [r (net.sf.picard.util.IntervalTree.)]
    (doseq [[start end] coll]
      (.put r start end :ignore))
    r))

(deftest test-kmer-mutations
  (let [read (create-read "sequence1" "ACCTT" "4M1S" 1)]
    (is (= [["AC" "AC"] ["CG" "CC"] ["GT" "CT"]] (kmer-mutations 2 read SEQUENCE1)))
    (is (= [["AC" "AC"] ["CG" "CC"]]
           (kmer-mutations
            2
            read
            SEQUENCE1
            :exclude (interval-tree [[4 4]]))))
    (is (= [] (kmer-mutations 2 read SEQUENCE1 :exclude (interval-tree [[1 4]]))))))

(deftest test-kmer-mutations-uncertain
  (let [read (partial create-read "sequence1" "ACCTT" "4M1S" 1)]
    (is (= [["AC" "AC"] ["CG" "CC"] ["GT" "CT"]] (kmer-mutations 2 (read) SEQUENCE1
                                                                 :drop-uncertain? true)))
    (is (= [["AC" "AC"] ["CG" "CC"]]
           (kmer-mutations
            2
            (doto (read) (.setAttribute "bq" (byte-array (map byte [100 100 100 75]))))
            SEQUENCE1
            :drop-uncertain? true)))
    (is (= []
           (kmer-mutations
            2
            (doto (read) (.setAttribute "bq" (byte-array (map byte [99 99 99 99]))))
            SEQUENCE1
            :drop-uncertain? true)))))

(deftest test-kmer-mutations-frame
  (let [read (create-read "sequence1" "ACCTT" "4M1S" 1)]
    (is (= [["AC" "AC"] ["GT" "CT"]] (kmer-mutations 2 read SEQUENCE1 :frame 0)))
    (is (= [["CG" "CC"]] (kmer-mutations 2 read SEQUENCE1 :frame 1)))))
