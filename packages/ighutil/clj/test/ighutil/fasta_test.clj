(ns ighutil.fasta-test
  (:require [clojure.test :refer :all]
            [ighutil.fasta :refer :all]))

(def ^:private FASTA "data/ighv.fasta")

(defn- er [x]
  (letfn [(tostr [[name seq]] [name (String. ^bytes seq)])]
    (->> x extract-references (mapv tostr))))

(deftest test-read-fasta
  (is (= (er FASTA)
         (er (java.io.File. ^String FASTA))))
  (let [r (er FASTA)]
    (is (= ["IGHV1-18*01" "IGHV1-18*03" "IGHV1-2*01"]
           (->> r (take 3) (mapv first))))
    (is (= "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA"
           (-> r first second)))))
