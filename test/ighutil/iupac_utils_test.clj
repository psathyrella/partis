(ns ighutil.iupac-utils-test
  (:import [io.github.cmccoy.dna IUPACUtils])
  (:require [clojure.test :refer :all]))

(defn roundtrip [^String s]
  (let [b (.getBytes s)
        packed (IUPACUtils/packBytes b)
        unpacked (IUPACUtils/unpackBytes packed)
        unpacked-str (String. unpacked)]
    (is (= s unpacked-str))))

;; All IUPAC
(def ALL "ACGTMRWSYKVHDBN")

(deftest test-roundtrip
  (roundtrip "NNNN")
  (roundtrip ALL))

(defn compatible? [^String s1 ^String s2]
  (let [p1 (IUPACUtils/packBytes (.getBytes s1))
        p2 (IUPACUtils/packBytes (.getBytes s2))]
    (IUPACUtils/areCompatible p1 p2)))

(deftest test-are-compatible
  (are [x y] (compatible? x y)
       "NNNN" "NNNN"
       ALL ALL
       "NNNN" "AAAA"
       "AAAA" "NNNN"
       "R" "D")
  (are [x y] (not (compatible? x y))
       "ACGT" "TGCA"
       "V" "T"
       "C" "D"
       "R" "Y"))

(defn subset? [^String s1 ^String s2]
  (let [p1 (IUPACUtils/packBytes (.getBytes s1))
        p2 (IUPACUtils/packBytes (.getBytes s2))]
    (IUPACUtils/isSubset p1 p2)))

(deftest test-is-subset
  (are [x y] (subset? x y)
       "ACGT" "ACGT"
       "NNNN" "ACGT"
       "WCGB" "ACGT")
  (are [x y] (not (subset? x y))
       "M" "R"
       "Y" "N"))
