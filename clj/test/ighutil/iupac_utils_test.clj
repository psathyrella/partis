(ns ighutil.iupac-utils-test
  (:import [io.github.cmccoy.dna IUPACUtils])
  (:require [clojure.test :refer :all]))


;; All IUPAC
(def ALL "ACGTMRWSYKVHDBN")

(deftest test-roundtrip
  (letfn [(roundtrip [^String s]
            (let [b (.getBytes s)
                  packed (IUPACUtils/packBytes b)
                  unpacked (IUPACUtils/unpackBytes packed)
                  unpacked-str (String. unpacked)]
              (is (= s unpacked-str))))]
    (roundtrip "NNNN")
    (roundtrip ALL)))

(deftest test-are-compatible
  (letfn [(compatible? [^String s1 ^String s2]
            (let [p1 (IUPACUtils/packBytes (.getBytes s1))
                  p2 (IUPACUtils/packBytes (.getBytes s2))]
              (IUPACUtils/areCompatible p1 p2)))]
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
         "R" "Y")))

(deftest test-is-subset
  (letfn [(subset? [^String s1 ^String s2]
            (let [p1 (IUPACUtils/packBytes (.getBytes s1))
                  p2 (IUPACUtils/packBytes (.getBytes s2))]
              (IUPACUtils/isSubset p1 p2)))]
    (are [x y] (subset? x y)
         "ACGT" "ACGT"
         "NNNN" "ACGT"
         "WCGB" "ACGT")
    (are [x y] (not (subset? x y))
         "M" "R"
         "Y" "N")))

(deftest test-kmer-pack-roundtrip
  (letfn [(roundtrip [^String s]
            (-> s
                .getBytes
                IUPACUtils/packBytes
                IUPACUtils/packKmer
                IUPACUtils/unpackKmer
                IUPACUtils/unpackBytes
                String.))]
    (are [x] (= x (roundtrip x))
         "ACGT"
         "NNNN"
         "WCGB"
         "AAAAA")))

(deftest test-disambiguate
  (letfn [(disambiguate [^String s]
            (->> s
                .getBytes
                IUPACUtils/packBytes
                IUPACUtils/disambiguate
                seq
                (mapv #(String. ^bytes %))))]
    (are [x y] (= x (disambiguate y))
         ["ACGT"] "ACGT"
         ["CCGT" "TCGT"] "YCGT"
         ["CCAT" "CCGT" "TCAT" "TCGT"] "YCRT")
    (is (= (int (Math/pow 4 4)) (count (disambiguate "NNNN"))))))
