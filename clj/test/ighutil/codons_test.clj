(ns ighutil.codons-test
  (:import [io.github.cmccoy.dna Codons])
  (:require [clojure.test :refer :all]))

(deftest test-translation
  (letfn [(translate-codon [^String s]
            (let [[a b c] (.getBytes s)]
              (Codons/translateCodon a b c)))]
    (are [x y] (= (byte x) (translate-codon y))
         \Y "TAT"
         \* "TGA"
         \K "AAR"
         \I "ATH"
         \X "NNN")))
