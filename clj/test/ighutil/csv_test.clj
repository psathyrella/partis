(ns ighutil.csv-test
  (:require [clojure.test :refer :all]
            [clojure.java.io :as io]
            [schema.test]
            [ighutil.csv :refer :all]))

(use-fixtures :once schema.test/validate-schemas)

(deftest test-read-csv
  (with-open [r (io/reader "data/naive_mutated.csv")]
    (let [rows (->> r csv-to-maps (take 20) vec)]
      (is (= {:reference "IGHV1-18" :position "1"}
             (first rows)))
      (is (= 20 (count rows))))))

(deftest test-read-typed-csv
  (with-open [r (io/reader "data/naive_mutated.csv")]
    (let [rows (->> (read-typed-csv
                     r
                     {:position int-of-string})
                    (take 20)
                    vec)]
      (is (= {:reference "IGHV1-18" :position 1}
             (first rows)))
      (is (= 20 (count rows))))))

(deftest test-of-string
  (are [x y] (= x y)
       0 (int-of-string "0")
       3.14 (double-of-string "3.14")
       1.0 (float-of-string "1.0")))
