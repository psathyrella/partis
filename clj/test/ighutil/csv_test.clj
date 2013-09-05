(ns ighutil.csv-test
  (:require [clojure.test :refer :all]
            [clojure.java.io :as io]
            [ighutil.csv :refer :all]))

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
