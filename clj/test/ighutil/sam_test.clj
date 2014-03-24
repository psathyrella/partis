(ns ighutil.sam-test
  (:require [clojure.test :refer :all]
            [schema.test]
            [ighutil.sam :refer :all]))

(use-fixtures :once schema.test/validate-schemas)

(def ^:private BAM "testdata/test.bam")

(defmacro with-sam-reader [bindings & body]
  (assert (vector? bindings))
  (assert (= 2 (count bindings)))
  `(with-open [~(first bindings) (bam-reader ~(second bindings))]
     ~@body))

(deftest test-ig-segment
  (with-sam-reader [reader BAM]
    (let [records (-> reader .iterator iterator-seq)]
      (doseq [record (take 10 records)]
        (is (= "IGHV" (ig-locus-segment record)))))))

(deftest test-accessors
  (with-sam-reader [reader BAM]
    (let [record (-> reader .iterator iterator-seq first)]
      (is (= 12 (nm record)))
      (is (= 39 (alignment-score record)))
      (is (= "04-A-M_0000013" (read-name record)))
      (is (= 231 (position record)))
      (is (= [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
              100 100 0 100 100 0 100 100 100 0 0 100 100 0 100 100
              100 0 100 100 100 100 0 100 100 100 0 100 100 100 0 100
              100 100 100 100 100 100 100 100 100 100 100 100 100 100
              100 0 100 100 100 100 0 100 100 100 0 100 100 100 100
              100 100 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
              -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
              -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
             ((comp vec exp-match) record)))
      (is (= "IGHV1-18*01" (reference-name record)))
      (is (= false (supplementary? record))))))

(deftest test-uncertain-sites
  (with-sam-reader [reader BAM]
    (let [records (-> reader .iterator iterator-seq)
          uncertain (map uncertain-sites records)
          card (map #(.cardinality %) uncertain)]
      (is (= [67 66] (take 2 card))))))
