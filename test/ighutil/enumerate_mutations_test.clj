(ns ighutil.enumerate-mutations-test
  (:require [clojure.test :refer :all]
            [ighutil.enumerate-mutations :refer :all]))

(deftest test-summarize-mutation-partition
  (testing "Simple case"
    (is (= {:v-gene nil
            :j-gene nil
            :cdr3-length nil
            :edges {"3" "1", "2" "1"}}
           (summarize-mutation-partition
               [{:name "1" :mutations [{:ref-idx 101 :ref \A :qry \C}]}
                {:name "2" :mutations [{:ref-idx 101 :ref \A :qry \C}]}
                {:name "3" :mutations [{:ref-idx 101 :ref \A :qry \C}
                                       {:ref-idx 105 :ref \G :qry \T}]}])))))
