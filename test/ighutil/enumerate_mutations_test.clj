(ns ighutil.enumerate-mutations-test
  (:import [io.github.cmccoy.sam Mutation Mutation$MutationType])
  (:require [clojure.test :refer :all]
            [ighutil.enumerate-mutations :refer :all]))

(defn- mut [ref-idx wt mut]
  (Mutation. Mutation$MutationType/MUTATION ref-idx wt mut))

(deftest test-summarize-mutation-partition
  (testing "Simple case"
    (is (= {:v-gene nil
            :j-gene nil
            :cdr3-length nil
            :edges {"3" '(["1" 1]),  "2" '(["1" 0])}
            :n-seqs 3}
           (summarize-mutation-partition
               [{:name "1" :mutations [(mut 101 "A" "C")]}
                {:name "2" :mutations [(mut 101 "A" "C")]}
                {:name "3" :mutations [(mut 101 "A" "C")
                                       (mut 105 "G" "T")]}])))))
