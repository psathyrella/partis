(ns ighutil.ubtree-test
  (:require [clojure.test :refer :all]
            [ighutil.ubtree :refer :all]))

(def tree (-> ubtree
              (insert [:a :b :f])
              (insert [:d :e])))

(deftest test-suffixes
  (are [x y] (= x (suffixes y))
       [] []
       [[:a]] [:a]
       [[:a :b] [:b]] '(:a :b)
       [[:a :c :e] [:c :e] [:e]] [:a :c :e]))

(deftest test-ubtree-construction
  (testing "Add one item"
    (are [x] (= [:a] (-> ubtree (insert x) :forest keys))
         [:a]
         [:a :b :c]
         [:a :d :f]))
  (testing "Add multiple items with same prefix"
    (is (= [:a] (-> ubtree (insert [:a :b]) (insert [:a :c]) :forest keys))))
  (testing "Add multiple keys with different prefixes"
    (is (= [:a :g] (-> ubtree
                       (insert [:a :b])
                       (insert [:g :z])
                       :forest
                       keys
                       sort)))))

(deftest test-ubtree-lookup-first
  (testing "Look up exact set"
    (are [x] (lookup-first tree x)
         [:a :b :f]
         [:d :e]))
  (testing "Look up super set"
    (are [x] (lookup-first tree x)
         [:a :b :c :f]
         [:a :b :c :f :z]
         [:a :b :d :e]))
  (testing "Look-up non super-set"
    (are [x] (not (lookup-first tree x))
         [:a :c]
         [:zzz]
         [:d :f])))

(deftest test-ubtree-lookup-subs
  (testing "Look up exact set"
    (are [x y] (= x (lookup-subs tree y))
         [[:a :b :f]] [:a :b :c :f]
         [[:a :b :f]] [:a :b :c :f :z]
         [[:d :e]] [:a :c :d :e :z :y])))
