(ns ighutil.ubtree
  (:require [clojure.set :as set]))

(defprotocol UBTreeProtocol
  (insert [this s]
    "Inserts ordered set *s* into the tree.")
  (lookup-first [this q]
    "Looks up whether *any* sets in `tree` are a subset of the ordered set `q`")
  (lookup-subs [this q]
    "Finds all sets which are subsets of ordered set `q` in *tree*"))

(defprotocol UBNodeProtocol
  (node-insert [this q] (assert false)))

(def ^{:private true} not-nil? (complement nil?))

(defrecord UBNode [element end-of-path sons]
  UBNodeProtocol
  (node-insert [_ q]
    (let [qi (first q)
          r (rest q)
          more (empty? r)
          c (get sons qi (UBNode. qi end-of-path {}))]
      (if (empty? q)
        (UBNode. element true sons)
        (UBNode. element end-of-path (assoc sons qi (node-insert c r)))))))

(defn- ensure-non-empty [s]
  (when-not (seq s)
    (throw (IllegalArgumentException. "Empty sequence"))))

(defrecord UBTree [forest]
  UBTreeProtocol
  (insert [_ s]
    (ensure-non-empty s)
    (let [qi (first s)
          r (rest s)
          more (not (empty? r))]
      (let [n (get forest qi (UBNode. qi false {}))]
        (UBTree. (assoc forest qi (node-insert n r))))))
  (lookup-first [_ q]
    (ensure-non-empty q)
    (assert false))
  (lookup-subs [_ q]
    (ensure-non-empty q)
    (assert false)))
