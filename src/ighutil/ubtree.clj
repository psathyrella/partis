(ns ighutil.ubtree
  "UBTree implementation
   See:
   Hoffman and Koeler 'A New Method to Index and Query Sets'
   http://www.ijcai.org/Past%20Proceedings/IJCAI-99-VOL-1/PDF/067.pdf"
  (:require [clojure.set :as set]))

(defprotocol UBTreeProtocol
  (insert [this s]
    "Inserts ordered set *s* into the tree.")
  (lookup-first [this q]
    "Looks up whether *any* sets in `tree` are a subset of the ordered set `q`")
  (lookup-subs [this q]
    "Finds all sets which are subsets of ordered set `q` in *tree*"))

(defprotocol UBNodeProtocol
  (node-insert [this q])
  (node-lookup-first [this q]
    "See lookup-first")
  (node-lookup-subs [this q]
    "See lookup-subs"))

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
        (UBNode. element end-of-path (assoc sons qi (node-insert c r))))))
  (node-lookup-first [_ q]
    (let [qi (first q)
          r (rest q)
          c (get sons qi)]
      (cond
       end-of-path true
       (empty? q) false
       (not c) (recur r)
       :else (node-lookup-first c r))))
  (node-lookup-subs [_ q]
    (assert false)))

(defn- ensure-non-empty [s]
  (when-not (seq s)
    (throw (IllegalArgumentException. "Empty sequence"))))

(defrecord UBTree [forest]
  UBTreeProtocol
  (insert [_ s]
    (ensure-non-empty s)
    (let [[qi & r] s
          more (not (empty? r))]
      (let [n (get forest qi (UBNode. qi false {}))]
        (UBTree. (assoc forest qi (node-insert n r))))))
  (lookup-first [_ q]
    (ensure-non-empty q)
    (loop [[qi & r] q]
      (let [n (get forest qi)]
        (cond
         (and n (node-lookup-first n r)) true
         (empty? r) false
         :else (recur r)))))
  (lookup-subs [_ q]
    (ensure-non-empty q)
    (let [[qi & r] q])
    (assert false)))

(def ubtree (UBTree. {}))
