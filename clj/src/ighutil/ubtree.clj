(ns ighutil.ubtree
  "Unlimited Branching Tree implementation
   See:
   Hoffman and Koeler 'A New Method to Index and Query Sets'
   http://www.ijcai.org/Past%20Proceedings/IJCAI-99-VOL-1/PDF/067.pdf

   **Note that everything here assumes sorted input.**"
  (:require [clojure.set :as set]))

(defprotocol UBTreeProtocol
  (insert [this s]
    "Inserts ordered set *s* into the tree.")
  (lookup-first [this q]
    "Looks up whether *any* sets in `tree` are a subset of the ordered set `q`")
  (lookup-subs [this q]
    "Finds all sets which are subsets of ordered set `q` in `tree`"))

(defprotocol UBNodeProtocol
  (node-insert [this q]
    "See UBTreeProtocol.insert")
  (node-lookup-first [this q]
    "See UBTreeProtocol.lookup-first")
  (node-lookup-subs [this q elems]
    "See UBTreeProtocol.lookup-subs"))

(defn suffixes [coll]
  "Returns all suffixes of coll, including coll itself."
  (let [c (vec coll)]
    (for [j (range (count c))]
      (subvec c j))))

(defrecord UBNode [element end-of-path sons]
  UBNodeProtocol
  (node-insert [_ q]
    (let [qi (first q)
          r (rest q)
          more (empty? r)
          c (get sons qi (->UBNode qi more {}))]
      (if (empty? q)
        (->UBNode element true sons)
        (->UBNode element end-of-path (assoc sons qi (node-insert c r))))))
  (node-lookup-first [_ q]
    (or end-of-path
        (let [s (suffixes q)
              f (fn [[qi & r]]
                  (if-let [n (get sons qi)]
                    (node-lookup-first n r)))]
          (some f s))))
  (node-lookup-subs [_ q elems]
    (concat
     (when end-of-path [elems])
     (let [s (suffixes q)
           f (fn [[qi & r]]
               (if-let [n (get sons qi)]
                 (node-lookup-subs n r (conj elems qi))))]
       (mapcat f s)))))

(defn- ensure-non-empty [s]
  (when-not (seq s)
    (throw (IllegalArgumentException. "Empty sequence"))))

(defrecord UBTree [forest]
  UBTreeProtocol
  (insert [_ s]
    (ensure-non-empty s)
    (let [[qi & r] s
          more (not (empty? r))]
      (let [n (get forest qi (->UBNode qi false {}))]
        (->UBTree (assoc forest qi (node-insert n r))))))
  (lookup-first [_ q]
    (ensure-non-empty q)
    (let [s (suffixes q)
          f (fn [[qi & r]]
              (if-let [n (get forest qi)]
                (node-lookup-first n r)))]
      (some f s)))
  (lookup-subs [_ q]
    (ensure-non-empty q)
    (let [s (suffixes q)
          f (fn [[qi & r]]
              (if-let [n (get forest qi)]
                (node-lookup-subs n r [qi])))]
      (mapcat f s))))

(def ubtree (->UBTree {}))
