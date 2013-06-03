(ns ighutil.imgt
  (:require [clojure.java.io :as io]
            [clojure.string :as string])
  (:import [net.sf.picard.reference FastaSequenceFile]))

(defn indexed
  "Returns a lazy sequence of [index, item] pairs, where items come
  from 's' and indexes count up from zero.

  (indexed '(a b c d))  =>  ([0 a] [1 b] [2 c] [3 d])"
  [s]
  (map vector (iterate inc 0) s))

(defn- nongap-indexed [sequence]
  (loop [result []
         sequence (indexed sequence)
         c 0]
    (if-let [[raw-idx base] (first sequence)]
      )
))

(defn- parse-fasta [^FastaSequenceFile f]
  (lazy-seq
   (when-let [ref (.nextSequence f)]
     (cons ref (parse-fasta f)))))

(defn read-fasta [file]
  (with-open [ref-file (FastaSequenceFile. (io/file file) false)]
    (doall (parse-fasta ref-file))))

(defn create-cysteine-map)
