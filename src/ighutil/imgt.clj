(ns ighutil.imgt
  (:require [clojure.java.io :as io]
            [clojure.string :as string])
  (:import [net.sf.picard.reference FastaSequenceFile ReferenceSequence]))

(defn indexed
  "Returns a lazy sequence of [index, item] pairs, where items come
  from 's' and indexes count up from zero.

  (indexed '(a b c d))  =>  ([0 a] [1 b] [2 c] [3 d])"
  [s]
  (map vector (iterate inc 0) s))

(def GAP-CHARS # {\. \-})

(defn- nongap-lookup [sequence]
  (loop [result []
         sequence (indexed sequence)
         c 0]
    (if-let [[raw-idx base] (first sequence)]
      (if (GAP-CHARS base)
        (recur result (rest sequence) c)
        (recur (conj result [raw-idx c]) (rest sequence) (inc c)))
      (into {} result))))

(defn- parse-fasta [^FastaSequenceFile f]
  (lazy-seq
   (when-let [ref (.nextSequence f)]
     (cons ref (parse-fasta f)))))

(defn read-fasta [file]
  (with-open [ref-file (FastaSequenceFile. (io/file file) false)]
    (doall (parse-fasta ref-file))))

(def CYSTEINE-POSITION (* 3 (- 104 1)))

(defn create-cysteine-map []
  (let [f (-> "ighutil/ighv_aligned.fasta" io/resource io/file)
        extract-position (fn [^ReferenceSequence s]
                           (let [name (second (.. s (getName) (split "\\|")))
                                 pos-map (->> s
                                          (.getBases)
                                          (String.)
                                          nongap-lookup)]
                             [name (get pos-map CYSTEINE-POSITION)]))]
    (into {} (map extract-position (read-fasta f)))))
