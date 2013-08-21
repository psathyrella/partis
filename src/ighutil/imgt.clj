(ns ighutil.imgt
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clojure.edn :as edn]
            [cheshire.core :as cheshire]
            [flatland.useful.seq :refer [indexed]])
  (:import [net.sf.picard.reference FastaSequenceFile ReferenceSequence]))

(def GAP-CHARS #{\. \-})

(defn strip-allele [^String s]
  "Remove the allele from a string"
  (let [idx (.lastIndexOf s 42)]  ; 42 == '*'
    (if (< idx 0)
      s
      (.substring s 0 idx))))

(defn- nongap-lookup [sequence]
  (loop [result []
         sequence (indexed sequence)
         c 0]
    (if-let [[raw-idx base] (first sequence)]
      (if (GAP-CHARS base)
        (recur result (rest sequence) c)
        (recur (conj result [c raw-idx]) (rest sequence) (inc c)))
      (into {} result))))


(defn- parse-fasta [lines]
  "Read a FASTA file"
  (let [fasta-header? (fn [^String s] (.startsWith s ">"))]
    (->> lines
         (remove string/blank?)
         (partition-by fasta-header?)
         (partition 2)
         (map (fn [[[^String name] seqs]]
                [(.substring name 1)
                 (string/join (apply concat seqs))])))))

(defn- parse-fasta-picard [^FastaSequenceFile f]
  (lazy-seq
   (when-let [ref (.nextSequence f)]
     (cons ref (parse-fasta f)))))

(defn read-fasta [file]
  (with-open [ref-file (FastaSequenceFile. (io/file file) false)]
    (doall (parse-fasta-picard ref-file))))

(def CYSTEINE-POSITION (* 3 (- 104 1)))

(defn create-v-meta
  ([] (create-v-meta "ighutil/ighv_aligned.fasta"))
  ([^String resource-path]
     (with-open [reader (->> resource-path io/resource io/reader)]
       (let [extract-position (fn [[^String name ^String sequence]]
                                (let [pos-map (nongap-lookup sequence)]
                                  [(second (.split name "\\|"))
                                   {:cysteine-position
                                    (get pos-map CYSTEINE-POSITION)
                                    :translation pos-map
                                    :sequence (.replaceAll sequence "[.-]" "")}]))]
         (->> reader
              line-seq
              parse-fasta
              (map extract-position)
              (into {}))))))

(defn slurp-edn-resource [resource-name]
  (with-open [reader (-> resource-name
                         io/resource
                         io/reader
                         (java.io.PushbackReader.))]
    (edn/read reader)))

(def v-gene-meta (delay (create-v-meta)))


(defn dump-v-gene-meta [fname]
  (->> @v-gene-meta
       (map #(update-in % [1 :translation]
                        (comp (partial apply mapv vector) sort)))
       (into {})
       cheshire/encode
       (spit fname)))
