(ns ighutil.subcommands.aligned-base-composition
  (:import [net.sf.samtools
            AlignmentBlock
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [cliopatra.command :refer [defcommand]]
            [hiphip.double :as dbl]
            [plumbing.core :refer [frequencies-fast distinct-fast]]
            [ighutil.io :as zio]
            [ighutil.imgt :refer [strip-allele]]
            [ighutil.sam :as sam]))

(def ^:private N 128)

(defn- block-composition [^bytes qbases ^AlignmentBlock b]
  (assert (> (alength qbases) 0) "No bases")
  (let [^doubles result (double-array N 0.0)
        start (dec (int (.getReadStart b)))
        length (int (.getLength b))]
    (doseq [i (range start (+ start length))]
      (dbl/ainc result (aget qbases i) 1.0))
    result))

(defn- pairwise-sum
  ([^doubles xs ^doubles ys]
     (dbl/afill! [x xs y ys] (+ x y))
     xs))

(defn- read-composition [^SAMRecord read]
  (let [qbases (.getReadBases read)]
    (->> read
         .getAlignmentBlocks
         (map (partial block-composition qbases))
         (reduce pairwise-sum))))

(defn- dbl-proportion [^doubles xs]
  (let [tot (dbl/asum xs)]
    (dbl/amap [x xs] (/ x tot))))

(defcommand aligned-base-composition
  "Aligned base composition"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn io/file]
               ["-o" "--out-file" "Target file"
                :parse-fn zio/writer
                :default (zio/writer "-")]]}
  (assert (.exists ^java.io.File in-file))
  (with-open [reader (SAMFileReader. ^java.io.File in-file)]
    (.setValidationStringency
     reader
     SAMFileReader$ValidationStringency/SILENT)
    (let [composition (->> reader
                           .iterator
                           iterator-seq
                           (filter (every-pred sam/mapped? sam/primary?))
                           (map read-composition)
                           (reduce pairwise-sum))
          props (dbl-proportion composition)
          rows (->> props
                    (map vector (map char (range N)) composition)
                    (remove (comp zero? second)))]
      (with-open [out-file ^java.io.Closeable out-file]
        (csv/write-csv out-file (concat
                                 [["base" "count" "proportion"]]
                                 rows))))))
