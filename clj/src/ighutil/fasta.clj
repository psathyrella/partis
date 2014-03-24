(ns ighutil.fasta
  "Simple FASTA parsing - loads all sequences into memory."
  (:require [clojure.java.io :as io]
            [schema.core :as s])
  (:import [net.sf.picard.reference
            ReferenceSequenceFile
            ReferenceSequenceFileFactory
            ReferenceSequence]))

(defprotocol ExtractRefable
  (extract-references [x]
    "Extract all references into a seq of [name, bases]
     Bases are encoded in a byte array."))

(extend-protocol ExtractRefable
  String
  (extract-references [^String x]
    (-> x io/file extract-references))
  java.io.File
  (extract-references [^java.io.File x]
    (-> x
        ReferenceSequenceFileFactory/getReferenceSequenceFile
        extract-references))
  ReferenceSequenceFile
  (extract-references [^ReferenceSequenceFile x]
    (.reset x)
    (->> (repeatedly #(.nextSequence x))
         (take-while identity)
         (map (fn [^ReferenceSequence s]
                [(.getName s) (.getBases s)])))))
