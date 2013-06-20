(ns ighutil.fasta
  (:import [net.sf.picard.reference
            ReferenceSequenceFile
            ReferenceSequence]))



(defn encode-bases [^bytes bases])

(defn extract-references [^ReferenceSequenceFile f]
  "Extract all references into a seq of [name, bases]
   Bases are encoded in a byte array."
  (.reset f)
  (->> (repeatedly #(.nextSequence f))
       (take-while identity)
       (map (fn [^ReferenceSequence s]
              [(.getName s) (.getBases s)]))))
