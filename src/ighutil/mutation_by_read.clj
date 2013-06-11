(ns ighutil.mutation-by-read
  (:import [net.sf.picard.reference
            ReferenceSequenceFileFactory
            ReferenceSequenceFile]
           [net.sf.samtools SAMRecord AlignmentBlock SAMFileReader])
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [clojure.core.memoize :refer (memo-lu)]
            [cliopatra.command :refer [defcommand]]))

;;; TODO: Better as for?
(defn- mutations-in-block [^bytes query-bases
                           ^bytes reference-bases
                           ^AlignmentBlock block]
  (let [l (.getLength block)]
    (loop [res []
           qi (dec (.getReadStart block))
           ri (dec (.getReferenceStart block))]
      (if (< qi l)
        (let [q (aget query-bases qi)
              r (aget reference-bases ri)
              eq (= q r)]
          (recur (if (not= q r)
                   (conj res {:ref-idx ri :qry-idx qi :ref (char r) :qry (char q)})
                   res)
                 (inc qi)
                 (inc ri)))
        res))))

(defn identify-mutations-in-read [^SAMRecord read reference-bases]
  (let [aligned-blocks (.getAlignmentBlocks read)
        s (.getReadBases read)]
    (into [] (mapcat
              (partial mutations-in-block s reference-bases)
              aligned-blocks))))

(defn- ref-getter [^ReferenceSequenceFile f]
  (memo-lu
   (fn [^String contig]
     (.getBases (.getSequence f contig)))))

(defn- identify-mutations-in-sam [^SAMFileReader sam ref-get]
  (let [sam-records (-> sam (.iterator) iterator-seq)]
    (for [^SAMRecord r (remove (fn [^SAMRecord r]
                                 (or (.getReadUnmappedFlag r)
                                     (.getNotPrimaryAlignmentFlag r)))
                               sam-records)]
      (let [ref-bases (ref-get (.getReferenceName r))]
        (identify-mutations-in-read r ref-bases)))))

(defcommand enumerate-mutations
  "Enumerate mutations by sample"
  {:opts-spec [["-j" "--jobs" "Number of processors" :default 1
                :parse-fn #(Integer. ^String %)]]
   :bind-args-to [reference sam-path]}
  (time (let [ref (-> reference
                 io/file
                 ReferenceSequenceFileFactory/getReferenceSequenceFile)
         ref-get (ref-getter ref)]
     (with-open [sam (SAMFileReader. (io/file sam-path))]
       (let [muts (identify-mutations-in-sam sam ref-get)]
         (doall muts))))))

(defn test-enumeration []
  (enumerate-mutations ["ighv.fasta" "AJ_memory_001_v.bam"]))
