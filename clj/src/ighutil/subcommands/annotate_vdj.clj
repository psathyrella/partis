(ns ighutil.subcommands.annotate-vdj
  (:require [cliopatra.command :refer [defcommand]]
            [clojure.data.csv :as csv]
            [plumbing.core :refer [?>]]
            [ighutil.gff3 :as gff3]
            [ighutil.imgt :as imgt]
            [ighutil.io :as zio]
            [ighutil.sam :as sam])
  (:import [net.sf.samtools SAMRecord]
           [net.sf.picard.util IntervalTree]
           [io.github.cmccoy.sam SAMUtils AlignedPair
            AlignedPair$MatchesReference]))

(defn- annotate-read [^SAMRecord read ^IntervalTree tree]
  (let [read-name (sam/read-name read)
        reference-name (sam/reference-name read)
        aligned-pairs (remove (fn [^AlignedPair ap]
                                (or (.isUnknown ap)
                                    (.isIndel ap)))
                              (SAMUtils/getAlignedPairs read))
        f (fn [m ^AlignedPair ap]
            (let [rpos (.getReferencePosition ap)
                  qpos (.getQueryPosition ap)
                  is-mismatch (.isMutation ap)
                  overlaps (conj (when tree (gff3/overlapping-vals
                                             tree
                                             (inc rpos)))
                                 reference-name)
                  add-base (fn [acc o]
                             (assoc acc o
                                    (-> (get acc o {:aligned 0 :mismatch 0
                                                    :type (if (= o reference-name)
                                                            "reference"
                                                            "feature")
                                                    :minqpos qpos
                                                    :maxqpos qpos})
                                        (update-in [:aligned] inc)
                                        (update-in [:maxqpos] #(min qpos %))
                                        (update-in [:maxqpos] #(max qpos %))
                                        (?> is-mismatch update-in [:mismatch] inc))))]
              (reduce add-base m overlaps)))]
    {:name read-name
     :reference reference-name
     :cdr3 (.getAttribute read "XL")
     :alignment-score (.getAttribute read "AS")
     :nm (.getAttribute read "NM")
     :annotations (reduce f {} aligned-pairs)} ))

(defn- build-rows [{:keys [name annotations cdr3 alignment-score nm reference]}]
  (for [[ann-name {:keys [aligned mismatch type minqpos maxqpos]}] annotations]
    [name cdr3 alignment-score nm reference ann-name type aligned mismatch
     minqpos maxqpos]))

(defcommand annotate-vdj
  "Annotate VDJ alignments"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn sam/bam-reader]
               ["-o" "--out-file" "Destination path" :required true
                :parse-fn zio/writer]]}
  (with-open [^net.sf.samtools.SAMFileReader in-file in-file
              ^java.io.Writer out-file out-file]
    (let [feature-tree (gff3/gff3-to-interval-map @imgt/ighvj-gff)]
      (csv/write-csv out-file [["read_name" "cdr3_length" "alignment_score"
                                "nm" "reference_name" "feature" "feature_type"
                                "n_aligned" "n_mutated" "minqpos" "maxqpos"]])
      (->> in-file
           .iterator
           iterator-seq
           (map (fn [^SAMRecord r]
                  (annotate-read r
                                 (->> r
                                      sam/reference-name
                                      (get feature-tree)))))
           (mapcat build-rows)
           (csv/write-csv out-file)))))
