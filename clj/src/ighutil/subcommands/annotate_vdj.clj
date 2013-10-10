(ns ighutil.subcommands.annotate-vdj
  (:require [cliopatra.command :refer [defcommand]]
            [plumbing.core :refer [?>]]
            [ighutil.gff3 :as gff3]
            [ighutil.io :as zio]
            [ighutil.sam :as sam]
            [ighutil.sam-tags :refer [TAG-EXP-MATCH]])
  (:import [net.sf.samtools SAMRecord]
           [net.sf.picard.util IntervalTree]
           [io.github.cmccoy.sam SAMUtils AlignedPair
            AlignedPair$MatchesReference]))

(def ^:private MISMATCH AlignedPair$MatchesReference/FALSE)
(def ^:private UNKNOWN AlignedPair$MatchesReference/FALSE)

(defn- annotate-read [^SAMRecord read ^IntervalTree tree]
  (let [read-name (sam/read-name read)
        reference-name (sam/reference-name read)
        aligned-pairs (remove (fn [^AlignedPair ap]
                                (or (= UNKNOWN (.getMatchesReference ap))
                                    (.isIndel ap)))
                              (SAMUtils/getAlignedPairs read))
        f (fn [m ^AlignedPair ap]
            (let [rpos (.getReferencePosition ap)
                  is-mismatch (= MISMATCH (.getMatchesReference ap))
                  overlaps (conj (when tree (gff3/overlapping-vals (inc rpos)))
                                 reference-name)
                  inc-nil (fnil inc 0)
                  add-base (fn [acc o]
                             (-> acc
                                 (get-in [o])
                                 (update-in [:aligned] inc-nil)
                                 (?> is-mismatch update-in [:mismatch] inc)))]
              (reduce add-base m overlaps)))]
    (reduce f {} aligned-pairs)))

(defcommand annotate-vdj
  "Annotate VDJ alignments"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn sam/bam-reader]
               ["-o" "--out-file" "Destination path" :required true
                :parse-fn zio/writer]
               ["-g" "--gff" "GFF3 path" :required true
                :parse-fn zio/reader]]}
  (with-open [^net.sf.samtools.SAMFileReader in-file in-file]
    (let [feature-tree (with-open [^java.io.Reader gff gff]
                         (-> gff
                             line-seq
                             gff3/parse-gff3
                             gff3/gff3-to-interval-map))]
      (->> in-file
           .iterator
           iterator-seq
           (map (fn [^SAMRecord r]
                  (annotate-read r
                                 (->> r
                                      sam/reference-name
                                      (get feature-tree)))))))))
