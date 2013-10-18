(ns ighutil.subcommands.annotate-vdj
  (:require [cliopatra.command :refer [defcommand]]
            [clojure.data.csv :as csv]
            [plumbing.core :refer [?>]]
            [ighutil.gff3 :as gff3]
            [ighutil.imgt :as imgt]
            [ighutil.io :as zio]
            [ighutil.sam :as sam])
  (:import [net.sf.samtools SAMRecord]
           [net.sf.picard.util IntervalTreeMap]
           [io.github.cmccoy.sam SAMUtils AlignedPair
            AlignedPair$MatchesReference]))

(defn- annotate-read [^SAMRecord read ^IntervalTreeMap tree]
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
                  overlaps (conj (gff3/overlapping
                                  tree
                                  reference-name
                                  (inc rpos))
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
                                        (update-in [:minqpos] #(min qpos %))
                                        (update-in [:maxqpos] #(max qpos %))
                                        (?> is-mismatch update-in [:mismatch] inc))))]
              (reduce add-base m overlaps)))]
    {:name read-name
     :reference reference-name
     :alignment-score (.getAttribute read "AS")
     :nm (.getAttribute read "NM")
     :annotations (reduce f {} aligned-pairs)}))

(defn- annotate-reads [reads ^IntervalTreeMap tree]
  "Annotate read using (possibly) multiple alignment segments"
  (let [^SAMRecord f (first reads)]
    {:name (sam/read-name f)}))



(defn- features-for-references [reference-names
                                ^IntervalTreeMap tm]
  (let [reference-names (into #{} reference-names)]
    (->> tm
         gff3/all-feature-names
         (filter (comp reference-names first))
         (map second)
         (into #{}))))

(def ^:private feature-keys [:type :reference :aligned :mutated
                             :minqpos :maxqpos])

(defn- names-for-feature [feature-name]
  (map #(str feature-name "_" (name %)) feature-keys))

(defn- header-row [feature-names]
  (->> feature-names
       (map names-for-feature)
       (apply concat ["read_name" "cdr3_length" "alignment_score" "nm"])))

(defn- to-row [{:keys [name cdr3-length alignment-score nm annotations]}
               feature-names]
  (->> feature-names
       (map #(for [k feature-keys]
               (get-in annotations [% k])))
       (apply concat [name cdr3-length alignment-score nm])))

(defcommand annotate-vdj
  "Annotate VDJ alignments"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn sam/bam-reader]
               ["-o" "--out-file" "Destination path" :required true
                :parse-fn zio/writer]]}
  (with-open [^net.sf.samtools.SAMFileReader in-file in-file
              ^java.io.Writer out-file out-file]
    (let [reference-names (sam/reference-names in-file)
          feature-tree (gff3/gff3-to-interval-map @imgt/ighvj-gff)
          feature-names (features-for-references reference-names feature-tree)]
      (csv/write-csv out-file [(header-row feature-names)])
      (->> in-file
           .iterator
           iterator-seq
           sam/partition-by-name
           (map #(annotate-read % feature-tree))
           (map #(to-row % feature-names))
           (csv/write-csv out-file)))))
