(ns ighutil.subcommands.annotate-vdj
  (:require [cliopatra.command :refer [defcommand]]
            [clojure.data.csv :as csv]
            [clojure.string :as string]
            [plumbing.core :refer [?>]]
            [ighutil.gff3 :as gff3]
            [ighutil.imgt :as imgt]
            [ighutil.io :as zio]
            [ighutil.sam :as sam])
  (:import [net.sf.samtools SAMRecord]
           [net.sf.picard.util IntervalTreeMap]
           [io.github.cmccoy.sam SAMUtils AlignedPair
            AlignedPair$MatchesReference]
           [io.github.cmccoy.dna Codons]))

(defn- translate-read [^SAMRecord read frame start]
  (-> read
      .getReadBases
      (Codons/translateSequence (+ (int frame) (int start)))
      (String.)))

(defn- annotate-read [^SAMRecord read ^IntervalTreeMap tree]
  (let [reference-name (sam/reference-name read)
        nm (sam/nm read)
        as (sam/alignment-score read)
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
                                 (sam/ig-segment read))
                  add-base (fn [acc o]
                             (assoc
                                 acc o
                                 (-> (get acc o {:aligned 0
                                                 :mismatch 0
                                                 :reference reference-name
                                                 :minqpos qpos
                                                 :maxqpos qpos
                                                 :alignment-score as
                                                 :nm nm})
                                     (update-in [:aligned] inc)
                                     (update-in [:minqpos] #(min qpos %))
                                     (update-in [:maxqpos] #(max qpos %))
                                     (?> is-mismatch update-in [:mismatch] inc))))]
              (reduce add-base m overlaps)))]
    (reduce f {} aligned-pairs)))

(defn- annotate-reads [reads ^IntervalTreeMap tree]
  "Annotate read using (possibly) multiple alignment segments"
  (let [^SAMRecord f (first reads)
        annotations (->> reads
                         (map #(annotate-read % tree))
                         (apply merge))
        v-start (get-in annotations [\V :minqpos])
        cys-start (get-in annotations ["Cys" :minqpos])
        tryp-end (get-in annotations ["J_Tryp" :maxqpos])
        frame (when (not (some nil? [cys-start tryp-end]))
                (mod (- (int cys-start) (int v-start)) 3))
        cdr3-length (when (not (some nil? [cys-start tryp-end]))
                      (inc (- tryp-end cys-start)))
        translation (when frame (translate-read f frame v-start))]
    {:name (sam/read-name f)
     :cdr3-length cdr3-length
     :frame frame
     :has-stop (when (not (nil? translation))
                 (.contains ^String translation "*"))
     :translation translation
     :in-frame (when cdr3-length (= 0 (mod (int cdr3-length) 3)))
     :annotations annotations}))

(defn- features-for-references [reference-names
                                ^IntervalTreeMap tm]
  (let [reference-names (into #{} reference-names)]
    (->> (gff3/all-feature-names tm :key identity)
         (filter (comp reference-names first))
         (map second)
         (into #{}))))

(defn- keys-for-feature [feature-name]
  (let [k [:aligned :mismatch :minqpos :maxqpos]]
    (if (#{\V \D \J} feature-name)
      (conj k :reference :alignment-score :nm)
      k)))

(defn- names-for-feature [feature-name]
  (->> feature-name
       keys-for-feature
       (map #(str feature-name "_" (string/replace (name %) \- \_)))))

(defn- header-row [feature-names]
  (->> feature-names
       (map names-for-feature)
       (apply concat ["read_name" "cdr3_length" "frame" "in_frame"
                      "has_stop" "translation"])))

(defn- to-row [{:keys [name cdr3-length in-frame frame
                       has-stop annotations translation] :as m}
               feature-names]
  (->> feature-names
       (map (fn [x] (for [k (keys-for-feature x)]
                      (get-in annotations [x k]))))
       (apply concat [name cdr3-length frame in-frame has-stop translation])))

(defcommand annotate-vdj
  "Annotate VDJ alignments"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn sam/bam-reader]
               ["-o" "--out-file" "Destination path" :required true
                :parse-fn zio/writer]
               ["-g" "--gff3" "GFF3 file [default: internal]"
                :parse-fn gff3/slurp-gff3
                :default @imgt/ighvj-gff]]}
  (with-open [^net.sf.samtools.SAMFileReader in-file in-file
              ^java.io.Writer out-file out-file]
    (let [reference-names (sam/reference-names in-file)
          feature-tree (gff3/gff3-to-interval-map
                        gff3
                        :key (comp :Name :attributes))
          feature-names (vec
                         (concat
                          (features-for-references reference-names feature-tree)
                          "VDJ"))]
      (csv/write-csv out-file [(header-row feature-names)])
      (->> in-file
           .iterator
           iterator-seq
           sam/partition-by-name
           (map #(annotate-reads % feature-tree))
           (map #(to-row % feature-names))
           (csv/write-csv out-file)))))
