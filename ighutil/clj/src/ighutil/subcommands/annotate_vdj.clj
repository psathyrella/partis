(ns ighutil.subcommands.annotate-vdj
  (:require [cliopatra.command :refer [defcommand]]
            [clojure.data.csv :as csv]
            [clojure.string :as string]
            [plumbing.core :refer [map-vals]]
            [ighutil.gff3 :as gff3]
            [ighutil.imgt :as imgt]
            [ighutil.io :as zio]
            [ighutil.sam :as sam]
            [lonocloud.synthread :as ->])
  (:import [net.sf.samtools AlignmentBlock SAMRecord]
           [net.sf.picard.util IntervalTreeMap]
           [io.github.cmccoy.sam SAMUtils AlignedPair
            AlignedPair$MatchesReference]
           [io.github.cmccoy.dna Codons]
           [io.github.cmccoy.sam VDJAnnotator
            VDJAnnotator$RegionAnnotation]
           [java.util Arrays]))

(defn- translate-read
  "Translate read sequence"
  [^SAMRecord read frame start]
  (-> read
      .getReadBases
      (Codons/translateSequence (+ (int frame) (int start)))
      (String.)))

(defn- aligned-boundaries [^SAMRecord read]
  (let [ap (.getAlignmentBlocks read)
        ^AlignmentBlock f (first ap)
        ^AlignmentBlock l (last ap)
        start (.getReadStart f)
        end (+ (.getReadStart l) (.getLength l))]
    [start end]))

(defn- annotate-read [^SAMRecord read ^IntervalTreeMap tree]
  (let [annot (.annotateVDJ (VDJAnnotator. read tree))
        to-map (fn [^VDJAnnotator$RegionAnnotation r]
                 {:reference (.name r)
                  :aligned (.aligned r)
                  :eroded5p (.eroded5P r)
                  :eroded3p (.eroded3P r)
                  :mismatch (.mismatch r)
                  :qstart (.qstart r)
                  :qend (inc (.qend r))
                  :alignment-score (.alignmentScore r)
                  :ties (.getAttribute read sam/TAG-TIES)
                  :nm (.nm r)})]
    (map-vals to-map annot)))

(defn- identify-insertions [reads]
  (let [boundaries (into
                    {}
                    (for [^SAMRecord r reads]
                      [(sam/ig-segment r)
                       (aligned-boundaries r)]))
        read-bases (.getReadBases ^SAMRecord (first reads))
        v-end (get-in boundaries [\V 1])
        d-start (get-in boundaries [\D 0])
        d-end (get-in boundaries [\D 1])
        j-start (get-in boundaries [\J 0])
        nnil? (complement nil?)]
    (-> {}
        (->/when (and (nnil? v-end) (nnil? d-start) (>= d-start v-end))
          (assoc :vd
            (String. (Arrays/copyOfRange
                      read-bases
                      (int v-end)
                      (int d-start)))))
        (->/when (and (nnil? d-end) (nnil? j-start) (>= j-start d-end))
          (assoc :dj
            (String. (Arrays/copyOfRange
                      read-bases
                      (int d-end)
                      (int j-start))))))))

(defn- annotate-reads
  "Annotate read using (possibly) multiple alignment segments"
  [reads ^IntervalTreeMap tree & {:keys [full-translation]}]
  (let [^SAMRecord f (first reads)
        annotations (->> reads
                         (map #(annotate-read % tree))
                         (apply merge))
        v-start (get-in annotations ["V" :qstart])
        cys-start (get-in annotations ["Cys" :qstart])
        tryp-start (get-in annotations ["J_Tryp" :qstart])
        tryp-end (get-in annotations ["J_Tryp" :qend])
        in-frame (when (not (some nil? [cys-start tryp-start]))
                   (= (mod (int cys-start) 3) (mod (int tryp-start) 3)))
        frame (when (not (some nil? [cys-start tryp-end]))
                (mod (- (int cys-start) (int v-start)) 3))
        cdr3-length (when (not (some nil? [cys-start tryp-end]))
                      (inc (- (int tryp-end) (int cys-start))))
        translation (when in-frame
                      (if full-translation
                        (translate-read f (mod cys-start 3) 0)
                        (translate-read f frame v-start)))]
    (merge
     {:name (sam/read-name f)
      :cdr3-length cdr3-length
      :frame frame
      :read-bases (String. (.getReadBases f))
      :has-stop (when (not (nil? translation))
                  (.contains ^String translation "*"))
      :translation translation
      :in-frame in-frame
      :annotations annotations}
     (identify-insertions reads))))

(defn- features-for-references [reference-names
                                ^IntervalTreeMap tm]
  (let [reference-names (into #{} reference-names)]
    (->> (gff3/all-feature-names tm :key identity)
         (filter (comp reference-names first))
         (map second)
         (into #{}))))

(defn- keys-for-feature [feature-name]
  (let [k [:aligned :mismatch :qstart :qend]]
    (if (#{"V" "D" "J"} feature-name)
      (conj k :reference :alignment-score :nm :ties :eroded5p :eroded3p)
      k)))

(defn- names-for-feature [feature-name]
  (->> feature-name
       keys-for-feature
       (map #(str feature-name "_" (string/replace (name %) \- \_)))))

(defn- header-row
  "Generate a header row for a collection of feature names"
  [feature-names]
  (->> feature-names
       (map names-for-feature)
       (apply concat ["read_name" "cdr3_length" "frame" "in_frame"
                      "has_stop" "read-bases" "translation" "vd" "dj"])))

(defn- to-row
  "Convert the results of annotate-reads to a row for CSV output"
  [{:keys [name cdr3-length in-frame frame
           has-stop annotations read-bases translation vd dj] :as m}
   feature-names]
  (->> feature-names
       (map (fn [x] (for [k (keys-for-feature x)]
                      (get-in annotations [x k]))))
       (apply concat [name cdr3-length frame in-frame has-stop read-bases translation vd dj])))

(defcommand annotate-vdj
  "Annotate VDJ alignments"
  {:opts-spec [["-i" "--in-file" "Source file" :required true
                :parse-fn sam/bam-reader]
               ["-o" "--out-file" "Destination path" :required true
                :parse-fn zio/writer]
               ["-g" "--gff3" "GFF3 file [default: internal]"
                :parse-fn gff3/slurp-gff3]
               ["--[no-]full-translation"
                "Completely translate sequence (outside of optimal alignment)?"
                :default false]]}
  (with-open [^net.sf.samtools.SAMFileReader in-file in-file
              ^java.io.Writer out-file out-file]
    (let [reference-names (sam/reference-names in-file)
          gff3 (or gff3 @imgt/ighvj-gff)
          feature-tree (gff3/gff3-to-interval-map
                        gff3
                        :key (comp :Name :attributes))
          feature-names (vec
                         (concat
                          (features-for-references
                           reference-names
                           feature-tree)
                          ["V" "D" "J"]))]
      (csv/write-csv out-file [(header-row feature-names)])
      (->> in-file
           .iterator
           iterator-seq
           sam/partition-by-name
           (map #(annotate-reads % feature-tree
                                 :full-translation full-translation))
           (map #(to-row % feature-names))
           (csv/write-csv out-file)))))
