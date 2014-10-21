(ns ighutil.subcommands.reset-primary
  "Resets the primary record for each sequence,
  first adding sequence data to each segment (IGHV, IGHD, IGHJ),
  then finding alleles which may be present in the data,
  finally choosing randomly among best-scoring alignments."
  (:require [clojure.java.io :as io]
            [clojure.set :as set]
            [clojure.string :as string]
            [cliopatra.command :refer [defcommand]]
            [ighutil.sam :refer [primary?
                                 alignment-score
                                 read-name
                                 reference-name
                                 partition-by-name
                                 partition-by-name-segment
                                 ig-segment
                                 bam-writer]
             :as sam]
            [ighutil.imgt :refer [strip-allele]]
            [me.raynes.fs :as fs]
            [plumbing.core :refer [frequencies-fast]]
            [taoensso.timbre :as timbre])
  (:import [net.sf.samtools
            SAMRecord
            SAMFileReader
            SAMFileReader$ValidationStringency
            SAMFileHeader$SortOrder]))

(defn- make-primary!
  "Make a record primary"
  [^SAMRecord read
   ^bytes read-bases
   ^bytes read-qualities &
   {:keys [supplementary] :or {supplementary false}}]
  (doto read
    (.setReadBases read-bases)
    (.setBaseQualities read-qualities)
    (.setNotPrimaryAlignmentFlag false)
    (.setSupplementaryAlignmentFlag supplementary)))

(defn- make-secondary!
  "Make a record secondary"
  [^SAMRecord read]
  (doto read
    (.setNotPrimaryAlignmentFlag true)
    (.setReadBases SAMRecord/NULL_SEQUENCE)
    (.setBaseQualities SAMRecord/NULL_QUALS)
    (.setSupplementaryAlignmentFlag false)))

(defn- set-supp-and-bases-per-gene
  "Set read bases and quals for each gene primary record"
  [reads]
  (when (some primary? reads)
    (assert (= 1 (count (filter primary? reads))))
    (let [^SAMRecord primary (first (filter primary? reads))
          ^bytes read-seq (.getReadBases primary)
          ^bytes read-qual (.getBaseQualities primary)
          update-first (fn [r]
                         (when (->> r (filter primary?) count (= 0))
                           (let [^SAMRecord f (first r)]
                             (make-primary!
                              f
                              read-seq
                              read-qual
                              :supplementary (not= (ig-segment f) \V))))
                         r)]
      (->> reads
           (partition-by ig-segment)
           (mapcat update-first)))))

(defn- unlikely-alleles
  "List alleles present at less than 'min-prop' proportion of reads for a gene"
  [alignments & {:keys [min-prop] :or {min-prop 0.1}}]
  (letfn [(filter-max-score
            [alignments]
            (let [sorted (sort-by alignment-score #(compare %2 %1) alignments)
                  max-score (-> sorted first alignment-score)
                  to-keep (vec (take-while #(= max-score (alignment-score %))
                                           sorted))
                  n-keep (count to-keep)]
              (into
               {}
               (for [record to-keep]
                 [(reference-name record) (/ 1.0 n-keep)]))))
          (drop-group
            [freqs]
            (let [all (->> freqs (map first) (into #{}))
                  tot (reduce + (map second freqs))
                  minimum (* min-prop (float tot))
                  to-keep (->> freqs
                               (filter (comp #(>= % minimum) second))
                               (sort-by (comp vec reverse) #(compare %2 %1))
                               (map first)
                               (take 2))]
              (set/difference all to-keep)))]
    (->> alignments
         partition-by-name-segment
         (map filter-max-score)
         (reduce #(merge-with + %1 %2))
         seq
         sort
         (partition-by (comp strip-allele first))
         (map drop-group)
         (reduce set/union))))


(defn- max-score-tiebreak
  "alignments should all have the same query

   Break the tie between SAM records with identical maximum alignment
   score using rand-nth.

   Drops any alignments to *to-remove*"
  [alignments & {:keys [record-selector to-remove]
                 :or {record-selector rand-nth
                      to-remove #{}}}]
  (assert (set? to-remove))
  (let [is-unlikely-allele? (comp to-remove reference-name)]
    (if (= 1 (count alignments))
     (remove is-unlikely-allele? alignments) ; Only one alignment
     (let [^SAMRecord primary (first (filter primary? alignments))
           sorted (sort-by alignment-score #(compare %2 %1) alignments)
           dropped-unlikely (remove is-unlikely-allele? sorted)]
       (assert (= 1 (count (filter primary? alignments))))
       (if (seq dropped-unlikely)
         (let [max-score (-> dropped-unlikely
                             first
                             alignment-score)
               max-records (vec (take-while
                                 (comp #(= max-score %) alignment-score)
                                 dropped-unlikely))
            ^SAMRecord selection (record-selector max-records)]
           (when (and primary (not= selection primary))
             (do
               ;; Swap out the primary read and the new selection
            (make-primary!
             selection
             (.getReadBases primary)
             (.getBaseQualities primary)
             :supplementary (.getSupplementaryAlignmentFlag primary))
            (make-secondary! primary)))
           ;; Add number of ties
           (let [ties (->> max-records
                           (remove #(= selection %))
                           (map reference-name)
                           sort
                           (string/join "|"))]
             (.setAttribute selection sam/TAG-TIES ties))
           dropped-unlikely)
         ;; drop all alignments for the read
         [])))))

(defn- order-by-vdj
  "Place reads in V-D-J order"
  [reads]
  (let [segment-order {\V 0 \D 1 \J 2}
        k (juxt (comp segment-order ig-segment) (complement primary?))]
    (sort-by k reads)))

(defn- set-frame-for-group
  "Set frame based on the first alignment"
  [reads]
  (let [reads (vec reads)
        ^SAMRecord v-alignment (first reads)
        frame (-> v-alignment
                  sam/position
                  (mod 3)
                  int)]
    (assert (primary? v-alignment))
    (when (= \V (sam/ig-segment v-alignment))
      (doseq [^SAMRecord r reads]
        (.setAttribute r sam/TAG-FRAME frame))))
  reads)

(defn reset-primary-record
  "Given a sequence of SAMRecord objects, sorted by name, assigns a
   new primary record for each partition (if :randomize is true
   [default]), copies bases and qualities to each primary segment
   record, sets the supplementary alignment flag for each non-V segment"
  [sam-records & {:keys [randomize to-remove random-gen]
                  :or {randomize true
                       to-remove #{}}}]
  (when  (and randomize (not random-gen))
    (timbre/error "Missing random number generator")
    (throw (IllegalArgumentException. "Missing random number generator")))
  (let [sel (if randomize
              (fn [reads]
                (let [reads (vec reads)
                      n (count reads)]
                  (nth reads (.nextInt ^java.util.Random random-gen n))))
              first)]
    (->> sam-records
         partition-by-name
         (mapcat set-supp-and-bases-per-gene)
         partition-by-name-segment
         (map vec)
         (mapcat #(max-score-tiebreak
                   %
                   :record-selector sel
                   :to-remove to-remove))
         partition-by-name
         (mapv order-by-vdj)
         (mapcat set-frame-for-group))))

(defcommand reset-primary
  "Reset primary alignment, breaking ties randomly"
  {:opts-spec [["-i" "--in-file"
                "Source BAM - must be sorted by *name*"
                :required true :parse-fn io/file]
               ["-o" "--out-file" "Destination path"
                :required true
                :parse-fn io/file]
               ["--seed" "Random number seed"
                :parse-fn #(Long/parseLong %)
                :default 0]
               ["--[no-]filter-alleles"
                "[Do not] filter to a maximum of 2 alleles per gene"
                :default true]
               ["--[no-]randomize"
                "[do not] randomize primary record among top-scoring alignments."
                :default true]
               ["--[no-]compress" "Compress output?" :default true]]}
  (assert (not= in-file out-file) "Same infile and outfile")
  (assert (fs/exists? in-file) "Non-extant input file.")
  (let [to-remove (if filter-alleles
                    (with-open [reader (SAMFileReader. ^java.io.File in-file)]
                      (->> reader
                           .iterator
                           iterator-seq
                           unlikely-alleles))
                    #{})
        random-gen (java.util.Random. ^Long seed)]
    (timbre/info "Removing " (sort to-remove))
    (with-open [reader (SAMFileReader. ^java.io.File in-file)]
      (.setValidationStringency
       reader
       SAMFileReader$ValidationStringency/SILENT)
      (let [header (.getFileHeader reader)
            read-iterator (->> reader
                               .iterator
                               iterator-seq)]
        (.setSortOrder header SAMFileHeader$SortOrder/unsorted)
        (timbre/info "Filtering / resetting to"
                     (.getName ^java.io.File out-file))
        (with-open [writer (bam-writer
                            out-file
                            reader
                            :compress compress)]
          (doseq [^SAMRecord read (reset-primary-record
                                   read-iterator
                                   :randomize randomize
                                   :random-gen random-gen
                                   :to-remove to-remove)]
            (assert (not (nil? read)))
            (.addAlignment writer read)))))))
