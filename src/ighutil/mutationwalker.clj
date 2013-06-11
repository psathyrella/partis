(ns ighutil.mutationwalker
  "Simple walker"
  (:import [org.broadinstitute.sting.gatk CommandLineGATK]
           [org.broadinstitute.sting.gatk.contexts AlignmentContext ReferenceContext]
           [org.broadinstitute.sting.utils.pileup
            ReadBackedPileup
            PileupElement])
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [cliopatra.command :refer [defcommand]])
  (:gen-class
   :name io.github.cmccoy.mutationwalk.MutationWalker
   :extends io.github.cmccoy.mutationwalk.BasePositionWalker
   :implements [org.broadinstitute.sting.gatk.walkers.TreeReducible]))

(defn- count-of-name [name]
  "Extract count from reads labeled 'name_count'"
  (let [[_ count] (re-matches #".*_(\d+)$" name)]
    (Integer/parseInt (or count 1))))

(defn- extract-ops-counts [^PileupElement x]
  (let [count (count-of-name (.. x (getRead) (getReadName)))]
    (cons (if (.isDeletion x)
            [:del count]
            [(-> (.getBase x) char str keyword) count])
          (when (> (.getLengthOfImmediatelyFollowingIndel x) 0)
            [[:ins count]])) ))

(defn- process-pileup [^ReadBackedPileup pileup]
  (let [coverage (.depthOfCoverage pileup)
        ops-counts (mapcat extract-ops-counts
                           (iterator-seq (.iterator pileup)))]
    (persistent!
     (reduce
      (fn [m [op c]]
        (assoc! m op (+ (get m op 0) c)))
      (transient {:A 0 :C 0 :G 0 :ins 0 :del 0})
      ops-counts))))

(defn- split-pileup-by-sample [^ReadBackedPileup pileup]
  (let [samples (.getSamples pileup)]
    (into (hash-map) (.getPileupsForSamples pileup samples))))

(defn -map
  "Count!"
  [this tracker ^ReferenceContext ref ^AlignmentContext context]
  (let [contig (.. context (getLocation) (getContig))
        loc (.. context (getLocation) (getStart))
        ref-base (char (.getBase ref))
        pileups (split-pileup-by-sample (.getBasePileup context))
        global-pileup (.getBasePileup context)]
    (into
     (hash-map)
     (map
      (fn [[k v]] [[contig loc ref-base k] (process-pileup v)])
      pileups))))

(defn -reduceInit
  "Initialize to empty"
  [this]
  [])

(defn -reduce
  ""
  [this cur other]
  (conj other cur))

(defn -onTraversalDone
  "Print a report"
  [this result]
  (let [out (io/writer (.out this))
        m (apply (partial merge-with (partial merge-with +)) result)]
    (csv/write-csv
     out
     [["reference" "position" "wt" "sample" "A" "C" "G" "T" "Ins" "Del" "coverage"]])
    (doseq [[[contig loc ref-base sample] v] m]
      (let [[a c g t ins del cov] ((juxt :A :C :G :T :ins :del :coverage) v)
            row [contig loc ref-base sample a c g t ins del cov]]
        (csv/write-csv out [row])))))

(defn -treeReduce [lhs rhs]
  (concat lhs rhs))

(defcommand mutations-csv
  "Mutation rates to CSV file using GATK"
  {:opts-spec []
   :bind-args-to args}
  (let [args (into-array
              String
              (concat ["-T" "MutationWalker" "-dt" "NONE"] args))]
    (CommandLineGATK/start (CommandLineGATK.) args)))
