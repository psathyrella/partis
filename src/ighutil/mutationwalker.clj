(ns ighutil.mutationwalker
  "Simple walker"
  (:import [org.broadinstitute.sting.gatk CommandLineGATK]
           [org.broadinstitute.sting.gatk.contexts AlignmentContext ReferenceContext]
           [org.broadinstitute.sting.utils.pileup ReadBackedPileup])
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io])
  (:gen-class
   :name io.github.cmccoy.mutationwalk.MutationWalker
   :extends io.github.cmccoy.mutationwalk.BasePositionWalker
   :implements [org.broadinstitute.sting.gatk.walkers.TreeReducible]))

(defn- process-pileup [^ReadBackedPileup pileup]
  (let [[a c g t] (vec (.getBaseCounts pileup))
        del (.getNumberOfDeletions pileup)
        ins (.getNumberOfInsertionsAfterThisElement pileup)]
    {:A a :C c :G g :T t :ins ins :del del}))

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
  [this])

(defn -reduce
  ""
  [this cur other]
  (conj other cur))

(defn -onTraversalDone
  ""
  [this result]
  (let [out (io/writer (.out this))
        m (apply (partial merge-with +) result)]
    (csv/write-csv
     out
     [["contig" "position" "wt" "sample" "A" "C" "G" "T" "Ins" "Del"]])
    (doseq [[[contig loc ref-base sample] v] m]
      (let [[a c g t ins del] ((juxt :A :C :G :T :ins :del) v)
            row [contig loc ref-base sample a c g t ins del]]
        (csv/write-csv out [row])))))

(defn -treeReduce [lhs rhs]
  (concat lhs rhs))

(defn -main [& args]
  (let [args (into-array
              String
              (concat ["-T" "MutationWalker" "-dt" "NONE"] args))]
    (CommandLineGATK/start (CommandLineGATK.) args)))
