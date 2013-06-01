(ns ighutil.mutationwalker
  "Simple walker"
  (:import [org.broadinstitute.sting.gatk CommandLineGATK]
           [org.broadinstitute.sting.gatk.contexts AlignmentContext ReferenceContext]
           [org.broadinstitute.sting.utils.pileup ReadBackedPileup])
  (:gen-class
   :name io.github.cmccoy.mutationwalk.MutationWalker
   :extends io.github.cmccoy.mutationwalk.BasePositionWalker))

(defn -map
  "Count!"
  [this tracker ^ReferenceContext ref ^AlignmentContext context]
  (let [contig (.. context (getLocation) (getContig))
        loc (.. context (getLocation) (getStart))
        ref-base (char (.getBase ref))
        global-pileup (.getBasePileup context)
        samples (.getSamples global-pileup)
        sample-pileups (into {} (.getPileupsForSamples global-pileup samples))
        base-count (fn [^ReadBackedPileup pileup] (vec (.getBaseCounts pileup)))]
    (into {} (map
              (fn [[k v]] [[contig loc ref-base k] (base-count v)])
              sample-pileups))))

(defn -reduceInit
  "Initialize to empty"
  [this]
  (java.util.HashMap.))

(defn -reduce
  ""
  [this cur ^java.util.HashMap other]
  (doseq [[k v] (seq cur)]
    (.put other k v))
  other)

(defn -onTraversalDone
  ""
  [this result]
  (let [^java.io.PrintStream out (.out this)]
    (println result)))

(defn -main [& args]
  (let [args (into-array
              String
              (concat ["-T" "MutationWalker" "-l" "ERROR" "-dt" "NONE"] args))]
    (CommandLineGATK/start (CommandLineGATK.) args)))
