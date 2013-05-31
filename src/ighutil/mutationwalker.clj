(ns ighutil.mutationwalker
  "Simple walker"
  (:import [org.broadinstitute.sting.gatk CommandLineGATK])
  (:gen-class
   :name io.github.cmccoy.mutationwalk.MutationWalker
   :extends io.github.cmccoy.mutationwalk.BasePositionWalker))

(defn -map
  "Count!"
  [this tracker ref context]
  (if-not (nil? tracker) 1))

(defn -reduceInit
  "Initialize to 0"
  [this]
  0)

(defn -reduce
  ""
  [this cur other]
  (+ cur other))

(defn -onTraversalDone
  ""
  [this result]
  (let [^java.io.PrintStream out (.out this)]
    (.println out result)))

(defn -main [& args]
  (let [args (into-array
              String
              (concat ["-T" "MutationWalker"] args))]
    (CommandLineGATK/start (CommandLineGATK.) args)))
