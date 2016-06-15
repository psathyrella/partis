(ns ighutil.main
  "Main entry point - dispatches to subcommands"
  (:require [cliopatra.command :as command])
  (:gen-class
   :name ighutil.main))

(defn- subcommand-name-to-symbol [^String s]
  (->> s
       (str "ighutil.subcommands.")
       symbol))

(def command-names (sort ["annotate-vdj"
                          "aligned-base-composition"
                          "calculate-match-probability"
                          "count-alignment-ties"
                          "count-mismatches"
                          "distinct-by-vdjcdr3"
                          "enumerate-mutations"
                          "identify-motif"
                          "kmer-matrix"
                          "mutation-correlation"
                          "mutations-by-site"
                          "reset-primary"]))

(defn dispatch-to-subcommand [& args]
  (let [command-syms (map subcommand-name-to-symbol command-names)]
    (apply require command-syms)
    (command/dispatch command-syms args)))

(defn -main [& args]
  (try
    (apply dispatch-to-subcommand args)
    (finally (shutdown-agents))))
