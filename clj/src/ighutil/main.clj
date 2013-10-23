(ns ighutil.main
  (:require [cliopatra.command :as command])
  (:gen-class
   :name ighutil.main))

(defn- subcommand-name-to-symbol [^String s]
  (->> s
       (str "ighutil.subcommands.")
       symbol))

(def command-names (sort ["annotate-vdj"
                          "calculate-match-probability"
                          "count-alignment-ties"
                          "count-mismatches"
                          "enumerate-mutations"
                          "identify-motif"
                          "identify-subsets"
                          "kmer-matrix"
                          "mutation-correlation"
                          "aligned-base-composition"
                          "mutations-by-site"
                          "reset-primary"]))

(defn -main [& args]
  (let [command-syms (map subcommand-name-to-symbol command-names)]
    (apply require command-syms)
    (command/dispatch command-syms args)))
