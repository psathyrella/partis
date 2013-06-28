(ns ighutil.main
  (:require [cliopatra.command :as command]
            [ighutil.subcommands
             add-quality-scores
             calculate-match-probability
             count-mismatches
             enumerate-mutations
             identify-motif
             identify-subsets
             mutations-by-site
             reset-primary])
  (:gen-class
   :name ighutil.main))

(defn -main [& args]
  (let [commands ['ighutil.subcommands.enumerate-mutations
                  'ighutil.subcommands.identify-subsets
                  'ighutil.subcommands.count-mismatches
                  'ighutil.subcommands.identify-motif
                  'ighutil.subcommands.mutations-by-site
                  'ighutil.subcommands.add-quality-scores
                  'ighutil.subcommands.reset-primary
                  'ighutil.subcommands.calculate-match-probability]]
    (command/dispatch (sort commands) args)))
