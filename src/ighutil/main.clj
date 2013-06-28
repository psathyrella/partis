(ns ighutil.main
  (:require [cliopatra.command :as command]
            [ighutil.subcommands
             enumerate-mutations
             identify-subsets
             count-mismatches
             mutations-by-site
             add-quality-scores
             identify-motif
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
                  'ighutil.subcommands.reset-primary]]
    (command/dispatch (sort commands) args)))
