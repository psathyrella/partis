(ns ighutil.main
  (:require [cliopatra.command :as command]
            [ighutil
             enumerate-mutations
             identify-subsets
             count-mismatches
             mutations-by-site
             add-quality-scores
             identify-motif])
  (:gen-class
   :name ighutil.main))

(defn -main [& args]
  (let [commands ['ighutil.enumerate-mutations
                  'ighutil.identify-subsets
                  'ighutil.count-mismatches
                  'ighutil.identify-motif
                  'ighutil.mutations-by-site
                  'ighutil.add-quality-scores]]
    (command/dispatch (sort commands) args)))
