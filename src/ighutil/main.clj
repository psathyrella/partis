(ns ighutil.main
  (:require [cliopatra.command :as command]
            [ighutil enumerate-mutations
             identify-subsets count-mismatches
             mutations-by-site])
  (:gen-class
   :name ighutil.main))

(defn -main [& args]
  (command/dispatch
   ['ighutil.enumerate-mutations
    'ighutil.identify-subsets
    'ighutil.count-mismatches
    'ighutil.mutations-by-site]
   args))
