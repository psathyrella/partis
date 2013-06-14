(ns ighutil.main
  (:require [cliopatra.command :as command]
            [ighutil enumerate-mutations
             identify-subsets count-mismatches])
  (:gen-class
   :name ighutil.main))

(defn -main [& args]
  (command/dispatch
   ['ighutil.enumerate-mutations
    'ighutil.identify-subsets
    'ighutil.count-mismatches]
   args))
