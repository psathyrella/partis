(ns ighutil.main
  (:require [cliopatra.command :as command]
            [ighutil.enumerate-mutations]
            [ighutil.mutationwalker])
  (:gen-class
   :name ighutil.main))

(defn -main [& args]
  (command/dispatch
   ['ighutil.enumerate-mutations
    'ighutil.mutationwalker]
   args))
