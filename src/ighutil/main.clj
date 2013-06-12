(ns ighutil.main
  (:require [cliopatra.command :as command]
            [ighutil.enumerate-mutations])
  (:gen-class
   :name ighutil.main))

(defn -main [& args]
  (command/dispatch
   ['ighutil.enumerate-mutations]
   args))
