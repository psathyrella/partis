(ns ighutil.mutation-by-read
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [cliopatra.command :as command :refer [defcommand]]))

(defn identify-mutations-in-read [read reference-bases])

(defcommand enumerate-mutations
  "Enumerate mutations by sample"
  {:opts-spec [] :bind-args-to args}
  (println args))
