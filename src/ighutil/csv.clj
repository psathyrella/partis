(ns ighutil.csv
  (:require [clojure.data.csv :as csv]
            [plumbing.core :refer [?>>]]))

(defn int-of-string [^String s]
  (Integer/parseInt s))

(defn float-of-string [^String s]
  (Float/parseFloat s))

(defn double-of-string [^String s]
  (Double/parseDouble s))

(defn csv-to-maps [fp & {:keys [keywordize?]
                         :or {keywordize? true}
                         :as options}]
  (let [rows (apply csv/read-csv
                    fp
                    (dissoc options :keywordize?))
        header (->> rows
                    first
                    (?>> keywordize? map keyword))
        make-row (partial zipmap header)]
    (map make-row (rest rows))))

(defn read-typed-csv [fp types & csv-opts]
  (let [rows (apply csv-to-maps fp csv-opts)
        type-entry (fn [[k v]] [k ((get types k identity) v)])
        type-row (fn [m] (->> m (map type-entry) (into {})))]
    (map type-row rows)))
