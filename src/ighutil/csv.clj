(ns ighutil.csv
  (:require [clojure.data.csv :as csv]))

(defn csv-to-maps [fp & {:keys [keywordize]
                         :or {keywordize true}
                         :as options}]
  (let [rows (apply (partial csv/read-csv fp)
                    (dissoc options :keywordize))
        header ((if keywordize (partial map keyword) identity)
                (first rows))
        make-row (partial zipmap header)]
    (map make-row (rest rows))))

(defn read-typed-csv [fp types & csv-opts]
  (let [rows (apply (partial csv-to-maps fp) csv-opts)
        type-entry (fn [[k v]] [k ((get types k identity) v)])
        type-row (fn [m] (->> m (map type-entry) (into {})))]
    (map type-row rows)))
