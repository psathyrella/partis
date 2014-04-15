(ns ighutil.csv
  (:require [clojure.data.csv :as csv]
            [plumbing.core :refer [?>>]]
            [schema.core :as s]
            [schema.macros :as sm]))

(sm/defn int-of-string :- s/Int [s :- s/Str]
  (Integer/parseInt s))

(sm/defn float-of-string :- s/Num [s :- s/Str]
  (Float/parseFloat s))

(sm/defn double-of-string :- s/Num [s :- s/Str]
  (Double/parseDouble s))

(defn csv-to-maps
  "Read records from a CSV file, creating a map for each row"
  [fp & {:keys [keywordize?]
         :or {keywordize? true}
         :as options}]
  (let [rows (apply csv/read-csv
                    fp
                    (apply concat (dissoc options :keywordize?)))
        header (->> rows
                    first
                    (?>> keywordize? map keyword))
        make-row #(zipmap header %)]
    (map make-row (rest rows))))

(sm/defn read-typed-csv :- s/Any
  "Read a CSV file, converting types"
  [fp types :- {s/Keyword s/Any} & csv-opts]
  (let [rows (apply csv-to-maps fp csv-opts)
        type-entry (fn [m [k v]] (assoc! m k (v (get m k))))
        type-row (fn [row] (persistent!
                            (reduce
                             type-entry
                             (transient row)
                             types)))]
    (map type-row rows)))
