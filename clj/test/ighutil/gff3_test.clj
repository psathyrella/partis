(ns ighutil.gff3-test
  (:require [ighutil.gff3 :refer :all]
            [clojure.test :refer :all]))

(deftest record-parsing
  (let [s "ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN
"]
    (is (= {:seqid "ctg123"
            :source nil
            :start 1000
            :end 9000
            :attributes {:ID "gene00001" :Name "EDEN"}}
           (parse-gff3-record s)))))
