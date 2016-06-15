(ns ighutil.gff3-test
  (:require [ighutil.gff3 :refer :all]
            [clojure.test :refer :all]
            [schema.test]))

(use-fixtures :once schema.test/validate-schemas)

(deftest record-parsing
  (let [s "ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN
"]
    (is (= {:seqid "ctg123"
            :source nil
            :start 1000
            :end 9000
            :attributes {:ID "gene00001" :Name "EDEN"}}
           (parse-gff3-record s)))))

(deftest interval-tree-map
  (let [records [{:seqid "ctg123" :source nil :start 1000 :end 9000 :attributes {:ID "gene00001"
                                                                                 :Name "Coding!"}}
                 {:seqid "ctg123" :source nil :start 10 :end 51 :attributes {:ID "feat00002"
                                                                             :Name "Feature!"}}]
        tm (gff3-to-interval-map records)]
    (is (= [] (overlapping tm "ctg123" 9)))
    (is (= [(second records)] (overlapping tm "ctg123" 51)))
    (is (= [(first records)] (overlapping tm "ctg123" 2500)))
    (is (= ["Feature!" "Coding!"] (map second
                                       (all-feature-names tm))))))
