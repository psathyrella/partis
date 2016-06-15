(require '[clojure.xml :as xml])

(require '[clojure.string :as string])

(def pth "/home/cmccoy/git/ighutil/python/vdjalign/imgt/data/ighj_adaptive.xml")

(def xdoc (xml/parse pth))

(defn create-gff-records [{:keys [immunoseq1 length cdr3index]}]
  (let [tryp-start (- (+ length cdr3index) 2)
        tryp-end (+ tryp-start 2)
        base {:seqid immunoseq1 :source "Adaptive"}]
    (mapv #(merge base %)
          [{:id (str "J_Tryp-" immunoseq1)
            :name "J_Tryp"
            :start tryp-start
            :end tryp-end}
           {:id (str "PreTryp-" immunoseq1)
            :name "PreTryp"
            :start 1
            :end (dec tryp-end)}
           {:id (str "PostTryp-" immunoseq1)
            :name "PostTryp"
            :start (inc tryp-end)
            :end length}])))

(defn load-xml [p]
  (let [acc (fn [acc i] (assoc acc i (Integer/parseInt (get acc i))))
        convert-types (fn [i] (reduce acc i [:cdr3index :length]))]
    (->> p
         xml/parse
         :content
         (map :attrs)
         (mapv convert-types))))

(defn create-gff-row [{:keys [source id name start end seqid]}]
  (format "%s\t%s\tregion\t%d\t%d\t.\t+\t0\tID=%s;Name=%s"
          seqid source start end id name))

(defn write-gff3 []
  (->> pth
       load-xml
       (mapcat create-gff-records)
       (map create-gff-row)
       (cons "##gff-version 3")
       (string/join "\n")
       (spit "ighj_adaptive.gff3")))
