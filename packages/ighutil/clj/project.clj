(defproject ighutil "0.1.0-SNAPSHOT"
  :description "Tools for working with immunoglobulin sequences at high throughput"
  :license {:name "GPL v3"
            :url "https://www.gnu.org/licenses/gpl.txt"}
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [org.clojars.chapmanb/picard "1.107"]
                 [org.clojars.chapmanb/sam "1.107"]
                 [org.clojure/data.csv "0.1.2"]
                 [org.clojure/math.combinatorics "0.0.7"]
                 [me.raynes/fs "1.4.5"]
                 [org.apache.commons/commons-compress "1.8"]
                 [com.google.guava/guava "16.0.1"]
                 [org.clojars.runa/cliopatra "1.1.0"]
                 [prismatic/hiphip "0.2.0"]
                 [prismatic/plumbing "0.2.2"]
                 [prismatic/schema "0.2.1"]
                 ;[org.flatland/useful "0.11.2"]
                 [primitive-math "0.1.3"]
                 [cheshire "5.3.1"]
                 [com.taoensso/timbre "3.1.6"]
                 [lonocloud/synthread "1.0.4"]]
  :profiles {:java6 {:dependencies [[org.codehaus.jsr166-mirror/jsr166y "1.7.0"]]}
             :uberjar {:aot [ighutil.fasta
                             ighutil.gff3
                             ighutil.imgt
                             ighutil.io
                             ighutil.main
                             ighutil.sam
                             ighutil.subcommands.reset-primary]}}
  :java-source-paths ["src/java"]
  :omit-source false
  :global-vars {*warn-on-reflection* true}
  :main ighutil.main)
