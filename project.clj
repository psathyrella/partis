(defproject ighutil "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 ;; GATK Requirements from bcbio.variation
                 ;;[org.clojars.chapmanb/gatk-lite "2.5.2"]
                 [org.clojars.chapmanb/picard "1.90"]
                 [org.clojars.chapmanb/sam "1.90"]
                 ;;[org.clojars.chapmanb/tribble "1.90"]
                 ;;[org.clojars.chapmanb/variant "1.90"]
                 ;;[org.simpleframework/simple-xml "2.0.4"]
                 ;;[colt/colt "1.2.0"]
                 ;;[it.unimi.dsi/fastutil "6.5.3"]
                 ;;[log4j/log4j "1.2.17"]
                 ;;[commons-lang/commons-lang "2.5"]
                 ;;[org.apache.servicemix.bundles/org.apache.servicemix.bundles.jets3t "0.8.1_1"]
                 ;; End GATK
                 [org.clojure/data.csv "0.1.2"]
                 [me.raynes/fs "1.4.3"]
                 [org.apache.commons/commons-compress "1.5"]
                 [com.google.guava/guava "14.0.1"]
                 [org.clojars.runa/cliopatra "1.1.0"]
                 [org.codehaus.jsr166-mirror/jsr166 "1.7.0"]]
  :java-source-paths ["src/java"]
  :javac-options ["-target" "1.6" "-source" "1.6"]
  :omit-source false
  :global-vars {*warn-on-reflection* true}
  :aot [ighutil.main]
  :main ighutil.main)
