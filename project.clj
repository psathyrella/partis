(defproject ighutil "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojars.chapmanb/picard "1.90"]
                 [org.clojars.chapmanb/sam "1.90"]
                 [org.clojure/data.csv "0.1.2"]
                 [me.raynes/fs "1.4.3"]
                 [org.apache.commons/commons-compress "1.5"]
                 [com.google.guava/guava "14.0.1"]
                 [org.clojars.runa/cliopatra "1.1.0"]]
  :profiles {:1.6 {:dependencies [[org.codehaus.jsr166-mirror/jsr166y "1.7.0"]]}}
  :java-source-paths ["src/java"]
  :javac-options ["-target" "1.6" "-source" "1.6"]
  :omit-source false
  :global-vars {*warn-on-reflection* true}
  :aot [ighutil.main]
  :main ighutil.main)
