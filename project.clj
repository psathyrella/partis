(defproject ighutil "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 ;; GATK Requirements from bcbio.variation
                 [org.clojars.chapmanb/gatk-lite "2.5.2"]
                 [org.clojars.chapmanb/picard "1.90"]
                 [org.clojars.chapmanb/sam "1.90"]
                 [org.clojars.chapmanb/tribble "1.90"]
                 [org.clojars.chapmanb/variant "1.90"]
                 ;; [org.clojars.chapmanb/cofoja "1.0-r139"]
                 [org.clojars.chapmanb/jama "1.0.2"]
                 [org.apache.commons/commons-jexl "2.1.1"]
                 [org.apache.commons/commons-math "2.2"]
                 [org.reflections/reflections "0.9.8"]
                 [org.simpleframework/simple-xml "2.0.4"]
                 [colt/colt "1.2.0"]
                 [it.unimi.dsi/fastutil "6.5.3"]
                 [log4j/log4j "1.2.17"]
                 [commons-lang/commons-lang "2.5"]
                 [org.apache.servicemix.bundles/org.apache.servicemix.bundles.jets3t "0.8.1_1"]
                 ;; End GATK
                 [org.clojure/data.csv "0.1.2"]]
  :java-source-paths ["src/java"]
  :javac-options ["-target" "1.6" "-source" "1.6"]
  :omit-source false
  :global-vars {*warn-on-reflection* true}
  :aot [ighutil.mutationwalker]
  :main io.github.cmccoy.mutationwalk.MutationWalker)
