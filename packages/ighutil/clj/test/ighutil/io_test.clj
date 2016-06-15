(ns ighutil.io-test
  (:require [clojure.test :refer :all]
            [me.raynes.fs :as fs]
            [schema.test]
            [ighutil.io :refer :all]))

(use-fixtures :once schema.test/validate-schemas)

(def ^:dynamic *tmpdir* nil)

(defn temp-dir-fixture [f]
  (binding [*tmpdir* (fs/temp-dir "ighutil-test")]
    (try
      (f)
      (finally
        (fs/delete-dir *tmpdir*)))))

(use-fixtures :each temp-dir-fixture)

(deftest test-input-streams
  (are [pth cls]
    (= (with-open [istream (input-stream pth)]
         (-> istream class .getSimpleName))
       cls)
    "testdata/test.txt" "BufferedInputStream"
    "testdata/test.txt.gz" "GzipCompressorInputStream"
    "testdata/test.txt.bz2" "BZip2CompressorInputStream"))

(deftest test-output-streams
  (are [pth cls]
    (= (with-open [ostream (->> pth (fs/file *tmpdir*) output-stream)]
         (-> ostream class .getSimpleName))
       cls)
    "test.txt" "BufferedOutputStream"
    "test.txt.gz" "GzipCompressorOutputStream"
    "test.txt.bz2" "BZip2CompressorOutputStream"))

(deftest test-readers
  (are [pth]
    (=
     (with-open [istream (reader pth)] (slurp istream))
     "Test!\n")
    "testdata/test.txt"
    "testdata/test.txt.gz"
    "testdata/test.txt.bz2"))

(deftest test-writers
  (are [pth cls]
    (= (with-open [ostream (->> pth (fs/file *tmpdir*) writer)]
         (-> ostream class .getSimpleName))
       cls)
    "test.txt" "BufferedWriter"
    "test.txt.gz" "BufferedWriter"
    "test.txt.bz2" "BufferedWriter"))
