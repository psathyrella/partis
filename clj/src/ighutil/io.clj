(ns ighutil.io
  "Mostly mirrors clojure.java.io, but supporting gzip- and bzip2-compressed
  inputs and outputs for files ending in .gz and .bz2 respectively, and
  stdin/stdout for '-'."
  (:require [clojure.java.io :as io]
            [me.raynes.fs :as fs]
            [schema.core :as s])
  (:import [org.apache.commons.compress.compressors.bzip2
            BZip2CompressorOutputStream BZip2CompressorInputStream]
           [org.apache.commons.compress.compressors.gzip
            GzipCompressorOutputStream GzipCompressorInputStream]))

(s/defschema Fileable
  "types supported by clojure.java.io/file"
  (s/either s/Str java.io.File java.net.URL))

(defn- output-stream-for-name [file-name]
  (case (fs/extension file-name)
    ".bz2" #(BZip2CompressorOutputStream. %)
    ".gz" #(GzipCompressorOutputStream. %)
    identity))

(defn- input-stream-for-name [file-name]
  (case (fs/extension file-name)
    ".bz2" #(BZip2CompressorInputStream. %)
    ".gz"  #(GzipCompressorInputStream. %)
    identity))

(s/defn output-stream :- java.io.OutputStream
  "Create an output stream for a file. '-' is mapped to stdout"
  [file-name :- Fileable]
  (if (= file-name "-")
    (io/output-stream System/out)
    (let [cls (output-stream-for-name file-name)]
      (-> file-name io/file io/output-stream cls))))

(s/defn writer :- java.io.Writer
  "Create a writer for a file. '-' is mapped to stdout"
  [file-name :- Fileable]
  (-> file-name output-stream io/writer))

(s/defn input-stream :- java.io.InputStream
  "Create an input stream for a file. '-' is mapped to stdin"
  [file-name :- Fileable]
  (if (= file-name "-")
    (io/input-stream System/in)
    (let [cls (input-stream-for-name file-name)]
      (-> file-name io/file io/input-stream cls))))

(s/defn reader :- java.io.Reader
  "Create a reader for a file. '-' is mapped to stdin"
  [file-name :- Fileable]
  (-> file-name input-stream io/reader))
