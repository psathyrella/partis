(ns ighutil.io
  (:require [me.raynes.fs :as fs]
            [clojure.java.io :as io])
  (:import [org.apache.commons.compress.compressors.bzip2
            BZip2CompressorOutputStream BZip2CompressorInputStream]
           [org.apache.commons.compress.compressors.gzip
            GzipCompressorOutputStream GzipCompressorInputStream]))

(defn- output-stream-for-name [file-name]
  (condp = (fs/extension file-name)
    ".bz2" #(BZip2CompressorOutputStream. %)
    ".gz"  #(GzipCompressorOutputStream. %)
    identity))

(defn- input-stream-for-name [file-name]
  (condp = (fs/extension file-name)
    ".bz2" #(BZip2CompressorInputStream. %)
    ".gz"  #(GzipCompressorInputStream. %)
    identity))

(defn ^java.io.OutputStream output-stream [file-name]
  (if (= file-name "-")
    (io/output-stream System/out)
  (let [cls (output-stream-for-name file-name)]
    (-> file-name io/file io/output-stream cls))))

(defn ^java.io.Writer writer [file-name]
  (-> file-name output-stream io/writer))

(defn ^java.io.InputStream input-stream [file-name]
  (if (= file-name "-")
    (io/input-stream System/in)
    (let [cls (input-stream-for-name file-name)]
      (-> file-name io/file io/input-stream cls))))

(defn ^java.io.Reader reader [file-name]
  (-> file-name input-stream io/reader))
