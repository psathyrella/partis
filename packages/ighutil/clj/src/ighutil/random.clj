(ns ighutil.random
  (:import [java.util Random]))

(defprotocol RandomGenerator
  (random-int [this] [this ^Integer n]))

(defrecord RNG [^Random random]
  RandomGenerator
  (random-int [_] (.nextInt random))
  (random-int [_ n] (.nextInt random n)))

(defn rng
  ([] (RNG. (Random.)))
  ([^Long seed] (RNG. (Random. seed))))
