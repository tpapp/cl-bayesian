(defpackage cl-bayesian
  (:nicknames mcmc)
  (:use common-lisp iterate let-plus anaphora cl-num-utils alexandria
        lla cl-random lla)
  (:shadowing-import-from cl-num-utils
                          mean variance xor) ; also in alexandria
  (:export
   
   ;; utilities - nothing is exported
   
   overdisperse

   ;; layout

   extract layout-length parse-layout flatten-into named-layout
   scalar-layout 

   ;; mcmc

   initialize-chain parameters-layout scalar-parameters-layout model
   common-model update-chain state sample-chain reset-chain parameters
   scalar-parameters  mcmc-sample burn-in 
   
   ;; slice-sample

   slice-sample-so

   ;; samplers
   
   variance-distribution lr-kv-dummies lr-kv multivariate-normal-model kappa
   inverse-scale
   
   ;; chains

   calculate-psrf psrf-r psrf-v psrf-w calculate-psrf-ranges
   
   mcmc-statistics accumulators autocovariance-accumulators sse-ranges
   sse-accumulators
   
   mcmc-summary  psrf accumulators mean-autocorrelations psrf-ranges
   summarize-mcmc-statistics pool-samples

   ))
