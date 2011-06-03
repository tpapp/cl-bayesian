(defpackage cl-bayesian (:nicknames mcmc)
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
   update-chain state sample-chain reset-chain parameters scalar-parameters 
   mcmc-sample burn-in 
   
   ;; slice-sample

   slice-sample-so

   ;; samplers
   
   variance-distribution lr-kv-dummies lr-kv multivariate-normal-model kappa
   inverse-scale
   
   ;; chains

   psrf psrf-r psrf-v psrf-w psrf-ranges sample-column-statistics
   autocovariance-accumulators partial-ranges partial-accumulators
   column-statistics

   ))
