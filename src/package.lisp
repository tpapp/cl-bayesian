(defpackage cl-bayesian
  (:nicknames mcmc)
  (:use common-lisp iterate let-plus anaphora cl-num-utils alexandria
        lla cl-random)
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
   
   univariate-normal-error univariate-normal-model lr-kv-dummies lr-kv
   multivariate-normal-model kappa inverse-scale
   
   ;; chains

   calculate-psrf psrf-r psrf-v psrf-w calculate-psrf-ranges
   
   mcmc-statistics accumulators autocovariance-accumulators sse-ranges
   sse-accumulators
-   
   mcmc-summary  psrf accumulators mean-autocorrelations psrf-ranges
   summarize-mcmc-statistics pool-samples
   
   ;; validation

   calculate-empirical-ranks calculate-p-statistics calculate-abs-z-statistics

   ;; dlm

   dlm1-parameters make-dlm1-parameters dlm1-parameters-G dlm1-parameters-mu
   dlm1-parameters-W dlm1-parameters-F dlm1-parameters-V dlm1-step dlm1-filter
   dlm1-sample dlm1-forward-filtering dlm1-backward-sampling dlm1-ff-bs
   dlm1-errors dlm1-simulate
   
   uddu uddu-u uddu-d uddu-update uddu-multiply-update dlm-parameters
   make-dlm-parameters dlm-parameters-G dlm-parameters-mu dlm-parameters-W
   dlm-parameters-F dlm-parameters-V dlm-forward-filtering
   dlm-backward-sampling dlm-ff-bs dlm-errors dlm-simulate

   ))
