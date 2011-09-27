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

   start-chain scalar-parameters-layout model
   common-model scalar-parameters draw-chain
   
   ;; slice-sample

   slice-sample-so

   ;; samplers
   
   univariate-normal-error univariate-normal-model lr-kv-dummies lr-kv
   multivariate-normal-model kappa inverse-scale
   
   ;; chains

   *suggested-minimum-burn-in**default-burn-in-fraction* discard-burn-in

   calculate-psrf psrf-r psrf-v psrf-w calculate-psrf-ranges
   
   mcmc-statistics accumulators autocovariance-accumulators sse-ranges
   sse-accumulators

   mcmc-summary psrf accumulators mean-autocorrelations psrf-ranges
   summarize-mcmc-statistics pool-samples
   
   ;; validation

   calculate-empirical-ranks calculate-p-statistics calculate-abs-z-statistics

   ;; dlm

   dlm-evolution1 dlm-evolution1-G dlm-evolution1-mu dlm-evolution1-W

   dlm-observation1 dlm-observation1-F dlm-observation1-V

   uddu uddu-u uddu-d uddu-update uddu-multiply-update

   ;; dlm-evolution dlm-evolution-G dlm-evolution-mu dlm-evolution-W

   ;; dlm-observation dlm-observation-F dlm-observation-V

   dlm-forward-filtering dlm-backward-sampling dlm-ff-bs
   dlm-errors dlm-simulate

   ))
