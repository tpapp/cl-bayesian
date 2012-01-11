(defpackage #:cl-bayesian
  (:nicknames #:mcmc)
  (:use #:common-lisp #:iterate #:let-plus #:anaphora #:cl-num-utils
        #:alexandria #:lla #:cl-random)
  (:shadowing-import-from #:cl-num-utils  ; also in alexandria
                          #:mean #:variance #:xor #:median)
  (:export
   ;; utilities
   #:overdisperse
   ;; layout
   #:extract
   #:layout-length
   #:parse-layout
   #:flatten-into
   #:named-layout
   #:scalar-layout
   ;; mcmc
   #:start-chain
   #:scalar-parameters-layout
   #:model
   #:common-model
   #:scalar-parameters
   #:draw-chain
   ;; slice-sample
   #:slice-sample-so
   ;; samplers
   #:univariate-normal-error
   #:univariate-normal-model
   #:lr-kv-dummies
   #:lr-kv
   #:multivariate-normal-model
   #:kappa
   #:inverse-scale
   ;; chains
   #:*suggested-minimum-burn-in**default-burn-in-fraction*
   #:discard-burn-in
   #:calculate-psrf
   #:psrf-r
   #:psrf-v
   #:psrf-w
   #:calculate-psrf-ranges
   #:mcmc-statistics
   #:accumulators
   #:autocovariance-accumulators
   #:sse-ranges
   #:sse-accumulators
   #:mcmc-summary
   #:psrf
   #:accumulators
   #:mean-autocorrelations
   #:psrf-ranges
   #:summarize-mcmc-statistics
   #:pool-samples
   ;; validation
   #:calculate-empirical-ranks
   #:calculate-p-statistics
   #:calculate-abs-z-statistics
   ;; dlm
   #:dlm-evolution1
   #:dlm-evolution1-G
   #:dlm-evolution1-mu
   #:dlm-evolution1-W
   #:dlm-observation1
   #:dlm-observation1-F
   #:dlm-observation1-V
   #:dlm
   #:make-dlm
   #:dlm-length
   #:sub-dlm
   #:uddu
   #:uddu-u
   #:uddu-d
   #:uddu-update
   #:uddu-multiply-update
   ;; dlm-evolution dlm-evolution-G dlm-evolution-mu dlm-evolution-W
   ;; dlm-observation dlm-observation-F dlm-observation-V
   #:dlm-forward-filtering
   #:dlm-backward-sampling
   #:dlm-ff-bs
   #:dlm-errors
   #:dlm-simulate
   #:pack-slots
   #:unpack-slots
   #:mosaic
   #:make-mosaic
   #:template-mosaic
   #:template-mosaic-symbols
   #:make-mosaic-matrix
   #:mcmc-mosaic-matrix
   #:mosaic-vector
   #:make-mosaic-vector))
