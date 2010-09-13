(in-package #:cl-bayesian-asd)

(defpackage #:cl-bayesian
    (:use :common-lisp :iterate :bind :anaphora :tpapp-utils :cl-num-utils
          :alexandria :lla :cl-random)
  (:shadowing-import-from :iterate :collecting :collect)
  (:shadowing-import-from :cl-random :variance)
  (:shadowing-import-from :cl-num-utils :mean :xor) ; also in alexandria
  (:export
   
   ;; utilities - nothing is exported

   ;; mcmc

   acceptance-counter acceptance-ratio increment-counter mcmc reset-counters
   update update-parameter current-parameters parameters-ix 
   define-current-parameters define-mcmc define-updater
   define-metropolis-updater log-posterior-ratio metropolis-step*
   metropolis-step run-mcmc *stop-mcmc*

   ;; slice-sample

   slice-sample-so

   ;; samplers
   
   variance-distribution

   ;; chains

   mcmc-chains mcmc-class initargs parameters-ix chains burn-in 
   pooled-parameters run-mcmc-chains psrf chains-psrf
   calculate-pooled-parameters

   ))
