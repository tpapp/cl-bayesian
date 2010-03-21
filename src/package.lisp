(in-package #:cl-bayesian-asd)

(defpackage #:cl-bayesian
    (:use :common-lisp
          :iterate
          :bind
          :cl-utilities
          :anaphora
          :xarray
          :lla
          :tpapp-utils)
  (:shadowing-import-from :iterate :collecting :collect)
  (:export
   
   ;; mcmc

   acceptance-counter acceptance-ratio increment-counter mcmc reset-counters
   update update-parameter update-parameter-in-vector current-parameters define-mcmc
   define-updater define-metropolis-updater log-posterior-ratio metropolis-step*
   metropolis-step run-mcmc

   ))
