(in-package #:cl-bayesian-asd)

(defpackage #:cl-bayesian
  (:use :common-lisp
        :iterate
        :bind
        :cl-utilities
        :anaphora
        :xarray
        :lla)
  (:shadowing-import-from :iterate :collecting :collect))
