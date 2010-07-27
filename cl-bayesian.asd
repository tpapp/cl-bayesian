(defpackage #:cl-bayesian-asd
  (:use :cl :asdf))

(in-package #:cl-bayesian-asd)

(defsystem #:cl-bayesian
  :description ""
  :author "Tamas K Papp"
  :license "LLGPL"
  :serial t 
  :components
  ((:module 
    "package-init"
    :pathname #P "src/"
    :components
    ((:file "package")))
   (:module
    "basics"
    :pathname #P"src/"
    :depends-on ("package-init")
    :serial t
    :components
    ((:file "utilities")
     (:file "mcmc")
     (:file "slice-sampling")
     (:file "samplers")
     ;; (:file "polynomials")
     )))
  :depends-on
  (:iterate :metabang-bind :anaphora :lla :cl-random :tpapp-utils :cl-num-utils
            :alexandria))
