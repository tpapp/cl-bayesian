(defsystem cl-bayesian
  :description ""
  :author "Tamas K Papp"
  :license "LLGPL"
  :in-order-to ((test-op (test-op cl-bayesian-tests)))
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
     (:file "chains")
     (:file "validation")
     (:file "dlm")
     ;; (:file "polynomials")
     )))
  :depends-on
  (iterate let-plus anaphora alexandria cl-num-utils lla cl-random))

(defsystem cl-bayesian-tests
  :description "Unit tests for the CL-BAYESIAN library."
  :author "Tamas K Papp"
  :license "Same as CL-BAYESIAN -- this is part of the latter."
  :serial t 
  :components
  ((:module 
    "package-init"
    :pathname #P "tests/"
    :components
    ((:file "package")))
   (:module
    "setup"
    :pathname #P"tests/"
    :depends-on ("package-init")
    :serial t
    :components
    ((:file "setup")))
   (:module
    "tests"
    :pathname #P"tests/"
    :components
    ((:file "chains")
     (:file "dlm"))))
  :depends-on
  (iterate let-plus anaphora alexandria lift cl-num-utils lla cl-random
   cl-random-tests))
