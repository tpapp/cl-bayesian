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
    ((:file "chains"))))
  :depends-on
  (iterate metabang-bind anaphora alexandria lift cl-num-utils lla cl-random
   cl-random-tests))
