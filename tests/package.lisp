(defpackage cl-bayesian-tests
  (:use common-lisp iterate let-plus anaphora alexandria cl-num-utils lla
        cl-random cl-bayesian lift)
  ;; also in alexandria
  (:shadowing-import-from cl-num-utils mean variance xor)
  (:export run-cl-bayesian-tests))
