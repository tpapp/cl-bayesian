(defpackage #:cl-bayesian-tests
    (:use :common-lisp :iterate :bind :anaphora :tpapp-utils :cl-num-utils
          :alexandria :lla :cl-random :cl-bayesian :lift)
  (:shadowing-import-from :iterate :collecting :collect)
  (:shadowing-import-from :cl-random :variance) ; also in alexandria
  (:shadowing-import-from :cl-num-utils :mean :xor) ; also in alexandria
  (:export run-cl-bayesian-tests))
