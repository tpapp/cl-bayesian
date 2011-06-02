;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; -*-

(in-package #:cl-bayesian)

(defgeneric overdisperse (distribution factor &key &allow-other-keys)
  (:documentation "Return an overdispersed distribution.  Usually variance is
  blown up by FACTOR and the mean is preserved, but methods may behave
  differently.  The semantics is only define heuristically: use this for
  generating overdispersed distributions for MCMC.")
  (:method (distribution factor &key nu)
           (let+ (((&accessors-r/o mean variance) distribution)
                  (variance (e* factor variance)))
             (if (numberp mean)
                 (let ((sd (sqrt variance)))
                   (if nu
                       (r-t mean (/ sd (t-scale-to-variance-coefficient nu))
                            nu)
                       (r-normal mean sd)))
                 (if nu
                     (r-multivariate-t 
                      mean
                      (e/ variance (t-scale-to-variance-coefficient nu))
                      nu)
                     (r-normal mean variance)))))
  (:method ((distribution r-gamma) factor &key)
    (let+ (((&accessors-r/o alpha beta) distribution))
      (r-gamma (/ alpha factor) (* beta factor))))
  (:method ((distribution r-inverse-gamma) factor &key)
    (let+ (((&accessors-r/o alpha beta) distribution))
      (r-inverse-gamma (1+ (* factor (1- alpha))) (* factor beta))))
  (:method ((distribution r-multivariate-t) factor &key)
    (let+ (((&accessors-r/o multivariate-normal scaling-factor)
            distribution))
      (r-multivariate-t nil nil nil
                        :multivariate-normal multivariate-normal
                        :scaling-factor (overdisperse
                                         scaling-factor factor)))))
