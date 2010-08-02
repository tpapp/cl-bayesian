;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian-tests)

(deftestsuite diagnostics-tests (cl-bayesian-tests) ())

(addtest (diagnostics-tests)
  (bind ((n 100)
         (zeroes (make-array n :initial-element 0d0))
         (ones (make-array n :initial-element 1d0))
         (twos (make-array n :initial-element 2d0))
         (s1 (concat zeroes ones))
         (s2 (concat zeroes twos))
         ((:values r r-upper) (potential-scale-reduction (list s1 s2))))
    (ensure-same r 1.284474)
    ;; (ensure-same r-upper 2.190318)
    ))
