;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian-tests)

(deftestsuite cl-bayesian-tests () ()
  (:equality-test #'==))

(defun small-zs? (z-statistics &key (threshold 2) (margin 4))
  "Test if the number of Z statistics larger than THRESHOLD is smaller than
the expected number when corrected by MARGIN."
  (format t "z-statistics: ~A~%" z-statistics)
  (let* ((count (count-if (curry #'< threshold) z-statistics))
         (p (* 2 (- 1 (cdf (r-normal) threshold))))
         (allowed-large-z (ceiling (* margin p (length z-statistics)))))
    (<= count allowed-large-z)))

;; EXTERNAL

(defun run ()
  "Run all the tests."
  (run-tests :suite 'cl-bayesian-tests))
