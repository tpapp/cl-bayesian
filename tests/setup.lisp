;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian-tests)

(deftestsuite cl-bayesian-tests () ()
  (:equality-test #'==))

;; EXTERNAL

(defun run-cl-bayesian-tests ()
  "Run all the tests for LLA."
  (run-tests :suite 'cl-bayesian-tests))
