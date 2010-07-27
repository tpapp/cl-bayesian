;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

(defun variance-distribution (residuals prior)
  "When residuals ~ NIID(0,variance), return a posterior distribution for the
variance with the given prior.  Possible priors are: 

  :REFERENCE  -- p(variance) \propto 1/variance

  :HIERARCHICAL -- p(variance) \propto (variance)^-1/2, recommended by
  Gelman (2006) for hierarchical models with at least 3 groups, as a first
  attempt.

  :NONE -- no prior, just gives the likelihood (may not be proper)."
  (let* ((ss (dot t residuals))
         (alpha (+ (/ (length residuals) 2d0)
                   (ecase prior
                     (:none -1d0)
                     (:hierarchical -0.5d0)
                     (:reference 0d0)))))
    (make-instance 'inverse-gamma :alpha alpha :beta (/ ss 2d0))))

