;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

(defun variance-distribution (ss n prior)
  "When residuals ~ NIID(0,variance), return a posterior distribution for the
variance with the given prior.  SS is the sum of squared residuals, and N is
the number of observations.  Possible priors are:

  :REFERENCE  -- p(variance) \propto 1/variance

  :HIERARCHICAL -- p(variance) \propto (variance)^-1/2, recommended by
  Gelman (2006) for hierarchical models with at least 3 groups, as a first
  attempt.

  :NONE -- no prior, just gives the likelihood (may not be proper)."
  (let* ((alpha (+ (/ n 2d0)
                   (ecase prior
                     (:none -1d0)
                     (:hierarchical -0.5d0)
                     (:reference 0d0)))))
    (r-inverse-gamma alpha (/ ss 2d0))))

;;; linear regression with known variance
;;; 
;;; Not used frequently in practice, but useful for Gibbs sampling.  Return a
;;; multivariate normal, which is the posterior of the coefficients.

(defun lr-kv-dummies (prior)
  "Return dummy observations as (Y . X) for the given prior, for use in a
linear regression with known variance (LR-KV)."
  (check-type prior r-multivariate-normal)
  (bind (((:accessors-r/o mean variance-left-sqrt) prior))
    (cons (solve variance-left-sqrt mean)
          (invert variance-left-sqrt))))

(defun lr-kv (y x variance &key prior)
  "Linear regression of Y on X with known VARIANCE for the errors (a single
scalar is accepted, in which case it is used as a diagonal matrix).  Use
LR-KV-DUMMIES to generate dummy observations from a prior."
  (bind ((x (as-regression-covariates x))
         ((:values y-transformed x-transformed) (transform-y-x y x variance))
         ((:values y-transformed x-transformed)
          (add-regression-dummies y-transformed x-transformed prior
                                  #'lr-kv-dummies))
         ((:values beta nil nil qr)
          (least-squares y-transformed x-transformed :method :qr)))
    (r-multivariate-normal beta (invert-xx qr))))
