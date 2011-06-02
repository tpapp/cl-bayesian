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
  (let+ (((&accessors-r/o mean variance-left-sqrt) prior))
    (cons (solve variance-left-sqrt mean)
          (invert variance-left-sqrt))))

(defun lr-kv (y x variance &key prior)
  "Linear regression of Y on X with known VARIANCE for the errors (a single
scalar is accepted, in which case it is used as a diagonal matrix).  Use
LR-KV-DUMMIES to generate dummy observations from a prior."
  (let+ ((x (as-regression-covariates x))
         ((&values y-transformed x-transformed) (transform-y-x y x variance))
         ((&values y-transformed x-transformed)
          (add-regression-dummies y-transformed x-transformed prior
                                  #'lr-kv-dummies))
         ((&values beta nil nil qr)
          (least-squares y-transformed x-transformed :method :qr)))
    (r-multivariate-normal beta (invert-xx qr))))

;;; multivariate normal model
;;;
;;; 

(defclass multivariate-normal-model ()
  ((inverse-scale :accessor inverse-scale :initarg :inverse-scale)
   (nu :reader nu :initarg :nu :documentation "Degrees of freedom.")
   (kappa :reader kappa :initarg :kappa :documentation "Number of
   observations (including dummies from prior, may be a fraction).")
   (mean :reader mean :initarg :mean :documentation "Posterior mean."))
  (:documentation "Random variable representing the posterior for a
  multivariate normal distribution estimated with unknown variance, reference
  or conjugate prior.  Second values return Sigma, the variance matrix."))

(defun multivariate-normal-model (y &key prior)
  "Estimate a multivariate normal model.  See p85-88 of Bayesian Data
Analysis, 2nd edition.  If prior is not given, it is the reference prior."
  (let+ (((n nil) (array-dimensions y))
         ((&values sse mean) (matrix-sse y))
         (kappa n)
         (nu n))
    (when prior
      (check-type prior multivariate-normal-model)
      (let ((kappa0 (kappa prior))
            (mean0 (mean prior)))
        (incf kappa kappa0)
        (incf nu (nu prior))
        (setf sse (e+ sse
                      (inverse-scale prior)
                      (mm (e- mean mean0) t
                          (/ (* kappa0 n) kappa))))
        (setf mean (e+ (e* (mean prior) (/ kappa0 kappa))
                       (e* mean (/ n kappa))))))
    (make-instance 'multivariate-normal-model :inverse-scale sse
                   :nu nu :kappa kappa :mean mean)))

(defmethod draw ((multivariate-normal-model multivariate-normal-model) &key
                 as-distribution?)
  (let+ (((&slots-r/o inverse-scale nu kappa mean) multivariate-normal-model)
         (sigma (draw (r-inverse-wishart nu inverse-scale)))
         (mean (draw (r-multivariate-normal mean sigma) :scale (/ kappa))))
    (if as-distribution?
        (r-multivariate-normal mean sigma)
        (values mean sigma))))
