;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

(defun univariate-normal-error (accumulator &optional (prior 0))
  "Univariate normal error model.  Sum of squared errors and the tally are
taken from accumulator.  The following priors are accepted:

  - r-inverse-gamma distributions

  - real numbers, corresponding to an inverse-gamma(prior,0) distribution in
    the limit (which, of course, may not be proper).  Eg prior=0 is usually
    used as a reference prior."
  (let+ (((&values a b) (aetypecase prior
                            (r-inverse-gamma (values (alpha it) (beta it)))
                            (real (values it 0)))))
    (r-inverse-gamma (+ (/ (tally accumulator) 2d0) a)
                     (+ (/ (sse accumulator 0d0) 2d0) b))))

;;;; ****************************************************************
;;;; Normal-inverse-chi-square distribution.
;;;;
;;;; 
;;;;
;;;;
;;;; ****************************************************************

(cl-random::define-rv r-univariate-normal-model (mu kappa alpha beta)
  (:documentation "Posterior distributions for univariate normal model -- the
   draws themselves are normal distributions, with sigma^2 ~
   inverse-gamma(alpha,beta) and m ~ N(mu,s^2/kappa).  This rv is also a
   conjugate prior for univariate normal data.  See BDA, Section 3.3.  The
   reparametrization here is alpha=nu/2 and beta=nu*s^2/2."
   :==-slots (mu kappa sigma-square))
  ((mu :type double-float :reader t)
   (kappa :type double-float :reader t)
   (sigma-square :type r-inverse-gamma))
  (cl-random::with-doubles (mu kappa alpha beta)
    (cl-random::make :mu mu :kappa kappa
                     :sigma-square (r-inverse-gamma alpha beta)))
  (draw (&key)
        (let ((s^2 (draw sigma-square)))
          (r-normal (draw (r-normal mu (/ s^2 kappa))) s^2))))

(define-indirect-accessors r-univariate-normal-model 
    r-univariate-normal-model-sigma-square s^2 nu alpha beta)

(defun univariate-normal-model (accumulator &optional prior)
  (let+ (((&accessors-r/o (n tally) mean sse) accumulator)
         ((&values mu0 kappa0 alpha0 beta0)
          (aetypecase prior
            (r-univariate-normal-model
             (values (mu it) (kappa it) (alpha it) (beta it)))
            (null
             (values 0d0 0d0 -0.5d0 0d0))))
         (kappa (+ kappa0 n))
         (mu (/ (+ (* kappa0 mu0) (* n mean)) kappa))
         (alpha (+ alpha0 (/ n 2)))
         (beta (+ beta0 (/ (+ sse (/ (* kappa0 n (expt (- mean mu0) 2))
                                     kappa))
                           2))))
    (r-univariate-normal-model mu kappa alpha beta)))



;; (defun variance-distribution (ss n prior)
;;   "When residuals ~ NIID(0,variance), return a posterior distribution for the
;; variance with the given prior.  SS is the sum of squared residuals, and N is
;; the number of observations.  Possible priors are:

;;   :REFERENCE  -- p(variance) \propto 1/variance

;;   :HIERARCHICAL -- p(variance) \propto (variance)^-1/2, recommended by
;;   Gelman (2006) for hierarchical models with at least 3 groups, as a first
;;   attempt.

;;   :NONE -- no prior, just gives the likelihood (may not be proper)."
;;   (let* ((alpha (+ (/ n 2d0)
;;                    (ecase prior
;;                      (:none -1d0)
;;                      (:hierarchical -0.5d0)
;;                      (:reference 0d0)))))
;;     (r-inverse-gamma alpha (/ ss 2d0))))

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
   observations (including dummies from prior, does not have to be an
   integer).")
   (mean :reader mean :initarg :mean :documentation "Posterior mean."))
  (:documentation "Random variable representing the posterior for a
  multivariate normal distribution estimated with unknown variance, reference
  or conjugate prior.  Second value returned is Sigma, the variance matrix."))

(defun multivariate-normal-model (y &key prior)
  "Estimate a multivariate normal model.  See p85-88 of Bayesian Data
Analysis, 2nd edition.  If prior is not given, the reference prior is used."
  (declare (optimize debug))
  (let+ (((n nil) (array-dimensions y))
         ((&values sse mean) (matrix-sse y))
         ((&values mu0 kappa0 nu0 inverse-scale0)
          (aetypecase prior
            (multivariate-normal-model
             (values (mean it) (kappa it) (nu it) (inverse-scale it)))
            (null
             (values (zero) 0d0 -1 (zero)))))
         (kappa (+ kappa0 n))
         (nu (+ nu0 n))
         (mu (e+ (e* (/ kappa0 kappa) mu0) (e* (/ n kappa) mean)))
         (inverse-scale (e+ inverse-scale0
                            sse
                            (outer (e- mean mu0) t (/ (* kappa0 n) kappa)))))
    (make-instance 'multivariate-normal-model :inverse-scale inverse-scale
                                              :nu nu :kappa kappa :mean mu)))

(defmethod draw ((multivariate-normal-model multivariate-normal-model) &key)
  (let+ (((&slots-r/o inverse-scale nu kappa mean) multivariate-normal-model)
         (sigma (draw (r-inverse-wishart nu inverse-scale)))
         (mean (draw (r-multivariate-normal mean sigma) :scale (/ kappa))))
    (r-multivariate-normal mean sigma)))
