;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package :cl-bayesian-tests)

(deftestsuite samplers-tests (cl-bayesian-tests)
  ())

(defun 2phase-posterior (n accumulator posterior-function)
  "Return two posteriors, the first accumulated in 2 sweeps, the first in
one."
  (let* ((v (filled-array (* 2 n) (generator (r-normal))))
         (acc1 (sweep (funcall accumulator) (subseq v 0 n)))
         (acc2 (sweep (funcall accumulator) (subseq v n)))
         (acc (sweep (funcall accumulator) v))
         (posterior1 (funcall posterior-function acc1))
         (posterior2 (funcall posterior-function acc2 posterior1)))
    (values posterior2 (funcall posterior-function acc))))

(defun simulate-ranks (prior simulate estimate parameters 
                       &key (n 500) (m 50))
  "Draw a value from PRIOR, call (SIMULATE DRAW) to simulate data,
then (ESTIMATE DATA PRIOR) to estimate the posterior.  Obtain N draws from the
latter, and return empirical ranks.  PARAMETERS is used to convert
prior/posterior draws to vectors."
  (filled-array m
                (lambda ()
                  (let+ ((theta0 (draw prior))
                         (y (funcall simulate theta0))
                         (posterior (funcall estimate y prior))
                         (theta+ (filled-array n (generator posterior))))
                    (calculate-empirical-ranks (funcall parameters theta0)
                                               (combine
                                                (map1 parameters theta+)))))))

(addtest (samplers-tests)
  univariate-normal-error-2phase
  (let+ (((&values p12 p) (2phase-posterior 10 #'mean-sse-accumulator
                                            #'univariate-normal-error))
         (ranks+
          (simulate-ranks (r-inverse-chi-square 9 2d0)
                          (lambda (v)
                            (filled-array 50 (generator (r-normal 0 v))))
                          (lambda (y prior)
                            (univariate-normal-error 
                             (sweep (mean-sse-accumulator) y)
                             prior))
                          #'vector))
         (z-statistics (calculate-abs-z-statistics ranks+)))
    (ensure-same p12 p :test #'==)
    (format t "z-statistics: ~A~%" z-statistics)
    (ensure (< (first* z-statistics) 3d0))))

(addtest (samplers-tests)
  univariate-normal-model-2phase
  (let+ (((&values p12 p) (2phase-posterior 10 #'mean-sse-accumulator
                                            #'univariate-normal-model))
         (prior (univariate-normal-model
                 (sweep (mean-sse-accumulator) (ivec 10))))
         (ranks+
          (simulate-ranks prior
                          (lambda+ (model-draw)
                            (filled-array 50 (generator model-draw)))
                          (lambda (y prior)
                            (univariate-normal-model
                             (sweep (mean-sse-accumulator) y)
                             prior))
                          (lambda (p)
                            (vector (mean p) (variance p)))))
         (z-statistics (calculate-abs-z-statistics ranks+)))
    (ensure-same p12 p :test #'==)
    (format t "z-statistics: ~A~%" z-statistics)
    (ensure (< (emax z-statistics) 3d0))))

(addtest (samplers-tests)
  lr-kv-dummy-2phase
  (let+ ((k 2)
         (n 10)
         ((&values y x) (cl-random-tests:random-y-x (* 2 n) k))
         (variance 7)
         ;; single step
         (p2 (lr-kv y x variance))
         ;; two steps, first half
         (h1 (cons 0 n))
         (p1 (lr-kv (sub y h1) (sub x h1 t) variance))
         ;; second half, using first half as prior
         (h2 (cons n nil))
         (p2-1 (lr-kv (sub y h2) (sub x h2 t) variance :prior p1)))
    (ensure-same (mean p2) (mean p2-1))
    (ensure-same (variance p2) (variance p2-1))))

(addtest (samplers-tests)
  lr-kv-small
  (let+ ((x (clo 1 1 :/
                 1 2
                 1 3
                 1 4
                 1 5
                 1 6
                 1 7))
         (y (clo 2 2 3 4 5 6 6))
         (sd 19d0)
         (lr (lr-kv y x (expt sd 2)))
         ((&accessors-r/o mean variance) lr)
         (x-t (e/ x sd))
         (y-t (e/ y sd)))
    (ensure-same mean (solve (mm t x-t) (mm (transpose x-t) y-t)))
    (ensure-same variance (invert (mm t x-t)))))

(addtest (samplers-tests)
  multivariate-normal-model
  (let+ ((*lift-equality-test* #'==)
         (k 2)
         (n 10)
         (y (filled-array (list (* 2 n) k) (curry #'random 10d0)
                          'double-float))
         ;; single step
         (p2 (multivariate-normal-model y))
         ;; ;; two steps
         (p1 (multivariate-normal-model (sub y (cons 0 n) t)))
         (p2-1 (multivariate-normal-model (sub y (cons n nil) t)
                                          :prior p1))
         (ranks
          (simulate-ranks p2
                          (lambda+ (model-draw)
                            (combine (filled-array 50
                                                   (generator model-draw))))
                          (lambda (y prior)
                            (multivariate-normal-model y :prior prior))
                          (lambda (p)
                            (concat 'double-float 
                                    (mean p)
                                    (flatten-array 
                                     (as-array (variance p))))))))
    (ensure-same (mean p2) (mean p2-1))
    (ensure-same (nu p2) (nu p2-1))
    (ensure-same (as-matrix (inverse-scale p2)) (inverse-scale p2-1))
    (ensure-same (kappa p2) (kappa p2-1))
    (ensure (small-zs? (calculate-abs-z-statistics ranks)))))
  
