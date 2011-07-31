;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package cl-bayesian-tests)

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

(defun simulate-ranks (prior simulate estimate &key (m 40) (n 500)
                                                    (draw-vector #'draw))
  (let+ ((theta0 (funcall draw-vector prior))
         (y (funcall simulate m theta0))
         (posterior (funcall estimate y prior))
         (theta+ (filled-array n (lambda ()
                                   (funcall draw-vector posterior)))))
    (calculate-empirical-ranks (ensure-vector theta0)
                               (ensure-matrix (combine theta+) :column))))



(addtest (samplers-tests)
  (let+ (((&values p12 p) (2phase-posterior 10 #'mean-sse-accumulator
                                            #'univariate-normal-error))
         (#(rank)
           (simulate-ranks (r-inverse-chi-square 9 2d0)
                           (lambda (m v)
                             (filled-array m (generator (r-normal 0 v))))
                           (lambda (y prior)
                             (univariate-normal-error 
                              (sweep (mean-sse-accumulator) y)
                              prior)))))
    (ensure-same p12 p :test #'==)
    (ensure (< 0.01 rank 0.99))))

(addtest (samplers-tests)
  (let+ (((&values p12 p) (2phase-posterior 10 #'mean-sse-accumulator
                                            #'univariate-normal-model))
         (prior (univariate-normal-model
                 (sweep (mean-sse-accumulator) #(ivec 10))))
         (#(rank1 rank2)
           (simulate-ranks prior
                           (lambda+ (m #(mu s^2))
                             (filled-array m (generator (r-normal mu s^2))))
                           (lambda (y prior)
                             (univariate-normal-model
                              (sweep (mean-sse-accumulator) y)
                              prior))
                           :draw-vector (lambda (p)
                                          (let+ (((&values mu s^2) (draw p)))
                                            (vector mu s^2))))))
    (ensure-same p12 p :test #'==)
    (ensure (< 0.01 rank1 0.99))
    (ensure (< 0.01 rank2 0.99))
    (vector rank1 rank2)))

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
                                          :prior p1)))
    (ensure-same (mean p2) (mean p2-1))
    (ensure-same (nu p2) (nu p2-1))
    (ensure-same (as-matrix (inverse-scale p2)) (inverse-scale p2-1))
    (ensure-same (kappa p2) (kappa p2-1))))
  
