;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package cl-bayesian-tests)

(deftestsuite samplers-tests (cl-bayesian-tests)
  ())

(addtest (samplers-tests)
  lr-kv-dummy-2phase
  (bind ((k 2)
         (n 10)
         ((:values y x) (cl-random-tests:random-y-x (* 2 n) k))
         (variance 7)
         ;; single step
         (p2 (lr-kv y x variance))
         ;; two steps, first half
         (h1 (si 0 n))
         (p1 (lr-kv (sub y h1) (sub x h1 t) variance))
         ;; second half, using first half as prior
         (h2 (si n 0))
         (p2-1 (lr-kv (sub y h2) (sub x h2 t) variance :prior p1)))
    (ensure-same (mean p2) (mean p2-1))
    (ensure-same (variance p2) (variance p2-1))))

(addtest (samplers-tests)
  lr-kv-small
  (bind ((x (clo 1 1 :/
                 1 2
                 1 3
                 1 4
                 1 5
                 1 6
                 1 7))
         (y (clo 2 2 3 4 5 6 6))
         (sd 19d0)
         (lr (lr-kv y x (expt sd 2)))
         ((:accessors-r/o mean variance) lr)
         (x-t (e/ x sd))
         (y-t (e/ y sd)))
    (ensure-same mean (solve (mm t x-t) (mm (transpose x-t) y-t)))
    (ensure-same variance (invert (mm t x-t)))))
