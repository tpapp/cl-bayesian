;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

(defun slice-sample-so (x g &key (w 1d0) (max-iter 100) lower upper
                        (gx (funcall g x)))
  "Slice sampling, starting at X, using the stepping out algorithm of
Neal (2003, especially Fig 3 and 5).  Return the new sample point and the
value of the log density at that point.  G is the log of the probability (+
constant), and may return NIL when the probability is zero.  W is the starting
width.  MAX-ITER, LOWER and UPPER define maximum number of iterations, and
lower/upper bounds for X.  When NIL, these three are considered as not
applicable and are not used to terminate the algorithm.  In particular, when
bounds are given, G is never called outside these."
  (assert gx () "p(x) = 0")
  (bind ((log-y (- gx (draw-standard-exponential)))
         (u (random 1d0))
         (left (- x (* w u)))
         (right (+ left w))
         ((:values j k) (if max-iter
                            (let* ((j (floor (random (float max-iter 1d0))))
                                   (k (- max-iter 1 j)))
                              (values j k))
                            (values nil nil)))
         ((:flet outside? (z)) (awhen (funcall g z) (<= it log-y))))
    ;; extend to the left
    (loop
      (when (or (and j (<= j 0)) (and lower (<= left lower)) (outside? left))
        (return))
      (decf left w)
      (when j (decf j)))
    ;; extend to the right
    (loop
      (when (or (and k (<= k 0)) (and upper (>= right upper)) (outside? right))
        (return))
      (incf right w)
      (when k (decf k)))
    ;; when bounds are given, shrink
    (when upper
      (minf right upper))
    (when lower
      (maxf left lower))
    ;; sample slice
    (loop
      (let* ((x1 (+ left (random (- right left))))
             (gx1 (funcall g x1)))
        (cond
          ;; note: we use <= for termination, following the code of Neal
          ;; instead of the paper, to ensure termination if the interval is
          ;; shrunk to the original point X
          ((and gx1 (<= log-y gx1)) (return (values x1 gx1)))
          ;; just shrink interval
          ((< x1 x) (setf left x1))
          (t (setf right x1)))))))

;; (defun test-slice-sample-so-univariate (g n &key (burn-in 1000)
;;                                         sampler-params (x 1d0))
;;   (bind ((g-counter 0)
;;          ((:flet g-c (x))
;;           (incf g-counter)
;;           (funcall g x))
;;          ((:flet update ())
;;           (setf x (apply #'slice-sample-so x #'g-c sampler-params)))
;;          (s (make-array n :element-type 'double-float)))
;;     (loop repeat burn-in do (update))
;;     (setf g-counter 0)
;;     (loop for i :below n do (setf (aref s i) (update)))
;;     (format t "~A calls/sample~%" (float (/ g-counter n)))
;;     s))

;; (defun test-sso-with-rv (rv &key (n 100000) (burn-in 10000)
;;                          sampler-params (x (mean rv)))
;;   (bind (((:flet g (x))
;;           (log-pdf rv x t))
;;          (s (test-slice-sample-so-univariate #'g n :x x :burn-in burn-in
;;                                              :sampler-params sampler-params)))
;;     (d:v (mean s) (mean rv) (variance s) (variance rv))
;;     (values (/ (- (mean s) (mean rv)) (sqrt (variance rv)))
;;             (/ (variance s) (variance rv)))))

;; (test-sso-with-rv (make-instance 'normal :mu 0d0 :sigma 1d0))

;; (test-sso-with-rv (make-instance 'gamma :alpha 9d0 :beta 3d0))
