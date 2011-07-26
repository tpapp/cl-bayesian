;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian-tests)

(deftestsuite dlm-tests (cl-bayesian-tests) ())

(addtest (dlm-tests)
  uddu
  (let* ((a (mm t (clo :double 1 2 :/ 3 4)))
         (a-uddu (uddu a))
         (h (mm t (clo :double 3 5 :/ 7 11)))
         (g (clo :double 9 2 :/ 17 4))
         (*lift-equality-test* #'==))
    (ensure-same (as-array a-uddu :hermitian? t) a)
    (ensure-same (as-array a-uddu) (as-array a))
    (ensure-same (as-array (uddu-update a-uddu h) :hermitian? t)
                 (e+ a h))
    (ensure-same (as-array (uddu-update a-uddu h t) :hermitian? t)
                 (e+ a h))
    (ensure-same (as-array (uddu-multiply-update a g h) :hermitian? t)
                 (convert-matrix :hermitian (e+ (mmm g a (transpose g)) h)))))

(addtest (dlm-tests)
  ;; checked against results from R/dlm, univariate case
  (let+ ((parameters (make-dlm-parameters :g 1 :W 0.7 :mu 0 :F 1 :V 0.5))
         (y (make-array 10 :initial-element #(1)))
         ((&values m C-inverse a)
          (dlm-forward-filtering #(0) (dlm-parameters-W parameters) y
                                 (make-array 10 :initial-element parameters)))
         (*lift-equality-test* #'==))
    (ensure-same (map1 #'first* m)
                 #(0.5833333 0.8603352 0.9544295 0.9851741 0.9951780
                   0.9984318 0.9994900 0.9998341 0.9999461 0.9999825))
    (ensure-same (map1 #'first* a)
                 #(0.0000000 0.5833333 0.8603352 0.9544295 0.9851741
                   0.9951780 0.9984318 0.9994900 0.9998341 0.9999461))
    (ensure-same (map1 (compose #'/ #'first* #'elements #'uddu-d) C-inverse)
                 #(0.5400617 0.5765434 0.5803942 0.5808015 0.5808446
                   0.5808491 0.5808496 0.5808497 0.5808497 0.5808497))
    (ensure (every (curry #'== #2A((1d0)))
                   (map1 #'uddu-u C-inverse)))))

(addtest (dlm-tests)
  ;; checked against results from R/dlm, bivariate case
  (let+ ((parameters (make-dlm-parameters :g (clo 1 1 :/
                                                  0 1)
                                          :W (clo :diagonal 0.7 0.9)
                                          :mu (clo 0 0)
                                          :F (clo 1 1)
                                          :V 0.5))
         (y (make-array 10 :initial-element #(1)))
         ((&values m C-inverse a)
          (dlm-forward-filtering #(0 0) (dlm-parameters-W parameters) y
                                 (make-array 10 :initial-element parameters)))
         (C (map1 #'invert C-inverse))
         (*lift-equality-test* #'==))
    (ensure-same (combine m)
                 (clo :double
                      0.3333333  0.4285714286 :/
                      0.6898470  0.3379694019
                      0.8907515  0.1594623492
                      0.9720953  0.0563039370
                      0.9968717  0.0145656823
                      1.0017594  0.0017499744
                      1.0015519 -0.0008421624
                      1.0007590 -0.0007769082
                      1.0002779 -0.0003851931
                      1.0000760 -0.0001424831))
    (ensure-same (combine a)
                 (clo :double
                   0.0000000  0.0000000000 :/
                   0.7619048  0.4285714286
                   1.0278164  0.3379694019
                   1.0502138  0.1594623492
                   1.0283992  0.0563039370
                   1.0114373  0.0145656823
                   1.0035094  0.0017499744
                   1.0007098 -0.0008421624
                   0.9999821 -0.0007769082
                   0.9998927 -0.0003851931))
    (ensure-same (combine (map1 (compose #'elements #'uddu-d) C))
                 (clo :double
                      0.4353537 0.8896176 :/
                      0.4612531 1.0088823
                      0.4634060 1.0358761
                      0.4637999 1.0404274
                      0.4638920 1.0409681
                      0.4639081 1.0410140
                      0.4639101 1.0410177
                      0.4639103 1.0410184
                      0.4639103 1.0410185
                      0.4639103 1.0410186))
    (ensure-same (combine (map1 #'as-array C))
                 (reshape '(10 2 2)
                          #(0.4666667 -0.3000000
                            -0.3000000  0.5142857
                            0.5909597 -0.4018081
                            -0.4018081  0.6396384
                            0.6155682 -0.4282114
                            -0.4282114  0.6722163
                            0.6189444 -0.4326606
                            -0.4326606  0.6786551
                            0.6192684 -0.4331623
                            -0.4331623  0.6795420
                            0.6192930 -0.4332006
                            -0.4332006  0.6796280
                            0.6192960 -0.4332034
                            -0.4332034  0.6796345
                            0.6192967 -0.4332040
                            -0.4332040  0.6796353
                            0.6192969 -0.4332042
                            -0.4332042  0.6796354
                            0.6192969 -0.4332042
                            -0.4332042  0.6796355)))))


;;; functions below are for testing

(defun dlm-generate-sample (a R parameters+)
  "Generate a sample for a DLM with given parametes, drawing the first state
from N(a,R).  Return (values theta+ y+)."
  (let* ((n (length parameters+))
         (theta+ (make-array n))
         (y+ (make-array n)))
    (dotimes (index n)
      (let+ (((&structure-r/o dlm-parameters- G W mu F V)
              (aref parameters+ index))
             (theta (if (zerop index)
                        (draw (r-multivariate-normal a R))
                        (e+ (mm G (aref theta+ (1- index)))
                            (draw (r-multivariate-normal mu W))))))
        (setf (aref y+ index)
              (draw (r-multivariate-normal (mm F theta) V))
              (aref theta+ index) theta)))
    (values theta+ y+)))

(defun dlm-flatten-theta (theta+)
  (flatten-array (combine theta+ 'double-float) :copy? t))

(defun dlm-simulated-ranks (n a R parameters+)
  "Empirical ranks of simulated data."
  (let+ (((&values theta0 y+) (dlm-generate-sample a R parameters+))
         (draws (combine
                 (filled-array n
                               (lambda ()
                                 (dlm-flatten-theta
                                  (dlm-ff-bs a R y+ parameters+)))))))
    (calculate-empirical-ranks (dlm-flatten-theta theta0) draws)))



;;; compare to R

;;; univariate
(defparameter *parameters*
  (make-dlm-parameters :g (clo :double 1 :/)
                       :mu (clo :double 0)
                       :W (clo :hermitian :double 0.05 :/)
                       :F (clo :double 1 :/)
                       :V (clo :hermitian 0.02 :/)))
(defparameter *prior-a* (clo :double 0))
(defparameter *prior-R* (clo :hermitian :double 2 :/))
(defparameter *parameters+* (make-array 10 :initial-element *parameters*))
(defparameter *parameters+* (make-array 1 :initial-element *parameters*))



;;; locally linear
(defparameter *prior-a* (clo :double 0 0))
(defparameter *prior-R* (clo :diagonal :double 1 1))
(defparameter *parameters+*
  (make-array 10
              :initial-element (make-dlm-parameters 
                                :g (clo :double
                                        1 1 :/ 
                                        0 1)
                                :mu (clo :double 0 0)
                                :w (clo :double :diagonal 0.05 0.03)
                                :f (clo :double 1 1 :/)
                                :v (clo :double :diagonal 0.02))))

(mean (filled-array 500 
                    (lambda ()
                      (first* (dlm-generate-sample *prior-a* *prior-r*
                                                   *parameters+*)))))

(dlm-g)

(defparameter *ranks* 
  (filled-array 50
                (lambda ()
                  (format t ".")
                  (dlm-simulated-ranks 200 *prior-a* *prior-r*
                                       *parameters+*))))

(defparameter *p* (calculate-p-statistics *ranks*))
(defparameter *z* (calculate-abs-z-statistics *ranks*))





(use-package '(:cl-2d :cl-colors))


(loop repeat 200 do
 (with-pdf-frame (frame "/tmp/foo.pdf")
   (defparameter *theta-s+*
     (dlm-ffbs *prior-a* *prior-R* *y+* *parameters+*))
   (plot-rows frame (ivec (length *theta+*))
              (map1 #'first* (combine (vector *theta+* *y+* *theta-s+*)))
              (vector (make-instance 'line-style :color +red+)
                      (make-instance 'line-style :color +green+)
                      (make-instance 'line-style :color +blue+))))
 (sleep 0.1))

(mm
 #S(UDDU
     :U #2A(#)
     :D (clo :diagonal 0.14072))
 #2A((1.0d0)))


