;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

;;; Validation method of Cook et al (2006).

(defun calculate-empirical-ranks (parameters draws)
  "Return a vector of empirical ranks for each parameter.  Note: ranks are
corrected by 0.5 and divided by n+1 to ensure that they are all in (0,1)."
  (let+ (((nrow ncol) (array-dimensions draws))
         (counts (make-array ncol :element-type 'fixnum :initial-element 0)))
    (assert (= ncol (length parameters)))
    (dotimes (row-index nrow)
      (dotimes (col-index ncol)
        (when (< (aref draws row-index col-index) (aref parameters col-index))
          (incf (aref counts col-index)))))
    (map 'vector (lambda (c) (/ (+ c 0.5d0) (1+ nrow))) counts)))

(defun calculate-p-statistic (ranks)
  "Given the RANKS for a single parameter, return a p-value.  Note that all
ranks have to be in (0,1)."
  (cdf (r-chi-square (length ranks))
       (iter
         (for r :in-vector ranks)
         (summing (expt (quantile (r-normal) r) 2)))))

(defun abs-z-transform (p)
  "Transform the probability p (in [0,1]) to the absolute value of a standard
  normal.  For 0 and 1, return NIL."
  (if (or (= p 1) (= p 0))
      nil
      (abs (quantile (r-normal) p))))

(defun calculate-p-statistics (ranks+)
  "Calculate the P statistics of a ranks (a vector of vectors, or equal
length)."
  (map1 #'calculate-p-statistic
        (subarrays 1 (transpose (combine ranks+)))))

(defun calculate-abs-z-statistics (ranks+)
  "Calculate the abs(z) statistics from ranks. "
  (map1 #'abs-z-transform (calculate-p-statistics ranks+)))

;;; testing the validation with a normal distribution

;; (defun simulate-linear-regression-y (prior x)
;;   "Return (valyes y parameters)."
;;   (let+ (((&values beta sigma) (draw prior)))
;;     (values
;;      (e+ (mm x beta) (filled-array (nrow x)
;;                                    (generator (r-normal 0 (sqrt sigma)))))
;;      (concat 'double-float (vector sigma) beta))))

;; (defun simulate-linear-regression-parameters (n y x prior)
;;   "Given Y, return a matrix of draws."
;;   (let ((lr (linear-regression y x :prior prior)))
;;     (combine (filled-array n (lambda ()
;;                                (let+ (((&values beta sigma) (draw lr)))
;;                                  (concat 'double-float (vector sigma) beta))))
;;              'double-float)))


;; (defparameter *q* (calculate-empirical-ranks *parameters* *draws*))
;; (map1 #'float *q*)

;; (defparameter *ranks+* 
;;   (let ((n 1000)
;;         (x *x*)
;;         (prior *prior*))
;;     (filled-array 100 (lambda ()
;;                         (let+ (((&values y parameters)
;;                                 (simulate-linear-regression-y prior x))
;;                                (draws (simulate-linear-regression-parameters
;;                                        n y x prior)))
;;                           (calculate-empirical-ranks parameters draws))))))

;; (histogram (sub (map1 #'float (combine *ranks+*)) t 0)
;;            (rcurry #'scott-rule :correction 0.25))

;; (defparameter *r* (sub (map1 #'float (combine *ranks+*)) t 0))
;; (range *r*)

;; (calculate-empirical-p-statistic (sub (map1 #'float (combine *ranks+*)) t 0))

;; (defparameter *p* (calculate-p-statistics *ranks+*))
