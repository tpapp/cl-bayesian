(in-package :bayesian-inference)

;;;; !!! all should be prepended with ts

(defclass date ()
  ((main :initarg :main :accessor main)
   (sub :initarg :sub :accessor sub :initform 0)
   (frequency :initarg :frequency :initform 1 :reader frequency)))

(defmethod initialize-instance :after ((date date) &key &allow-other-keys)
  (with-slots (sub frequency) date
    ;; check frequency and sub
    (assert (plusp frequency))
    (assert (and (<= 0 sub) (< sub frequency))))
  date)

(defmethod print-object ((date date) stream)
  (print-unreadable-object (date stream :type t)
    (with-slots (main sub frequency) date
      (format stream "~a ~a/~a" main sub frequency))))

(defun date->integer (date)
  (with-slots (main sub frequency) date
    (+ (* main frequency) sub)))

(defun integer->date (num frequency)
  (multiple-value-bind (main sub) (floor num frequency)
    (make-instance 'date :main main :sub sub :frequency frequency)))

(defun date+ (date offset)
  (integer->date (+ (date->integer date) offset) (frequency date)))

(defun date- (date offset)
  (integer->date (+ (date->integer date) offset) (frequency date)))

(defun frequency= (&rest dates)
  (apply #'= (mapcar #'frequency dates)))

(defun date-diff (date1 date2)
  (assert (frequency= date1 date2))
  (- (date->integer date1) (date->integer date2)))

(defun date-min (&rest dates)
  (assert (apply #'frequency= dates))
  (integer->date (apply #'min (mapcar #'date->integer dates))
		 (frequency (first dates))))

(defun date-max (&rest dates)
  (assert (apply #'frequency= dates))
  (integer->date (apply #'max (mapcar #'date->integer dates))
		 (frequency (first dates))))

;;;; !! time series
;;;; !!! document everything

(defclass ts (date)
  ((vec :initarg :vec :accessor vec)))

(defun end-date (ts)
  (date+ ts (1- (length (vec ts)))))

(defmethod print-object ((ts ts) stream)
  (print-unreadable-object (ts stream :type t)
    (with-slots (main sub frequency vec) ts
      (format stream "~a ~a/~a+~a: ~a" main sub frequency (length vec) vec))))

(defun intersect (&rest tss)
  (declare (optimize (debug 3)))
  (assert (apply #'frequency= tss))
  (bind (((max-start common-length)
	  (iter
	    (for ts :in tss)
	    (for start := (date->integer ts))
	    (maximize start :into max-start)
	    (minimize (+ start (length (vec ts))) :into min-end)
	    (finally (return (list max-start (- min-end max-start))))))
	 (frequency (frequency (car tss))))
    (with-slots (main sub) (integer->date max-start frequency)
      (mapcar (lambda (ts)
		(with-slots (vec) ts
		  (let ((start-index (- max-start (date->integer ts))))
		    (make-instance 'ts :main main :sub sub :frequency frequency
				   :vec (subseq vec start-index 
						(+ start-index common-length))))))
	      tss))))

;; (defparameter *d1* (make-instance 'date :main 2004 :sub 8 :frequency 12))
;; (defparameter *d2* (make-instance 'date :main 2003 :sub 4 :frequency 12))

;; (date+ *d1* 16)
;; (date-diff *d1* *d2*)

;; (defparameter *ts1* (make-instance 'ts 
;; 				   :main 2004 :sub 0 :frequency 4 
;; 				   :vec #(0 1 2 3 4 5 6 7 8 9 10 11)))

;; (defparameter *ts2* (make-instance 'ts 
;; 				   :main 2005 :sub 2 :frequency 4 
;; 				   :vec #(0 1 2 3 4 5 6 7 8 9 10 11)))

;; (intersect *ts1* *ts2*)
