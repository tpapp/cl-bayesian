(in-package :cl-bayesian)

;;;;  Polynomials for time series.
;;;;
;;;;  Coefficients are in increasing order.  These functions are
;;;;  intended for operation on lag polynomials where the first
;;;;  coefficient is 1, which is omitted.

(deftype polynomial (&optional n)
  `(simple-array double-float (,n)))

(defun make-poly (&rest coefficients)
  (make-array (length coefficients) :element-type 'double-float
              :initial-contents coefficients))

(defun poly* (a b)
  "Multiply two polynomials."
  (check-type a polynomial)
;;  (check-type b polynomial)
  (let* ((a-n (length a))
         (b-n (length b)))
    ;; take care of unit polynomials
    (when (zerop a-n)
      (return-from poly* b))
    (when (zerop b-n)
      (return-from poly* a))
    ;; actual multiplication
    (let ((c (make-array (+ a-n b-n) :element-type 'double-float)))
      ;; 1
    (dotimes (a-i a-n)
      (setf (aref c a-i) (aref a a-i)))
      ;; rest
      (dotimes (b-i b-n)
        (let ((b-coeff (aref b b-i)))
          (incf (aref c b-i) b-coeff)
          (iter
            (for c-i :from (1+ b-i))
            (for a-coeff :in-vector a)
            (incf (aref c c-i) (* a-coeff b-coeff)))))
      c)))

(defun filter (x phi)
  "Filter x through the lag polynomial phi."
  (check-type phi polynomial)
  (check-type x (simple-array double-float (*)))
  (let ((x-n (length x))
        (phi-n (length phi)))
    (unless (< (1+ phi-n) x-n)
      (error "series is too short for filter"))
    (when (zerop phi-n)
      (return-from filter x))
    (let ((result (make-array (- x-n phi-n) :element-type 'double-float)))
      ;; 1
      (iter
        (for result-i :from 0)
        (for x-i :from phi-n :below x-n)
        (setf (aref result result-i)
              (+ (aref x x-i)
                 (iter
                   (for x-elt :in-vector x :downfrom (1- x-i))
                   (for phi-elt :in-vector phi)
                   (summing (* x-elt phi-elt))))))
      result)))
