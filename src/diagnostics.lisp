;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

(defun potential-scale-reduction (sequences &optional (confidence 0.975d0))
  "Estimate the potential scale reduction factor.  Algorithm is from Gelman and
Rubin (1992), but the degrees of freedom correction is according to Brooks and
Gelman (1998)."
  ;; !!! should return the upper limit of the confidence interval as the second
  ;; value.  Since the F distribution is not implemented yet, this functionality
  ;; is not available now.
  (declare (ignore confidence))
  (let* ((sequences (coerce sequences 'list))
         (m (length sequences))
         (n (length (first sequences))))
    (assert (every (lambda (sequence) (= (length sequence) n)) (cdr sequences)))
    (iter
      (for sequence :in sequences)
      (let ((mean (mean sequence)))
        (collecting mean :into means :result-type vector)
        (collecting (sample-var sequence mean)
                    :into variances :result-type vector))
      (finally
       (let* ((mu (mean means))
              (b (* n (sample-var means mu)))
              (w (mean variances))
              (var-b (/ (* 2 (expt b 2)) (1- m)))
              (var-w (/ (sample-var variances w) m))
              (1+1/m (1+ (/ m)))
              (n-1 (1- n))
              (V (/ (+ (* n-1 w) (* 1+1/m b)) n))
              (var-V (/ (+ (* (expt n-1 2) var-w)
                           (* (expt 1+1/m 2) var-b)
                           (* 2 1+1/m n-1 (/ n m)
                              (- (sample-cov variances (eexpt means 2))
                                 (* 2 mu (sample-cov variances means)))))
                        (expt n 2)))
              (df (/ (* 2 (expt V 2)) var-V))
              (df-adj (/ (+ df 3) (1+ df)))
              ;; (b-df (1- m))
              ;; (w-df (/ (* 2 (expt w 2)) var-w))
              (R^2-fixed (/ n-1 n))
              (R^2-random (* (/ 1+1/m n) (/ b w))))
         (return (sqrt (* df-adj (+ R^2-fixed R^2-random)))))))))

