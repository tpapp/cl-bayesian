;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian-tests)

(deftestsuite diagnostics-tests (cl-bayesian-tests) ())

(addtest (diagnostics-tests)
  (let+ ((n 100)
         (zeroes (make-array n :initial-element 0d0))
         (ones (make-array n :initial-element 1d0))
         (twos (make-array n :initial-element 2d0))
         (s1 (concat 'double-float zeroes ones))
         (s2 (concat 'double-float zeroes twos))
         (psrf (psrf (mapcar (curry #'sweep 'sse)
                             (list s1 s2)))))
    (ensure-same (psrf-r psrf) 1.284474)
    ;; (ensure-same r-upper 2.190318)
    ))

(addtest (diagnostics-tests)
  (let+ ((nrow 200)
         (ncol 20)
         (matrix (filled-array (list nrow ncol) (curry #'random 1d0)
                               'double-float))
         (model (gensym))
         (burn-in (floor nrow 3))
         (mcmc-sample (make-instance 'mcmc-sample :model model
                                                  :elements matrix
                                                  :burn-in burn-in))
         (partial-ranges #((20 . 40) (60 . 120) (100 . 180)))
         (lags 5)
         ((&slots-r/o (model2 model) (partial-ranges2 partial-ranges)
                      autocovariance-accumulators partial-accumulators)
          (column-statistics mcmc-sample
                             :partial-ranges partial-ranges
                             :lags lags))
         ((&flet+ sweep-with-accumulators
              ((start . end) accumulator-generator)
            (aprog1 (filled-array ncol accumulator-generator)
              (loop for row-index from start below end do
                (loop for col-index below ncol do
                  (add (aref it col-index)
                       (aref matrix row-index col-index)))))))
         (*lift-equality-test* #'==)
         (partial-matrix (combine partial-accumulators)))
    (ensure-same autocovariance-accumulators
                 (sweep-with-accumulators (cons burn-in nrow)
                                          (curry #'autocovariance-accumulator
                                                 lags)))
    (iter
      (for partial-range in-vector partial-ranges with-index index)
      (ensure-same (sub partial-matrix t index)
                   (sweep-with-accumulators partial-range
                                            #'mean-sse-accumulator)))
    (ensure-same model model2 :test #'eq)
    (ensure-same partial-ranges partial-ranges2
                 :test #'equalp)))

;;; old implementation of psrf working directly with sequences saved here for
;;; comparison and testing purposes
;;; 
;; (defun psrf-direct (sequences &key (confidence 0.975d0) skip-length-check?)
;;   "Estimate the potential scale reduction factor.  Algorithm is from Gelman and
;; Rubin (1992), but the degrees of freedom correction is according to Brooks and
;; Gelman (1998)."
;;   ;; !!! should return the upper limit of the confidence interval as the second
;;   ;; value.  Since the F distribution is not implemented yet, this functionality
;;   ;; is not available now.
;;   (declare (ignore confidence))
;;   (let ((m (length sequences))
;;         (n (length (aref sequences 0))))
;;     (unless skip-length-check?
;;       (assert (every (lambda (sequence) (= (length sequence) n))
;;                      (subseq sequences 1))))
;;     (iter
;;       (for sequence :in-vector sequences)
;;       (let ((mean (mean sequence)))
;;         (collecting mean :into means :result-type vector)
;;         (collecting (variance sequence)
;;                     :into variances :result-type vector))
;;       (finally
;;        (let* ((mu (mean means))
;;               (b (* n (variance means)))
;;               (w (mean variances))
;;               (var-b (/ (* 2 (expt b 2)) (1- m)))
;;               (var-w (/ (variance variances) m))
;;               (1+1/m (1+ (/ m)))
;;               (n-1 (1- n))
;;               (V (/ (+ (* n-1 w) (* 1+1/m b)) n))
;;               (var-V (/ (+ (* (expt n-1 2) var-w)
;;                            (* (expt 1+1/m 2) var-b)
;;                            (* 2 1+1/m n-1 (/ n m)
;;                               (- (covariance-xy variances (eexpt means 2))
;;                                  (* 2 mu (covariance-xy variances means)))))
;;                         (expt n 2)))
;;               (df (/ (* 2 (expt V 2)) var-V))
;;               (df-adj (/ (+ df 3) (1+ df)))
;;               ;; (b-df (1- m))
;;               ;; (w-df (/ (* 2 (expt w 2)) var-w))
;;               (R^2-fixed (/ n-1 n))
;;               (R^2-random (* (/ 1+1/m n) (/ b w))))
;;          (d:p "we are in PSRF2~%")
;;          (d:v n means variances b var-b var-w var-V df-adj)
;;          (return (sqrt (* df-adj (+ R^2-fixed R^2-random)))))))))

