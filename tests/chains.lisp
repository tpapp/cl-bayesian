;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian-tests)

(deftestsuite diagnostics-tests (cl-bayesian-tests) ())

(addtest (diagnostics-tests)
  psrf-test
  ;; test the PSRF calculations on a known result (from R)
  (let+ ((n 100)
         (zeroes (make-array n :initial-element 0d0))
         (ones (make-array n :initial-element 1d0))
         (twos (make-array n :initial-element 2d0))
         (s1 (concat 'double-float zeroes ones))
         (s2 (concat 'double-float zeroes twos))
         (psrf (calculate-psrf (mapcar (curry #'sweep 'sse)
                                       (list s1 s2)))))
    (ensure-same (psrf-r psrf) 1.284474)
    ;; (ensure-same r-upper 2.190318)
    ))

(addtest (diagnostics-tests)
  mcmc-statistics-test
  ;; calculate statistic for random elements
  (let+ ((elements (filled-array '(50 5) (curry #'random 1d0) 'double-float))
         (burn-in 20)
         (columns (subarrays 1 (transpose (sub elements (cons burn-in nil) t))))
         (lags 4)
         (model (gensym))
         (sample (make-instance 'mcmc-sample :model model :elements elements))
         (accumulator-generator #'mean-sse-accumulator)
         (statistics 
          (mcmc-statistics sample :divisions 3 :minimum-length 0 :lags lags
                                  :accumulator-generator accumulator-generator
                                  :burn-in-fraction (/ burn-in 
                                                       (nrow elements))))
         ((&slots-r/o sse-ranges) statistics)
         (*lift-equality-test* #'==))
    ;; accumulators for columns
    (ensure-same (accumulators statistics)
                 (map1 (lambda (s) (sweep (funcall accumulator-generator) s))
                       columns))
    ;; sse accumulators
    (iter
      (for column :in-vector (subarrays 1 (transpose elements)))
      (for sse-accumulator :in-vector (sse-accumulators statistics))
      (let* ((accumulators (map 'vector (lambda+ ((start . end))
                                          (sweep (mean-sse-accumulator)
                                                 (subseq column start end)))
                                sse-ranges)))
        (ensure-same sse-accumulator accumulators)))
    ;; autocovariance accumulators
    (ensure-same (autocovariance-accumulators statistics)
                 (map1 (lambda (v) 
                         (sweep (autocovariance-accumulator lags) v))
                       columns))))

(addtest (diagnostics-tests)
  mcmc-statistics-test2
  ;; In this test the period after burn-in is not composed of the apparent
  ;; sse-ranges, the purpose of this test is to see if the statistics are
  ;; calculated correctly.
  (let+ ((nrow 200)
         (ncol 2)
         (matrix (filled-array (list nrow ncol) (curry #'random 1d0)
                               'double-float))
         (model (gensym))
         (burn-in (floor nrow 3))
         (mcmc-sample
          (make-instance 'mcmc-sample :model model :elements matrix))
         (sse-ranges #((20 . 40) (60 . 120) (100 . 180)))
         (lags 5)
         ((&slots-r/o (model2 model) (sse-ranges2 sse-ranges)
                      autocovariance-accumulators sse-accumulators)
          (mcmc-statistics mcmc-sample :sse-ranges sse-ranges :lags lags
                                       :burn-in-fraction (/ burn-in nrow)))
         ((&flet+ sweep-with-accumulators
              ((start . end) accumulator-generator)
            (aprog1 (filled-array ncol accumulator-generator)
              (loop for row-index from start below end do
                (loop for col-index below ncol do
                  (add (aref it col-index)
                       (aref matrix row-index col-index)))))))
         (*lift-equality-test* #'==)
         (partial-matrix (combine sse-accumulators)))
    (iter
      (for sse-range :in-vector sse-ranges :with-index index)
      (ensure-same (sub partial-matrix t index)
                   (sweep-with-accumulators sse-range
                                            #'mean-sse-accumulator)))
    (ensure-same model model2 :test #'eq)
    (ensure-same sse-ranges sse-ranges2
                 :test #'equalp)
    (ensure-same autocovariance-accumulators
                 (sweep-with-accumulators (cons burn-in nrow)
                                          (curry #'autocovariance-accumulator
                                                 lags)))
    partial-matrix))

(addtest (diagnostics-tests)
  mcmc-summary-test
  ;; testing mcmc summaries
  (let+ ((model (gensym))
         (m 2)                          ; number of variables
         (n 400)                        ; total length
         (burn-in 200)
         (burn-in-fraction (/ burn-in n))
         (lag 10)
         ((&flet make-sample (stencil &key (model model) (burn-in burn-in)
                                      (m m) (n n))
            (let ((elements (make-array (list n m)))
                  (stencil (coerce stencil 'vector))
                  (stencil-length (length stencil)))
              (dotimes (index (array-total-size elements))
                (setf (row-major-aref elements index)
                      (aref stencil (mod index stencil-length))))
              (make-instance 'mcmc-sample :elements elements :model model))))
         (samples (mapcar #'make-sample
                          '((0 -1 0 0)
                            (1 2 0 0)
                            (3 5 7 11 13 17))))
         (crude-mean 
          (mean (subarrays 1 (stack* t :v 
                                     (mapcar (lambda (s)
                                               (sub (elements s)
                                                    (cons burn-in nil) t))
                                             samples)))))
         ((&flet sample-autocorrelations (sample)
            (map1 (rcurry #'autocorrelations lag)
                  (subarrays 1 (transpose (sub (elements sample)
                                               (cons burn-in nil) t))))))
         (crude-autocorrelations
          (map1 #'mean
                (subarrays 1
                           (transpose
                            (combine (map 'vector #'sample-autocorrelations
                                          samples))))))
         (statistics 
          (mapcar (lambda (s) 
                    (mcmc-statistics s :burn-in-fraction burn-in-fraction))
                  samples))
         (summary (summarize-mcmc-statistics statistics))
         ((&flet summarize-incompatible-chains (&rest arguments)
            (summarize-mcmc-statistics
             (list (first statistics)
                   (scalar-statistics (apply #'make-sample '(1)
                                             arguments))))))
         (*lift-equality-test* #'==))
    ;; check mean and autocovariance
    (ensure-same (map1 #'mean (accumulators summary)) crude-mean)
    (ensure-same (mean-autocorrelations summary) crude-autocorrelations)
    ;; incompatible chains should give errors
    (ensure-error (summarize-incompatible-chains :model 'foo))
    ;; (ensure-error (summarize-incompatible-chains :burn-in 111))
    (ensure-error (summarize-incompatible-chains :m 9))
    (ensure-error (summarize-incompatible-chains :n 17))))

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

