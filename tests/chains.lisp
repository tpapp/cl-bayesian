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

;;; We test by calculating the statistics for IID random elements.

(defstruct iid-model
  "For testing MCMC diagnostics."
  (n 5)
  (generator (generator (r-normal))))

(defmethod scalar-parameters-layout ((model iid-model))
  (array-layout (iid-model-n model)))

(defmethod draw ((model iid-model) &key)
  ;; usually we there are no DRAW methods for models, but in the IID case it
  ;; makes sense
  (let+ (((&structure-r/o iid-model- n generator) model))
    (make-iid-state :model model :elements (generate-array n generator))))

(defstruct iid-state
  "For testing MCMC diagnostics.  Of course not a state in the actual sense."
  model
  elements)

(defmethod model ((state iid-state))
  (iid-state-model state))

(defmethod draw ((state iid-state) &key)
  (let+ (((&structure-r/o iid-state- model) state))
    (draw model)))

(defmethod scalar-parameters ((state iid-state) &key copy?)
  (maybe-copy-array (iid-state-elements state) copy?))

(addtest (diagnostics-tests)
  mcmc-statistics-test
  (let+ ((model (make-iid-model))
         (sample (draw-chain (draw model) 50 :stream nil))
         (burn-in 20)
         (columns*  (subarrays 1 (transpose 
                                (combine (map1 #'iid-state-elements
                                               sample)))))
         (columns (map1 (lambda (c) (subseq c burn-in)) columns*))
         (lags 4)
         (accumulator-generator #'mean-sse-accumulator)
         (statistics 
          (mcmc-statistics sample :divisions 3 :minimum-length 0 :lags lags
                                  :accumulator-generator accumulator-generator
                                  :burn-in-fraction (/ burn-in (length sample))))
         ((&slots-r/o sse-ranges) statistics)
         (*lift-equality-test* #'==))
    ;; accumulators for columns
    (ensure-same (accumulators statistics)
                 (map1 (lambda (s) (sweep (funcall accumulator-generator) s))
                       columns))
    ;; sse accumulators
    (iter
      (for column :in-vector columns*)
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
  (let+ ((n-sample 200)
         (n-parameters 2)
         (model (make-iid-model :n n-parameters))
         (sample (draw-chain (draw model) n-sample :stream nil))
         (burn-in (floor n-sample 3))
         (sse-ranges #((20 . 40) (60 . 120) (100 . 180)))
         (lags 5)
         ((&slots-r/o (model2 model) (sse-ranges2 sse-ranges)
                      autocovariance-accumulators sse-accumulators)
          (mcmc-statistics sample :sse-ranges sse-ranges :lags lags
                                  :burn-in-fraction (/ burn-in n-sample)))
         ((&flet+ sweep-with-accumulators
              ((start . end) accumulator-generator)
            (let ((acc (generate-array n-parameters accumulator-generator)))
              (loop for sample-index from start below end do
                (loop for p across (iid-state-elements
                                    (aref sample sample-index))
                      for a across acc
                      do (add a p)))
              acc)))
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
                 (sweep-with-accumulators (cons burn-in n-sample)
                                          (curry #'autocovariance-accumulator
                                                 lags)))
    partial-matrix))

(addtest (diagnostics-tests)
  mcmc-summary-test
  ;; testing mcmc summaries for univariate samples built from stencils
  (let+ ((burn-in-fraction 0.4)
         (lag 10)
         (model (make-iid-model :n 1))
         ((&flet make-sample (stencil &key (model model) (n 600))
            "Return a sample of length N by repeating STENCIL."
            (iter
              (with stencil := (map 'vector #'vector stencil))
              (with stencil-length := (length stencil))
              (for index :below n)
              (collect 
                  (make-iid-state :model model
                                  :elements (aref stencil (mod index stencil-length)))
                :result-type vector))))
         (samples (mapcar #'make-sample
                          '((0 -1 0 0)
                            (1 2 0 0)
                            (3 5 7 11 13 17))))
         (columns (mapcar (lambda (s)
                            (map1 (compose #'first* #'iid-state-elements)
                                  (subseq s
                                          (cl-bayesian::calculate-burn-in
                                           (length s) burn-in-fraction))))
                          samples))
         (crude-mean (mean (stack* t :v columns)))
         (crude-autocorrelations (mean (map1 (rcurry #'autocorrelations lag)
                                             columns)))
         (statistics 
          (mapcar (lambda (s) 
                    (mcmc-statistics s :burn-in-fraction burn-in-fraction))
                  samples))
         (summary (summarize-mcmc-statistics statistics))
         ((&flet summarize-incompatible-chains (&rest arguments)
            (summarize-mcmc-statistics
             (list (first statistics)
                   (mcmc-statistics (apply #'make-sample '(1) arguments)
                                    :burn-in-fraction burn-in-fraction)))))
         (*lift-equality-test* #'==))
    ;; check mean and autocovariance
    (ensure-same (map1 #'mean (accumulators summary)) (vector crude-mean))
    (ensure-same (mean-autocorrelations summary) (vector crude-autocorrelations))
    ;; incompatible chains should give errors
    (ensure-error (summarize-incompatible-chains :model (make-iid-model :n 2)))
    (ensure-error (summarize-incompatible-chains :n 977))))


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

