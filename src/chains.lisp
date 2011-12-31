;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

(defparameter *suggested-minimum-burn-in* 200
  "The suggested minimum number of observations discarded as burn-in.  Using
less than this generates a warning.")

(defparameter *default-burn-in-fraction* 0.5
  "The default fraction of observations discarded as burn-in.")

(defun check-burn-in (burn-in)
  "Check that burn-in is above the suggested minimum (when defined).  Does not
return a value, called for side effects (conditions)."
  (when (and *suggested-minimum-burn-in*
             (< burn-in *suggested-minimum-burn-in*))
    (warn "Burn-in ~A is below suggested minimum burn-in (~A)."
          burn-in *suggested-minimum-burn-in*))
  (values))

(defun calculate-burn-in (n &optional (burn-in-fraction *default-burn-in-fraction*))
  "Calculate burn-in from burn-in fraction and total number of samples.
Note: efined as a separate function for consistency of rounding."
  (assert (within? 0 burn-in-fraction 1))
  (ceiling (* n burn-in-fraction)))

(defun discard-burn-in (samples &optional (burn-in-fraction *default-burn-in-fraction*))
  "Discard the burn-in from the samples."
  (subseq samples (calculate-burn-in (length samples) burn-in-fraction)))

(defstruct psrf 
  "Potential scale reduction factor."
  r
  v
  w)

(defun calculate-psrf (accumulators &key (confidence 0.975d0))
  "Estimate the potential scale reduction factor.  Algorithm is from Gelman
and Rubin (1992), but the degrees of freedom correction is according to Brooks
and Gelman (1998)."
  ;; !!! should return the upper limit of the confidence interval as the
  ;; second value.  Since the F distribution is not implemented yet in
  ;; cl-random, this functionality is not available now.
  (declare (ignore confidence))
  (let+ ( ;; length and number of chains
         (m (length accumulators))
         (n (common accumulators :key #'tally))
         ;; means and variances for each
         (means (map1 #'mean accumulators))
         (variances (map1 #'variance accumulators))
         ;; calculate psrf
         ((&accessors (mu mean) (var-m variance)) (sweep 'sse means))
         (b (* n var-m))
         ((&accessors (w mean) (var-v variance)) (sweep 'sse variances))
         (var-b (/ (* 2 (expt b 2)) (1- m)))
         (var-w (/ var-v m))
         (1+1/m (1+ (/ m)))
         (n-1 (1- n))
         (V (/ (+ (* n-1 w) (* 1+1/m b)) n))
         (var-V (/ (+ (* (expt n-1 2) var-w)
                      (* (expt 1+1/m 2) var-b)
                      (* 2 1+1/m n-1 (/ n m)
                         (- (covariance-xy variances (eexpt means 2))
                            (* 2 mu (covariance-xy variances means)))))
                   (expt n 2)))
         (df (/ (* 2 (expt V 2)) var-V))
         (df-adj (/ (+ df 3) (1+ df)))
         ;; (b-df (1- m))
         ;; (w-df (/ (* 2 (expt w 2)) var-w))
         (R^2-fixed (/ n-1 n))
         (R^2-random (* (/ 1+1/m n) (/ b w))))
    (make-psrf :R (sqrt (* df-adj (+ R^2-fixed R^2-random)))
               :V V
               :W w)))

(defun calculate-psrf-ranges (n 
                              &key (divisions 20)
                                   (burn-in-fraction *default-burn-in-fraction*)
                                   (minimum-length 100))
  "Calculate ranges for PSRF.  Return as a list of (start . end) values.
Ranges narrower than MINIMUM-LENGTH are discarded."
  (iter
    (for division from 1 to divisions) 
    (let* ((end (ceiling (* division n) divisions))
           (start (floor (* end burn-in-fraction))))
      (when (<= (+ start minimum-length) end)
        (collect (cons start end))))))

(defclass mcmc-statistics ()
  ((model :accessor model :initarg :model)   
   (accumulators :accessor accumulators :initarg :accumulators
                 :documentation "Vector of accumulators for columns of
                 scalars.")
   (autocovariance-accumulators :initarg :autocovariance-accumulators
                                :accessor autocovariance-accumulators
                                :documentation "Vector of autocovariance
                                accumulators for each variable.")
   (sse-ranges :initarg :sse-ranges :accessor sse-ranges)
   (sse-accumulators :initarg :sse-accumulators :accessor sse-accumulators
                         :documentation "Vector of partial mean-sse
                         accumulators for each variable."))
  (:documentation "Statistics for the sample from a single MCMC chain."))

(defun mcmc-statistics (sample
                        &key (divisions 20) (minimum-length 100)
                             sse-ranges (lags 10)
                             (accumulator-generator #'mean-sse-accumulator)
                             (burn-in-fraction *default-burn-in-fraction*))
  "Helper function to calculate an MCMC-STATISTICS object from a sample."
  (let+ ((model (model (first* sample)))
         (n-parameters (layout-length (scalar-parameters-layout model)))
         (n-sample (length sample))
         (accumulators (generate-array n-parameters accumulator-generator))
         (autocovariance-accumulators
          (generate-array n-parameters (curry #'autocovariance-accumulator lags)))
         (sse-ranges (aif sse-ranges
                          it
                          (calculate-psrf-ranges 
                           n-sample
                           :divisions divisions
                           :minimum-length minimum-length
                           :burn-in-fraction burn-in-fraction)))
         (burn-in (calculate-burn-in n-sample burn-in-fraction))
         ((&values subranges index-lists)
          (subranges sse-ranges :shadow-ranges `((,burn-in . ,n-sample))))
         (sse-accumulators
          (combine
           (map 'vector
                (lambda+ ((start . end))
                  (let ((sse-accumulators
                         (generate-array n-parameters #'mean-sse-accumulator)))
                    (loop for sample-index from start below end do
                      (let+ (((&accessors-r/o scalar-parameters)
                              (aref sample sample-index)))
                        (iter
                          (for parameter :in-vector scalar-parameters
                               :with-index parameter-index)
                          (for sse-accumulator :in-vector sse-accumulators)
                          (add sse-accumulator parameter)
                          (when (<= burn-in sample-index)
                            (add (aref accumulators parameter-index)
                                 parameter)
                            (add (aref autocovariance-accumulators 
                                       parameter-index)
                                 parameter)))))
                    sse-accumulators))
                subranges)))
         (sse-accumulators
          (map1 (lambda (accumulators)
                  (iter
                    (for index-list :in-vector index-lists)
                    (let ((accumulators (sub accumulators
                                             (coerce index-list 'vector))))
                      (collect (pool* accumulators) :result-type vector))))
                (subarrays 1 (transpose sse-accumulators)))))
    ;; check burn-in
    (check-burn-in burn-in)
    ;; return results
    (make-instance 'mcmc-statistics
                   :model model
                   :accumulators accumulators
                   :autocovariance-accumulators autocovariance-accumulators
                   :sse-ranges sse-ranges
                   :sse-accumulators sse-accumulators)))

(defclass mcmc-summary ()
  ((model :accessor model :initarg :model)   
   (psrf :accessor psrf :initarg :psrf)
   (psrf-ranges :accessor psrf-ranges :initarg :psrf-ranges)
   (accumulators :accessor accumulators :initarg :accumulators)
   (mean-autocorrelations :accessor mean-autocorrelations
                          :initarg :mean-autocorrelations)))

(defun summarize-mcmc-statistics (mcmc-statistics)
  "Calculate summaries of column statistics."
  (let+ ((mcmc-statistics (coerce mcmc-statistics 'vector))
         (model (common-model mcmc-statistics))
         (sse-ranges (common mcmc-statistics :key #'sse-ranges :test #'equalp
                             :error "SSE ranges are not compatible."))
         ((&flet pool-chains (accessor transformation reduction)
            ;; Extract vector of statistics from each chain using ACCESSOR,
            ;; apply TRANSFORMATION, then pool them using REDUCTION
            (map1
             reduction
             (subarrays 1 
                        (map1 transformation
                              (transpose 
                               (combine (map1 accessor 
                                              mcmc-statistics))))))))
         ;; psrf
         (variance-accumulators
          (permute 
           (combine (map1 (compose #'combine #'sse-accumulators)
                          mcmc-statistics))
           '(1 2 0)))
         (psrf (map1 #'calculate-psrf
                     (subarrays 2 variance-accumulators)))
         ;; autocorrelations
         (mean-autocorrelations
          (pool-chains #'autocovariance-accumulators
                       #'autocorrelations
                       #'mean))
         ;; pooled accumulators
         (accumulators (pool-chains #'accumulators #'identity #'pool*)))
    (assert sse-ranges)
    (make-instance 'mcmc-summary
                   :model model :psrf psrf
                   :accumulators accumulators
                   :mean-autocorrelations mean-autocorrelations
                   :psrf-ranges sse-ranges)))

(defun pool-samples (samples &key (burn-in-fraction *default-burn-in-fraction*))
  "Pool MCMC samplers, discarding BURN-IN-FRACTION."
  ;; does not check model
  (stack* t :v
          (map 'list
               (lambda (sample)
                 (subseq sample
                         (calculate-burn-in (length sample) burn-in-fraction)))
               samples)))

;; (defun psrf-summary-quantiles (mcmc-statistics-summary
;;                                &key (quantiles #(0d0 0.025d0 0.25d0 0.5d0
;;                                                  0.75d0 0.975d0 1d0)))
;;   (let ((r (map1 #'psrf-r (psrf mcmc-statistics-summary)))
;;         (q (combine (map1 (quantile)))))
;;     )
;;   )

;;; - functions that pool samples should just verify that they point to the
;;;   same model and be done with it

;; (defclass mcmc-chains ()
;;   ((mcmc-class :accessor mcmc-class :initarg :mcmc-class :documentation
;;                "Class used for creating MCMC instances.")
;;    (initargs :accessor initargs :initarg :initargs :documentation
;;              "Initial arguments used for creating MCMC instances.")
;;    (parameters-ix :accessor parameters-ix :initarg :parameters-ix
;;                   :documentation "Index for the parameter vectors.")
;;    (chains :accessor chains :initarg :chains :documentation
;;            "Chains, always holding the current state.")
;;    (chain-results :accessor chain-results :initarg :chain-results
;;                   :type simple-vector
;;                   :documentation "Matrices holding the chain-results.")
;;    (burn-in :accessor burn-in :initarg :burn-in
;;             :documentation "Burn-in, used to discard start of the sequence
;;             before inference.")
;;    (pooled-parameters :accessor pooled-parameters :documentation
;;                       "Pooled parameters.")))

;; (defun run-mcmc-chains (m n mcmc-class initargs &key (burn-in (floor n 2))
;;                         (thin 1))
;;   "Run M MCMC chains, each of length N, with given class and initargs."
;;   (iter
;;     (with parameters-ix)
;;     (for chain-index :below m)
;;     (format t "Running chain ~A/~A~%" chain-index (1- m))
;;     (let ((mcmc (apply #'make-instance mcmc-class initargs)))
;;       (collecting (run-mcmc mcmc n :burn-in 0 :thin thin)
;;                   :result-type vector :into chain-results)
;;       (collecting mcmc :result-type vector :into chains)
;;       (when (first-iteration-p)
;;         (setf parameters-ix (parameters-ix mcmc)))
;;       (finally 
;;        (return
;;          (make-instance 'mcmc-chains
;;                         :chains chains
;;                         :chain-results chain-results
;;                         :initargs initargs
;;                         :parameters-ix parameters-ix
;;                         :mcmc-class mcmc-class
;;                         :burn-in burn-in))))))




;; (defun chains-psrf (mcmc-chains &key (divisions 20) (burn-in-fraction 2))
;;   "Calculate the potential scale reduction factor for "
;;   (let+ (((&slots-r/o chain-results) mcmc-chains)
;;          ((n k) (array-dimensions (aref chain-results 0)))
;;          (limits (iter
;;                    (for index :from 1 :to divisions)
;;                    (collecting (ceiling (* n index) divisions)
;;                                :result-type vector)))
;;          (psrf-matrix (make-array (list divisions k))))
;;     (dotimes (param-index k)
;;       (let ((sequences (map 'vector (lambda (chain) (sub chain t param-index))
;;                             chain-results)))
;;         (iter
;;           (for limit :in-vector limits :with-index limit-index)
;;           (let* ((start (floor limit burn-in-fraction))
;;                  (sequences (map 'vector (lambda (sequence)
;;                                            (subseq sequence start limit))
;;                                  sequences)))
;;             (setf (aref psrf-matrix limit-index param-index)
;;                   (psrf sequences))))))
;;     (values limits psrf-matrix)))

;; (defun calculate-pooled-parameters (mcmc-chains &key (start (burn-in mcmc-chains))
;;                                     (end (nrow (aref (chain-results mcmc-chains) 0))))
;;   "Combine MCMC chains into a single matrix, preserving column structure.
;; START and END mark the iterations used."
;;   (let+ ((chain-length (- end start))
;;          ((&slots-r/o chain-results) mcmc-chains)
;;          (m (length chain-results))
;;          (first-chain (aref chain-results 0))
;;          (pooled (make-array 
;;                         (list (* chain-length m) (ncol first-chain))
;;                         :element-type (array-element-type first-chain))))
;;     (iter
;;       (for chain :in-vector chain-results)
;;       (for end-row :from chain-length :by chain-length)
;;       (for start-row :previous end-row :initially 0)
;;       (setf (sub pooled (si start-row end-row) t)
;;             (sub chain (si start end) t)))
;;     pooled))

;; (defmethod slot-unbound (class (instance mcmc-chains)
;;                          (slot-name (eql 'pooled-parameters)))
;;   (setf (slot-value instance 'pooled-parameters)
;;         (calculate-pooled-parameters instance)))
