;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

(defclass mcmc-chains ()
  ((mcmc-class :accessor mcmc-class :initarg :mcmc-class :documentation
               "Class used for creating MCMC instances.")
   (initargs :accessor initargs :initarg :initargs :documentation
             "Initial arguments used for creating MCMC instances.")
   (parameters-ix :accessor parameters-ix :initarg :parameters-ix
                  :documentation "Index for the parameter vectors.")
   (chains :accessor chains :initarg :chains :documentation
           "Chains, always holding the current state.")
   (chain-results :accessor chain-results :initarg :chain-results
                  :type simple-vector
                  :documentation "Matrices holding the chain-results.")
   (burn-in :accessor burn-in :initarg :burn-in
            :documentation "Burn-in, used to discard start of the sequence
            before inference.")
   (pooled-parameters :accessor pooled-parameters :documentation
                      "Pooled parameters.")))

(defun run-mcmc-chains (m n mcmc-class initargs &key (burn-in (floor n 2))
                        (thin 1))
  "Run M MCMC chains, each of length N, with given class and initargs."
  (iter
    (with parameters-ix)
    (for chain-index :below m)
    (format t "Running chain ~A/~A~%" chain-index (1- m))
    (let ((mcmc (apply #'make-instance mcmc-class initargs)))
      (collecting (run-mcmc mcmc n :burn-in 0 :thin thin)
                  :result-type vector :into chain-results)
      (collecting mcmc :result-type vector :into chains)
      (when (first-iteration-p)
        (setf parameters-ix (parameters-ix mcmc)))
      (finally 
       (return
         (make-instance 'mcmc-chains
                        :chains chains
                        :chain-results chain-results
                        :initargs initargs
                        :parameters-ix parameters-ix
                        :mcmc-class mcmc-class
                        :burn-in burn-in))))))

(defun psrf (sequences &key (confidence 0.975d0) skip-length-check?)
  "Estimate the potential scale reduction factor.  Algorithm is from Gelman and
Rubin (1992), but the degrees of freedom correction is according to Brooks and
Gelman (1998)."
  ;; !!! should return the upper limit of the confidence interval as the second
  ;; value.  Since the F distribution is not implemented yet, this functionality
  ;; is not available now.
  (declare (ignore confidence))
  (let ((m (length sequences))
        (n (length (aref sequences 0))))
    (unless skip-length-check?
      (assert (every (lambda (sequence) (= (length sequence) n))
                     (subseq sequences 1))))
    (iter
      (for sequence :in-vector sequences)
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

(defun chains-psrf (mcmc-chains &key (divisions 20) (burn-in-fraction 2))
  "Calculate the potential scale reduction factor for "
  (bind (((:slots-r/o chain-results) mcmc-chains)
         ((n k) (array-dimensions (aref chain-results 0)))
         (limits (iter
                   (for index :from 1 :to divisions)
                   (collecting (ceiling (* n index) divisions)
                               :result-type vector)))
         (psrf-matrix (make-array (list divisions k))))
    (dotimes (param-index k)
      (let ((sequences (map 'vector (lambda (chain) (sub chain t param-index))
                            chain-results)))
        (iter
          (for limit :in-vector limits :with-index limit-index)
          (let* ((start (floor limit burn-in-fraction))
                 (sequences (map 'vector (lambda (sequence)
                                           (subseq sequence start limit))
                                 sequences)))
            (setf (aref psrf-matrix limit-index param-index)
                  (psrf sequences))))))
    (values limits psrf-matrix)))

(defun calculate-pooled-parameters (mcmc-chains &key (start (burn-in mcmc-chains))
                                    (end (nrow (aref (chain-results mcmc-chains) 0))))
  "Combine MCMC chains into a single matrix, preserving column structure.
START and END mark the iterations used."
  (bind ((chain-length (- end start))
         ((:slots-r/o chain-results) mcmc-chains)
         (m (length chain-results))
         (first-chain (aref chain-results 0))
         (pooled (make-array 
                        (list (* chain-length m) (ncol first-chain))
                        :element-type (array-element-type first-chain))))
    (iter
      (for chain :in-vector chain-results)
      (for end-row :from chain-length :by chain-length)
      (for start-row :previous end-row :initially 0)
      (setf (sub pooled (si start-row end-row) t)
            (sub chain (si start end) t)))
    pooled))

(defmethod slot-unbound (class (instance mcmc-chains)
                         (slot-name (eql 'pooled-parameters)))
  (setf (slot-value instance 'pooled-parameters)
        (calculate-pooled-parameters instance)))
