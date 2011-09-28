;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

;;;; DLM specification

(defstruct (dlm-evolution1 (:constructor dlm-evolution1))
  "Evolution equation for DLM, univariate case.  state' ~ N(G*state+mu,W)."
  (G 1d0 :type double-float)
  (mu 0d0 :type double-float)
  (W 1d0 :type double-float))

(defstruct (dlm-observation1 (:constructor dlm-observation1))
  "Observation equation for DLM, univariate case.  data ~ N(F*state,V)."
  (F 1d0 :type double-float)
  (V 1d0 :type double-float))

(defstruct (dlm (:constructor make-dlm%))
  "Specification for a dynamic linear model."
  evolution+ observation+)

(defun dlm-length (dlm)
  "Return the length of a DLM."
  (length (dlm-observation+ dlm)))

(defun checked-dlm-length (dlm &rest vectors)
  "Return the length of a DLM, checking consistency."
  (let ((n (dlm-length dlm)))
    (loop for v in vectors do
      (check-type v vector)
      (assert (= n (length v))))
    n))

(defun make-dlm (evolution+ observation+ 
                 &key (length
                       (cond
                         ((vectorp evolution+) (1+ (length evolution+)))
                         ((vectorp observation+) (length observation+))
                         (t (error "Can't infer length when the other ~
                                    arguments are not vectors.")))))
  "Create a dynamic linear model.  Objects other than vectors are accepted and
will be recycled appropriately.  LENGTH should be specified when none of the
other arguments are vectors."
  (flet ((ensure-vector (object length)
           (if (vectorp object)
               (prog1 object
                 (assert (= length (length object))))
               (make-array length :initial-element object))))
    (make-dlm% :evolution+ (ensure-vector evolution+ (1- length))
               :observation+ (ensure-vector observation+ length))))

(defun sub-dlm (dlm start &optional (end (dlm-length dlm)))
  "Return a contiguous part of a DLM."
  (let+ (((&structure-r/o dlm- evolution+ observation+) dlm))
    (make-dlm% :evolution+ (subseq evolution+ start (1- end))
               :observation+ (subseq observation+ start end))))

;;;; single-step building blocks for recursive methods

(defgeneric dlm-step (mC evolution)
  (:documentation "Return the updated distribution.")
  (:method ((mc r-normal) (evolution dlm-evolution1))
    (let+ (((&accessors-r/o (m mean) (C variance)) mC)
           ((&structure-r/o dlm-evolution1- G mu W) evolution))
      (r-normal (+ (* G m) mu)
                (+ W (* (expt G 2) C))))))

(defgeneric dlm-filter (aR observation data)
  (:documentation "Given a prior, an observation (equation) and the
  corresponding data point, return the posterior.")
  (:method ((aR r-normal) (observation dlm-observation1) (data real))
    (let+ ((data (coerce data 'double-float))
           ((&structure-r/o dlm-observation1- F V) observation)
           ((&accessors-r/o (a mean) (R variance)) aR)
           (forecast (* F a))
           (F/V (/ F V))
           (C (/ (+ (/ R) (* F F/V))))
           (gain (* F/V C)))
      (r-normal (+ a (* gain (- data forecast))) C)))
  (:method (prior observation (data null))
    prior))

(defgeneric dlm-sample (mC evolution next-draw next-aR)
  (:documentation "Sample backward, returning a draw.")
  (:method ((mC r-normal) (evolution dlm-evolution1) (next-draw real)
            (next-aR r-normal))
    (let+ (((&structure-r/o dlm-evolution1- G W) evolution)
           ((&accessors-r/o (m mean) (C variance)) mC)
           ((&accessors-r/o (next-a mean)) next-aR)
           (G/W (/ G W))
           (H (/ (+ (/ C) (* G G/W))))
           (B (* G/W H)))
      (draw (r-normal (+ m (* B (- next-draw next-a))) H)))))

;;;; user interface
;;; 
;;; The convention for argument order should be
;;;;  aR evolution+ observation+ data+ theta+ 

(defun dlm-forward-filtering (initial-distribution dlm data+)
  "Forward filtering for univariate dynamic linear models."
  (let+ (((&structure-r/o dlm- evolution+ observation+) dlm)
         (n (checked-dlm-length dlm data+))
         (distribution initial-distribution)
         (aR+ (make-array n))
         (mC+ (make-array n)))
    (iter
      (for observation :in-vector observation+ :with-index index)
      (for data :in-vector data+)
      (when (plusp index)
        (setf distribution (dlm-step distribution 
                                     (aref evolution+ (1- index)))))
      (setf (aref aR+ index) distribution
            distribution (dlm-filter distribution observation data)
            (aref mC+ index) distribution))
    (values mC+ aR+)))

(defun dlm-backward-sampling (dlm mC+ aR+)
  "Backward sampling.  Requires the output of dlm-forward-filtering and the
parameters.  Returns a vector of draws.  For univariate DLMs."
  (let+ (((&structure-r/o dlm- evolution+) dlm)
         (n (checked-dlm-length dlm mC+ aR+))
         (theta+ (make-array n))
         (last (1- n))
         (theta (draw (aref mC+ last))))
    (setf (aref theta+ last) theta)
    (iter
      (for evolution :in-vector evolution+ :from (1- last) :downto 0
                     :with-index index)
      (setf theta (dlm-sample (aref mC+ index) evolution theta
                              (aref aR+ (1+ index)))
            (aref theta+ index) theta))
    theta+))

(defun dlm-ff-bs (initial-distribution dlm data+)
  "Forward filtering and backward sampling.  Return (values theta+ m+
  C-inverse+ a+)."
  (let+ (((&values mC+ aR+)
          (dlm-forward-filtering initial-distribution dlm data+)))
    (values (dlm-backward-sampling dlm mc+ ar+) mC+ aR+)))

(defgeneric dlm-evolution-error (evolution state next-state)
  (:documentation "Calculate the evolution error (incudes MU).")
  (:method ((evolution dlm-evolution1) (state real) (next-state real))
    (- next-state (* (dlm-evolution1-G evolution) state))))

(defgeneric dlm-observation-error (observation state data)
  (:documentation "Calculate the observation error.")
  (:method ((observation dlm-observation1) (state real) (data real))
    (- data (* (dlm-observation1-F observation) state)))
  (:method (observation state (data null))
    nil))

(defun dlm-errors (dlm data+ theta+)
  "Return vectors of the errors of the state equation (omega, with the mean mu
which is *not* subtracted) and the observation equation (nu) as two values."
  (let+ (((&structure-r/o dlm- evolution+ observation+) dlm)
         (n (checked-dlm-length dlm data+ theta+))
         (omega+ (make-array (1- n)))
         (nu+ (make-array n)))
    (iter
      (for theta :in-vector theta+ :with-index index)
      (for theta-p :previous theta)
      (for observation :in-vector observation+)
      (for data :in-vector data+)
      (unless (zerop index)
        (let ((index-p (1- index)))
          (setf (aref omega+ index-p)
                (dlm-evolution-error (aref evolution+ index-p)
                                     theta-p theta))))
      (setf (aref nu+ index) (dlm-observation-error observation theta data)))
    (values omega+ nu+)))

(defmethod draw ((evolution dlm-evolution1) &key state)
  (let+ (((&structure-r/o dlm-evolution1- G mu W) evolution))
    (draw (r-normal (+ mu (* G state)) W))))

(defmethod draw ((observation dlm-observation1) &key state)
  (let+ (((&structure-r/o dlm-observation1- F V) observation))
    (draw (r-normal (+ (* F state)) V))))

(defun dlm-simulate (aR dlm)
  "Generate a sample for a DLM with given parametes, drawing the first state
from N(a,R).  Return (values state+ data+).  For univariate DLMs."
  (let+ (((&structure-r/o dlm- evolution+ observation+) dlm)
         (n (dlm-length dlm))
         (state (draw aR))
         (state+ (make-array n :element-type 'double-float))
         (data+ (make-array n :element-type 'double-float)))
    (dotimes (index n)
      (unless (zerop index)
        (setf state (draw (aref evolution+ (1- index)) :state state)))
      (setf (aref state+ index) state
            (aref data+ index) (draw (aref observation+ index) :state state)))
    (values state+ data+)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; multivariate dlm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; (defstruct uddu
;;   "Decomposition UD^2U^T where D is a diagonal and U is a unitary matrix."
;;   u d)

;; (defun uddu (a)
;;   "UDDU decomposition of hermitian matrix."
;;   (check-type a hermitian-matrix)
;;   (let+ (((&structure-r/o spectral-factorization- z w)
;;           (spectral-factorization a)))
;;     (make-uddu :u z :d (esqrt w))))

;; (defmethod left-square-root ((a uddu))
;;   (mm (uddu-u a) (uddu-d a)))

;; (defmethod as-matrix ((a uddu))
;;   (mm (left-square-root a) t))

;; (defmethod as-array ((a uddu) &key copy?)
;;   (declare (ignore copy?))
;;   (as-array (as-matrix a)))

;; (defmethod invert ((a uddu) &key)
;;   (let+ (((&structure-r/o uddu- u d) a))
;;     (make-uddu :u u :d (invert d))))

;; (defun uddu-left-svd (l)
;;   "Helper function to calculate the UDDU decomposition of LL^T, where L is not
;; necessarily square."
;;   (let+ (((&structure-r/o svd- u d) (svd l :thin)))
;;     (make-uddu :u u :d d)))

;; (defun uddu-update (a h &optional factorize?)
;;   "Return the UDDU decomposition of A+H, where H is a Hermitian matrix (or
;; anything else with a representation that yields a left square root)."
;;   (if factorize?
;;       ;; method which first factorizes out U
;;       (let+ (((&structure-r/o uddu- u d) a)
;;              ((&structure-r/o svd- (u+ u) (d+ d))
;;               (svd (stack 'double-float :h
;;                           d (mm (transpose u) (left-square-root h))) :thin)))
;;         (make-uddu :u (mm u u+) :d d+))
;;       ;; direct method
;;       (uddu-left-svd (stack 'double-float :h
;;                             (left-square-root a) (left-square-root h)))))

;; (defun uddu-multiply-update (a g h)
;;   "Return the UDDU decomposition of GAG^T+H, where A is an UDDU decomposition
;;   and H is a hermitian matrix (in any representation that can yield a left
;;   square root)."
;;   (uddu-left-svd (stack 'double-float :h
;;                         (mm g (left-square-root a)) (left-square-root h))))

;; (defstruct (dlm-parameters (:constructor make-dlm-parameters%))
;;   "We follow the notation of the Harrison and West (1997) book:

;;   theta_t = G_t theta_{t-1} + omega_t, where omega_t ~ N(mu_t, W_t)

;;   Y_t = F theta_t + nu_t, where nu_t ~ N(0,V_t)

;;   W and V can be factorizations."
;;     G mu W F V)

;; (defun make-dlm-parameters (&key G mu W F V)
;;   "Create a DLM parameters structure, checking dimensions.  Scalars are
;; converted to appropriate 1x1 matrices, and F may be a vector, interpeted as a
;; row matrix."
;;   (let+ (((&flet ensure-variance (m)
;;             (if (numberp m)
;;                 (clo :diagonal m)
;;                 m)))
;;          (G (ensure-matrix G))
;;          (mu (ensure-vector mu))
;;          (W (ensure-variance W))
;;          (F (ensure-matrix F :row))
;;          (V (ensure-variance V)))
;;     (assert (and (= (nrow G) (length mu) (nrow W) (ncol F))
;;                  (= (nrow F) (nrow V))
;;                  (square? G) (square? W) (square? V)))
;;     (make-dlm-parameters% :g G :mu mu :w W :F F :V V)))

;; (defun dlm-filter (a R Y parameters)
;;   "Given the prior N(a,R), the observation Y and the DLM parameters at that
;;   point, calculate the posterior N(m,C).  Return (values m C-inverse)."
;;   (let+ (((&structure-r/o dlm-parameters- F V) parameters)
;;          (F-transpose (transpose F))
;;          (forecast (mm F a)) ; mean forecast
;;          (V-inverse (invert V))
;;          (C-inverse (uddu-update (invert R)
;;                                  (xx (mm F-transpose
;;                                          (left-square-root V-inverse)))))
;;          (gain (mmm (invert C-inverse) F-transpose V-inverse)))
;;     (values (e+ a (mm gain (e- Y forecast)))
;;             C-inverse)))

;; (defun dlm-step (m C-inverse parameters)
;;   "Given the posterior N(m,C), calculate the prior N(a,R) for the next step.
;; Return (values a R).  This step should be called each point in time, before
;; incorporating observations (if any)."
;;   (let+ (((&structure-r/o dlm-parameters- G mu W) parameters))
;;     (values (e+ (mm G m) mu)
;;             (uddu-multiply-update (invert C-inverse) G W))))

;; (defun dlm-sample (m C-inverse next-theta next-a next-parameters)
;;   "Sample backward for DLM.  Return a draw."
;;   (let+ (((&structure-r/o dlm-parameters- G W) next-parameters)
;;          (H-inverse (uddu-update C-inverse 
;;                                  (xx (mm (transpose G)
;;                                          (left-square-root (invert W))))))
;;          (H (invert H-inverse))
;;          (B (mmm H (transpose G) (invert W))))
;;     (draw (r-multivariate-normal (e+ m (mm B (e- next-theta next-a)))
;;                                  H))))

;; (defun dlm-forward-filtering (a0 R0 y+ parameters+)
;;   "Forward filtering for dynamic linear models.
;; Prior on the state is N(a,R), Y+ is a vector of observations, PARAMETERS is a
;; vector of DLM-PARAMETERs.  Return the following values:

;;  - m+, a vector of means
;;  - C-inverse+, a vector of UDDU decompositions of C^{-1}
;;  - a+, the predicted means"
;;   (let* ((n (length y+))
;;          (m+ (make-array n))
;;          (C-inverse+ (make-array n))
;;          (a+ (make-array n)))
;;     ;; forward filtering
;;     (iter
;;       (for parameters :in-vector parameters+ :with-index index)
;;       (let+ (((&values a R) (if (zerop index)
;;                                 (values a0 R0)
;;                                 (let ((m (aref m+ (1- index)))
;;                                       (C-inverse (aref C-inverse+
;;                                                        (1- index))))
;;                                   (dlm-step m C-inverse parameters)))))
;;         (setf (values (aref m+ index) (aref C-inverse+ index))
;;               (aif (aref y+ index)
;;                    (dlm-filter a R it parameters)
;;                    (values a (invert R)))
;;               (aref a+ index) a)))
;;     (values m+ C-inverse+ a+)))

;; (defun dlm-backward-sampling (m+ a+ C-inverse+ parameters+)
;;   "Backward sampling.  Requires the output of dlm-forward-filtering and the
;; parameters.  Returns a vector of draws."
;;   (let* ((n (length m+))
;;          (theta+ (make-array n))
;;          (last (1- n)))
;;     (setf (aref theta+ last)
;;           (draw  (r-multivariate-normal (aref m+ last)
;;                                         (invert (aref C-inverse+ last)))))
;;     (iter
;;       (for parameters :in-vector parameters+ :downto 1 :with-index index)
;;       (setf (aref theta+ (1- index))
;;             (dlm-sample (aref m+ (1- index)) (aref C-inverse+ (1- index))
;;                         (aref theta+ index) (aref a+ index) parameters)))
;;     theta+))

;; (defun dlm-ff-bs (a0 R0 y+ parameters+)
;;   "Forward filtering and backward sampling.  Return (values theta+ m+
;;   C-inverse+ a+)."
;;   (let+ (((&values m+ C-inverse+ a+)
;;           (dlm-forward-filtering a0 R0 y+ parameters+)))
;;     (values (dlm-backward-sampling m+ a+ c-inverse+ parameters+)
;;             m+ C-inverse+ a+)))

;; (defun dlm-errors (theta+ y+ parameters+)
;;   "Return vectors of the errors of the state equation (omega, mean included)
;; and the observation equation (nu) as two values.  Note that the first element
;; of OMEGA is NIL as it is not identified."
;;   (let+ ((n (common-length theta+ y+ parameters+))
;;          ((&assert n))
;;          (omega (make-array n :initial-element nil))
;;          (nu (make-array n :initial-element nil)))
;;     (iter
;;       (for theta :in-vector theta+ :with-index index)
;;       (for theta-p :previous theta)
;;       (for parameters :in-vector parameters+)
;;       (for y :in-vector y+)
;;       (let+ (((&structure-r/o dlm-parameters- G F) parameters))
;;         (unless (zerop index)
;;           (setf (aref omega index) (e- theta (mm G theta-p))))
;;         (when y
;;           (setf (aref nu index) (e- y (mm F theta))))))
;;     (values omega nu)))

;; (defun dlm-simulate (a R parameters+)
;;   "Generate a sample for a DLM with given parametes, drawing the first state
;; from N(a,R).  Return (values theta+ y+)."
;;   (let* ((n (length parameters+))
;;          (theta+ (make-array n))
;;          (y+ (make-array n)))
;;     (dotimes (index n)
;;       (let+ (((&structure-r/o dlm-parameters- G W mu F V)
;;               (aref parameters+ index))
;;              (theta (if (zerop index)
;;                         (draw (r-multivariate-normal a R))
;;                         (e+ (mm G (aref theta+ (1- index)))
;;                             (draw (r-multivariate-normal mu W))))))
;;         (setf (aref y+ index)
;;               (draw (r-multivariate-normal (mm F theta) V))
;;               (aref theta+ index) theta)))
;;     (values theta+ y+)))
