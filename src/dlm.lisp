;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

(defstruct uddu
  "Decomposition UD^2U^T where D is a diagonal and U is a unitary matrix."
  u d)

(defun uddu (a)
  "UDDU decomposition of hermitian matrix."
  (check-type a hermitian-matrix)
  (let+ (((&structure-r/o spectral-factorization- z w)
          (spectral-factorization a)))
    (make-uddu :u z :d (esqrt w))))

(defmethod left-square-root ((a uddu))
  (mm (uddu-u a) (uddu-d a)))

(defmethod as-array ((a uddu) &key hermitian?)
  (lla::maybe-hermitian (mm (left-square-root a) t) hermitian?))

(defmethod invert ((a uddu) &key)
  (let+ (((&structure-r/o uddu- u d) a))
    (make-uddu :u u :d (invert d))))

(defun uddu-left-svd (l)
  "Helper function to calculate the UDDU decomposition of LL^T, where L is not
necessarily square."
  (let+ (((&structure-r/o svd- u d) (svd l :thin)))
    (make-uddu :u u :d d)))

(defun uddu-update (a h &optional factorize?)
  "Return the UDDU decomposition of A+H, where H is a Hermitian matrix (or
anything else with a representation that yields a left square root)."
  (if factorize?
      ;; method which first factorizes out U
      (let+ (((&structure-r/o uddu- u d) a)
             ((&structure-r/o svd- (u+ u) (d+ d))
              (svd (stack 'double-float :h
                          d (mm (transpose u) (left-square-root h))) :thin)))
        (make-uddu :u (mm u u+) :d d+))
      ;; direct method
      (uddu-left-svd (stack 'double-float :h
                            (left-square-root a) (left-square-root h)))))

(defun uddu-multiply-update (a g h)
  "Return the UDDU decomposition of GAG^T+H, where A is an UDDU decomposition
  and H is a hermitian matrix (in any representation that can yield a left
  square root)."
  (uddu-left-svd (stack 'double-float :h
                        (mm g (left-square-root a)) (left-square-root h))))

(defstruct (dlm-parameters (:constructor make-dlm-parameters%))
  "We follow the notation of the Harrison and West (1997) book:

  theta_t = G_t theta_{t-1} + omega_t, where omega_t ~ N(mu_t, W_t)

  Y_t = F theta_t + nu_t, where nu_t ~ N(0,V_t)

  W and V can be factorizations."
    G mu W F V)

(defun make-dlm-parameters (&key G mu W F V)
  "Create a DLM parameters structure, checking dimensions.  Scalars are
converted to appropriate 1x1 matrices, and F may be a vector, interpeted as a
row matrix."
  (let+ (((&flet ensure-variance (m)
            (if (numberp m)
                (clo :diagonal m)
                m)))
         (G (ensure-matrix G))
         (mu (ensure-vector mu))
         (W (ensure-variance W))
         (F (ensure-matrix F :row))
         (V (ensure-variance V)))
    (assert (and (= (nrow G) (length mu) (nrow W) (ncol F))
                 (= (nrow F) (nrow V))
                 (square? G) (square? W) (square? V)))
    (make-dlm-parameters% :g G :mu mu :w W :F F :V V)))

(defun dlm-filter (a R Y parameters)
  "Given the prior N(a,R), the observation Y and the DLM parameters at that
  point, calculate the posterior N(m,C).  Return (values m C-inverse)."
  (let+ (((&structure-r/o dlm-parameters- F V) parameters)
         (F-transpose (transpose F))
         (forecast (mm F a)) ; mean forecast
         (V-inverse (invert V))
         (C-inverse (uddu-update (invert R)
                                 (xx (mm F-transpose
                                         (left-square-root V-inverse)))))
         (gain (mmm (invert C-inverse) F-transpose V-inverse)))
    (values (e+ a (mm gain (e- Y forecast)))
            C-inverse)))

(defun dlm-step (m C-inverse parameters)
  "Given the posterior N(m,C), calculate the prior N(a,R) for the next step.
Return (values a R).  This step should be called each point in time, before
incorporating observations (if any)."
  (let+ (((&structure-r/o dlm-parameters- G mu W) parameters))
    (values (e+ (mm G m) mu)
            (uddu-multiply-update (invert C-inverse) G W))))

(defun dlm-sample (m C-inverse next-theta next-a next-parameters)
  "Sample backward for DLM.  Return a draw."
  (declare (optimize debug))
  (let+ (((&structure-r/o dlm-parameters- G W) next-parameters)
         (H-inverse (uddu-update C-inverse 
                                 (xx (mm (transpose G)
                                         (left-square-root (invert W))))))
         (H (invert H-inverse))
         (B (mmm H (transpose G) (invert W))))
    (draw (r-multivariate-normal (e+ m (mm B (e- next-theta next-a)))
                                 H))))

(defun dlm-forward-filtering (a0 R0 y+ parameters+)
  "Forward filtering for dynamic linear models.
Prior on the state is N(a,R), Y+ is a vector of observations, PARAMETERS is a
vector of DLM-PARAMETERs.  Return the following values:

 - m+, a vector of means
 - C-inverse+, a vector of UDDU decompositions of C^{-1}
 - a+, the predicted means"
  (let* ((n (length y+))
         (m+ (make-array n))
         (C-inverse+ (make-array n))
         (a+ (make-array n)))
    ;; forward filtering
    (iter
      (for parameters :in-vector parameters+ :with-index index)
      (let+ (((&values a R) (if (zerop index)
                                (values a0 R0)
                                (let ((m (aref m+ (1- index)))
                                      (C-inverse (aref C-inverse+
                                                       (1- index))))
                                  (dlm-step m C-inverse parameters)))))
        (setf (values (aref m+ index) (aref C-inverse+ index))
              (aif (aref y+ index)
                   (dlm-filter a R it parameters)
                   (values a (invert R)))
              (aref a+ index) a)))
    (values m+ C-inverse+ a+)))

(defun dlm-backward-sampling (m+ a+ C-inverse+ parameters+)
  "Backward sampling.  Requires the output of dlm-forward-filtering and the
parameters.  Returns a vector of draws."
  (let* ((n (length m+))
         (theta+ (make-array n))
         (last (1- n)))
    (setf (aref theta+ last)
          (draw  (r-multivariate-normal (aref m+ last)
                                        (invert (aref C-inverse+ last)))))
    (iter
      (for parameters :in-vector parameters+ :downto 1 :with-index index)
      (setf (aref theta+ (1- index))
            (dlm-sample (aref m+ (1- index)) (aref C-inverse+ (1- index))
                        (aref theta+ index) (aref a+ index) parameters)))
    theta+))

(defun dlm-ff-bs (a0 R0 y+ parameters+)
  "Forward filtering and backward sampling.  Return (values theta+ m+
  C-inverse+ a+)."
  (let+ (((&values m+ C-inverse+ a+)
          (dlm-forward-filtering a0 R0 y+ parameters+)))
    (values (dlm-backward-sampling m+ a+ c-inverse+ parameters+)
            m+ C-inverse+ a+)))

(defun dlm-errors (theta+ y+ parameters+)
  "Return vectors of the errors of the state equation (omega, mean included)
and the observation equation (nu) as two values.  Note that the first element
of OMEGA is NIL as it is not identified."
  (let* ((n (common-length theta+ y+ parameters+))
         (omega (make-array n :initial-element nil))
         (nu (make-array n)))
    (iter
      (for theta :in-vector theta+ :with-index index)
      (for theta-p :previous theta)
      (for parameters :in-vector parameters+)
      (for y :in-vector y+)
      (let+ (((&structure-r/o dlm-parameters- G F) parameters))
        (unless (zerop index)
          (setf (aref omega index) (e- theta (mm G theta-p))))
        (setf (aref nu index) (e- y (mm F theta)))))
    (values omega nu)))
