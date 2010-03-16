(in-package :cl-bayesian)

;;;  Counter for Metropolis (and Metropolis-Hastings) steps.
;;;
;;;  Counts the TOTAL and ACCEPTED number of steps.  To reset the
;;;  counter, simply create a new one.

(defclass acceptance-counter ()
  ((total :accessor total :initform 0)
   (accepted :accessor accepted :initform 0)))

(defun acceptance-ratio (counter)
  (coerce (/ (accepted counter) (total counter)) 'double-float))

(defun increment-counter (counter accepted-p)
  (when accepted-p
    (incf (accepted counter)))
  (incf (total counter))
  (values))

;;;  MCMC macros.
;;;
;;;  DEFINE-MCMC defines a class with the given parameter names, and
;;;  also methods for generic functions RESET-COUNTERS,
;;;  CURRENT-PARAMETERS (returns parameters for later processing) and
;;;  UPDATE, which performs an updating step on all variables.
;;;
;;;  The user needs to define (update-parameter class 'name) methods,
;;;  which return the updated parameter value.  Update calls these in
;;;  turn, and takes care of assigning this to the slots.
;;;
;;;  Slots can be atoms and vectors.  The difference is that atoms are
;;;  treated as an opaque object and updated as a block, while vectors
;;;  are updated elementwise.


(define-abstract-class mcmc ()
  ())

(defgeneric reset-counters (mcmc)
  (:documentation "reset all the counters in an MCMC object"))

(defgeneric update (mcmc)
  (:documentation "update parameters in an MCMC object"))

(defgeneric update-parameter (mcmc parameter)
  (:documentation "return the updated parameter in an MCMC object"))

(defgeneric update-parameter-in-vector (mcmc parameter index)
  (:documentation "return the updated parameter for index i in an MCMC object"))

(defgeneric current-parameters (mcmc)
  (:documentation "return the parameters"))

(defmacro define-mcmc (class-name direct-superclasses slots &rest options)
  "Example:
  (define-mcmc model ()
    ((x :parameter (atom :updater gibbs))
     (y :parameter (atom :updater metropolis))))"
  ;; NOTES: currently, metropolis vector updaters are sharing a counter
  (let (parameters                ; symbols
        vector-parameters         ; symbols
        counters                  ; slot name, original var name pairs
        propdists)                ; slot name, original var name pairs
    (labels ((process-parameter-specifier (name parameter-specifier)
               "Process parameter specifier."
               (bind (((type &key (updater :gibbs)
                             (counter 'counter counter-supplied-p)
                             (propdist 'propdist propdist-supplied-p))
                       (if (atom parameter-specifier)
                           (list parameter-specifier)
                           parameter-specifier)))
                 (case type
                   (atom)
                   (vector
                      (push name vector-parameters))
                   (otherwise (error "parameter type ~a not recognized" type)))
                 (push name parameters)
                 (case updater
                   ;; Gibbs: nothing needs to be done, just some sanity checks
                   (:gibbs 
                      (when (or counter-supplied-p propdist-supplied-p)
                        (error "a Gibbs sampler doesn't need a counter and/or ~
                                updater-parameters")))
                   ;; Metropolis
                   (:metropolis
                      (push (cons (make-symbol* name '- counter) name) counters)
                      (push (cons (make-symbol* name '- propdist) name) propdists))
                   (otherwise
                      (error "updater ~a not recognized" updater)))))
             (process-slot-specifier (slot-specifier)
               "Extract parameter definitions, return filtered slot
               specifier with MCMC-specific keyword pairs removed."
               (bind (((slot-name &rest options) slot-specifier)
                      (pairs (group options 2))
                      (parameter (find :parameter pairs :key #'first)))
                 (awhen (has-duplicates? pairs :key #'first)
                   (error "Key ~A occurs multiple times in slot specifier ~A."
                          (first pairs) slot-specifier))
                 (when parameter
                   (process-parameter-specifier slot-name (second parameter)))
                 (cons slot-name (mapcan (lambda (pair)
                                           (if (eq (first pair) :parameter)
                                               nil
                                               pair))
                                         pairs))))
             (generate-counter-slot (counter)
               "Generate the slot definition for a counter."
               (let ((name (car counter))
                     (documentation (format nil "counter for ~A" (cdr counter))))
               `(,name :accessor ,name :documentation ,documentation
                       :initform (make-instance 'acceptance-counter))))
             (generate-propdist-slot (propdist)
               "Generate the slot definition for a proposal distribution."
               (let ((name (car propdist))
                     (documentation (format nil "parameter(s) of the proposal ~
                                                 distribution for ~A" 
                                            (cdr propdist))))
               `(,name :accessor ,name :documentation ,documentation
                       :initarg ,(make-keyword name)))))
      (check-type class-name symbol)
      `(progn
         ;; class definition
         (defclass ,class-name (mcmc ,@direct-superclasses)
           ,(concatenate 'list
             (mapcar #'process-slot-specifier (reverse slots))
             (mapcar #'generate-counter-slot counters)
             (mapcar #'generate-propdist-slot propdists))
           ,@options)
         ;; reset
         (defmethod reset-counters ((mcmc ,class-name))
           ,@(mapcar (lambda (counter)
                       `(setf (,(car counter) mcmc)
                              (make-instance 'acceptance-counter)))
                     counters)
           (values))
         ;; update all variables
         (defmethod update ((mcmc ,class-name))
           (dolist (parameter ',parameters)
             (setf (slot-value mcmc parameter)
                   (update-parameter mcmc parameter)))
           (values))
         ;; updaters for vectors
         ,@(mapcar (lambda (name)
                     `(defmethod update-parameter ((mcmc ,class-name) (parameter (eql ',name)))
                        (with-slots (,name) mcmc
                          (dotimes (i (length ,name))
                            (setf (aref ,name i) (update-parameter-in-vector mcmc ',name i))))))
                   vector-parameters)))))


;;;;
;;;;  Utility functions for defining updaters.
;;;;

(defmacro define-updater ((class parameter &key (instance 'mcmc) (vector-index nil))
                          (&rest slots)  &body body)
  "Define an update-parameter (or update-parameter-in-vector, if
vector-index) method specialized to class and parameter.  The method will
be called with the given instance name.  Slots are expanded with bind
using :slots-read-only.  If vector-index, it will be used to index the vector."
  (check-type class symbol)
  (check-type parameter symbol)
  (check-type instance symbol)
  (check-type vector-index symbol)      ; nil is a symbol, too
  `(defmethod ,@(if vector-index
                    `(update-parameter-in-vector 
                        ((,instance ,class) (parameter (eql ',parameter))
                         ,vector-index))
                    `(update-parameter
                      ((,instance ,class) (parameter (eql ',parameter)))))
       (bind (,@(if slots
                    `(((:slots-read-only ,@slots) ,instance))
                    nil))
         ,@body)))

(defmacro define-metropolis-updater ((class parameter &key
                                            (instance 'mcmc)
                                            (vector-index nil)
                                            (counter (make-symbol* parameter
                                                                   '-counter))
                                            (propdist (make-symbol* parameter
                                                                    '-propdist)))
                                     (&rest slots) &body body)
  "Like define-updater, but with counter and proposal distribution
available with the given slot names (can be slot-name
or (variable-name slot-name)."
  `(define-updater (,class ,parameter :instance ,instance :vector-index ,vector-index) ,slots
     (bind (((:slots ,parameter ,counter ,propdist) ,instance))
       ,@body)))


(defun log-posterior-ratio (x xnext log-posterior/proposal)
  "Calculate the log posterior ratio by calling the
log-posterior/proposal function at x and xnext.  NIL is interpreted as
minus infinity, and evaluation is lazy."
  (let ((p-xnext (funcall log-posterior/proposal xnext)))
    (if p-xnext
        (let ((p-x (funcall log-posterior/proposal x)))
          (if p-x
              (- p-xnext p-x)
              (error "current point has zero likelihood: this should never happen")))
        nil)))

(defun metropolis-step* (x x-proposal l-p-ratio)
  "Return (values X-NEXT PROPOSAL-ACCEPTED-P).  X-NEXT is X or
X-PROPOSAL, based on L-P-RATIO (the log posterior-ratio)."
  (let ((accept-p (cond
                    ((null l-p-ratio) nil)
                    ((<= 0 l-p-ratio) t)
                    (t (< (random 1d0) (exp l-p-ratio))))))
    (if accept-p
        (values x-proposal t)
        (values x))))
  
(defun metropolis-step (x x-proposal l-p-ratio counter)
  "Perform a Metropolis(-Hastings) step, incrementing the counter if
necessary."
  (bind (((:values x-next accepted-p) (metropolis-step* x x-proposal l-p-ratio)))
    (when accepted-p
      (increment-counter counter accepted-p))
    x-next))

(defun run-mcmc (mcmc n &key (burn-in (max (floor n 10) 1000)))
  (dotimes (i burn-in)
    (update mcmc)) 
  (reset-counters mcmc)
  (let ((result (make-ra-matrix :double 0 0 :capacity n)))
    (dotimes (i n)
      (update mcmc)
      (add result (current-parameters mcmc)))
    result))
