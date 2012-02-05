;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

(defgeneric pack-slots (mapping target object)
  (:documentation ""))

(defgeneric unpack-slots (mapping source object)
  (:argument-precedence-order object source mapping)
  (:documentation "")
  (:method (mapping source (class standard-class))
    (unpack-slots mapping source (make-instance class)))
  (:method (mapping source (symbol symbol))
    (unpack-slots mapping source (make-instance symbol))))

(defstruct (mosaic (:constructor make-mosaic%))
  "A mosaic is a mapping that packs arrays of any dimension into a flat
vector (the size of which can be obtained with MOSAIC-SIZE)."
  (keys nil :type vector :read-only t)
  (table nil :type hash-table :read-only t)
  (size nil :type fixnum :read-only t))

(defun make-mosaic (keys-and-dimensions)
  "Create a mosaic from a sequence of (cons key dimensions)."
  (let* ((table (make-hash-table :test #'equal))
         (offset 0)
         (keys (map 'vector
                    (lambda (key-and-dimensions)
                      (let+ (((key . dimensions)
                              (ensure-list key-and-dimensions))
                             ((&values nil present?) (gethash key table)))
                        (assert (not present?) () "Duplicate key ~A." key)
                        (prog1 key
                          (setf (gethash key table) (cons offset dimensions))
                          (incf offset (product dimensions)))))
                    keys-and-dimensions)))
    (make-mosaic% :keys keys :table table :size offset)))

(defun template-mosaic (keys-and-objects)
  "Create a mosaic using objects as a template."
  (make-mosaic (map 'vector (lambda+ ((key . object))
                              (cons key
                                    (if (arrayp object)
                                        (array-dimensions object)
                                        nil)))
                    keys-and-objects)))

(defmacro template-mosaic-symbols (&rest keys)
  "Convenience macro for templating a mosaic on keys or (key variable)."
  `(template-mosaic
    (list ,@(mapcar (lambda (key)
                      (let+ (((key &optional (form key)) (ensure-list key)))
                        `(cons ',key ,form)))
                    keys))))

(defun mosaic-location (mosaic key &optional subscripts)
  "Return (CONS OFFSET DIMENSIONS).  Read-only, consequences are undefined if
modified."
  (let+ (((&values value present?) (gethash key (mosaic-table mosaic))))
    (assert present? () "Key ~A not found." key)
    (if subscripts
        (let+ (((&values offset dimensions)
                (subarray-location (cdr value) subscripts)))
          (cons (+ (car value) offset) dimensions))
        value)))

(defmethod print-object ((mosaic mosaic) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (mosaic stream :type t)
        (format stream "size: ~A" (mosaic-size mosaic))
        (map nil (lambda (key)
                   (let+ (((offset . dimensions)
                           (mosaic-location mosaic key)))
                     (format stream "~&~4T~A ~:A [~A]" key dimensions offset)))
             (mosaic-keys mosaic)))))

(defun mosaic-displace-vector (mosaic key vector)
  (let+ (((offset . dimensions) (mosaic-location mosaic key)))
    (displace-array vector dimensions offset)))

(defmethod pack-slots ((mosaic mosaic) (vector vector)
                       (object standard-object))
  (assert (length= vector (mosaic-size mosaic)))
  (map nil (lambda (key)
             (let+ (((offset . dimension) (mosaic-location mosaic key))
                    (value (slot-value object key)))
               (if dimension
                   (progn
                     (assert (equal dimension (array-dimensions value)))
                     (replace vector (flatten-array value) :start1 offset))
                   (setf (aref vector offset) value))))
       (mosaic-keys mosaic))
  object)

(defun mosaic-unpack-key (mosaic vector key)
  (let+ (((offset . dimension) (mosaic-location mosaic key)))
    (if dimension
        (clnu:maybe-copy-array
         (displace-array vector dimension offset) nil)
        (aref vector offset))))

(defmethod unpack-slots :before ((mosaic mosaic) (vector vector) object)
  (assert (length= vector (mosaic-size mosaic))))

(defmethod unpack-slots ((mosaic mosaic) (vector vector)
                         (object standard-object))
  (map nil (lambda (key)
             (setf (slot-value object key)
                   (mosaic-unpack-key mosaic vector key)))
       (mosaic-keys mosaic))
  object)

(defmethod unpack-slots ((mosaic mosaic) (vector vector)
                         (object (eql :alist)))
  (loop for key across (mosaic-keys mosaic)
        collect (cons key (mosaic-unpack-key mosaic vector key))))

(defgeneric mosaic (object)
  (:documentation "Return mosaic of object."))

(defclass mosaic-with-elements (simple-print-object-mixin)
  ((mosaic :initarg :mosaic :type mosaic :reader mosaic)
   (elements :initarg :elements :type matrix :reader elements))
  (:documentation  "Base structure for a mosaic and associated elements."))

(defclass mosaic-matrix (mosaic-with-elements)
  ()
  (:documentation "A matrix with a mosaic indexing the columns."))

(defmethod nrow ((matrix mosaic-matrix))
  (nrow (elements matrix)))

(defun make-mosaic-matrix (mosaic nrow-or-matrix
                           &rest make-array-arguments
                           &key (element-type t element-type?)
                                (initial-element nil initial-element?))
  "Make a mosaic matrix.  Keyword arguments are passed on to make-array.  When
the second argument is a matrix, it is used (without copying), but the
consistency with mosaic is checked."
  (declare (ignorable initial-element))
  (let ((size (mosaic-size mosaic)))
    (make-instance 'mosaic-matrix
                   :mosaic mosaic
                   :elements (aetypecase nrow-or-matrix
                               (array-length
                                (apply #'make-array
                                       (list it size)
                                       make-array-arguments))
                               (matrix
                                (assert (or (not element-type?)
                                            (equalp (upgraded-array-element-type element-type)
                                                    (array-element-type nrow-or-matrix))))
                                (assert (not initial-element?))
                                it)))))

(defclass mosaic-vector (mosaic-with-elements)
  ((elements :type vector)))

(defun make-mosaic-vector (mosaic &rest make-array-arguments
                                  &key (element-type t) initial-element
                                       initial-contents)
  "Make a mosaic vector.  Keyword arguments are passed on to make-array."
  (declare (ignorable element-type initial-element initial-contents))
  (make-instance 'mosaic-vector
                 :mosaic mosaic
                 ;; FIXME: use semantics like make-mosaic-matrix
                 :elements (apply #'make-array (mosaic-size mosaic)
                                  make-array-arguments)))


(defmethod pack-slots ((mosaic-matrix mosaic-matrix) (row fixnum) object)
  (let+ (((&slots-r/o mosaic elements) mosaic-matrix))
    (pack-slots mosaic (subarray elements row) object)))

(defmethod unpack-slots ((mosaic-matrix mosaic-matrix) (row fixnum) object)
  (let+ (((&slots-r/o mosaic elements) mosaic-matrix))
    (unpack-slots mosaic (subarray elements row) object)))

(defmethod sub ((mosaic-matrix mosaic-matrix) &rest selections)
  (let+ (((row-selection key-selection &rest subscripts) selections)
         ((&slots-r/o mosaic elements) mosaic-matrix)
         ((nrow ncol) (array-dimensions elements))
         (row-selection (sub-resolve-selection row-selection nrow t)))
    (if (eq key-selection t)
        (let ((elements (sub elements row-selection t)))
          (if (vectorp elements)
              (make-instance 'mosaic-vector :mosaic mosaic :elements elements)
              (make-instance 'mosaic-matrix :mosaic mosaic :elements elements)))
        (let+ (((offset . dimensions)
                (mosaic-location mosaic key-selection subscripts))
               ((&flet extract (row-index)
                  (if dimensions
                      (displace-array elements dimensions
                                      (+ offset (* ncol row-index)))
                      (aref elements row-index offset)))))
          (if (fixnum? row-selection)
              (extract row-selection)
              (if dimensions
                  (map 'vector #'extract row-selection)
                  (sub elements row-selection offset)))))))

(defmethod sub ((mosaic-vector mosaic-vector) &rest selections)
  (let+ (((key-selection &rest subscripts) selections)
         ((&slots-r/o mosaic elements) mosaic-vector))
    (if (eq key-selection t)
        mosaic-vector
        (let+ (((offset . dimensions)
                (mosaic-location mosaic key-selection subscripts)))
          (if dimensions
              (displace-array elements dimensions offset)
              (aref elements offset))))))

(defmethod map-columns (function (matrix mosaic-matrix)
                        &key element-type)
  (let+ (((&slots-r/o mosaic elements) matrix))
    (make-instance 'mosaic-matrix
                   :mosaic mosaic
                   :elements (map-columns function elements
                                          :element-type element-type))))

(defmethod map-rows (function (matrix mosaic-matrix)
                     &key (element-type t) (mosaic nil mosaic?))
  (let+ ((result (map-rows function (elements matrix)
                           :element-type element-type)))
    (when (eq mosaic t)
      (setf mosaic (mosaic matrix)))
    (if (and (typep result 'matrix) mosaic?)
        (make-mosaic-matrix mosaic result)
        result)))

(defmethod quantiles ((matrix mosaic-matrix) qs)
  (map-columns (lambda (c) (quantiles c qs)) matrix))

;; (defclass foo ()
;;   ((a :accessor a :initarg :a)
;;    (b :accessor b :initarg :b)
;;    (c :accessor c :initarg :c)))

;; (defparameter *m* (make-mosaic '((a 1) (b 2 3) c)))

;; (defparameter *f* (make-instance 'foo :a #(1) :b #2A((3 5 7) (13 11 19)) :c 9))

;; (defparameter *g* (make-instance 'foo))

;; (defparameter *v* (make-array (mosaic-size *m*)))

;; (pack-slots *m* *v* *f*)
;; (unpack-slots *m* *v* (find-class 'foo))
;; (unpack-slots *m* *v* ('foo))

;; (type-of (find-class 'foo))

;; (standard-object)
;; (standard-class)

;; (make-instance 'foo)
