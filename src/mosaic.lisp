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
  (let* ((table (make-hash-table))
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

(defun mosaic-location (mosaic key)
  "Return (CONS OFFSET DIMENSIONS).  Read-only, consequences are undefined if
modified."
  (let+ (((&values value present?) (gethash key (mosaic-table mosaic))))
    (assert present? () "Key ~A not found." key)
    value))

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

(defstruct (mosaic-matrix (:constructor make-mosaic-matrix%))
  (mosaic nil :type mosaic :read-only t)
  (elements nil :type matrix :read-only t))

(defun make-mosaic-matrix (mosaic nrow &optional (element-type t))
  "Make a mosaic matrix."
  (make-mosaic-matrix% :mosaic mosaic
                       :elements (make-array (list nrow (mosaic-size mosaic))
                                             :element-type element-type)))

(defmethod pack-slots ((mosaic-matrix mosaic-matrix) (row fixnum) object)
  (let+ (((&structure-r/o mosaic-matrix- mosaic elements) mosaic-matrix))
    (pack-slots mosaic (subarray elements row) object)))

(defmethod unpack-slots ((mosaic-matrix mosaic-matrix) (row fixnum) object)
  (let+ (((&structure-r/o mosaic-matrix- mosaic elements) mosaic-matrix))
    (unpack-slots mosaic (subarray elements row) object)))


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
