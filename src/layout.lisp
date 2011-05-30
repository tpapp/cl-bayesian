;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:cl-bayesian)

;;; Layouts are a means for flattening collections of objects into vectors.,
;;; then reconstructing them again.

(defun displace-subvector (vector start end)
  (displace-array vector (- end start) start))

(defgeneric extract (vector layout keys &key copy?)
  (:documentation "Return object from vector, extracted using KEYS and
  LAYOUT."))

(defgeneric layout-length (layout)
  (:documentation "Return the number of elements in layout."))

(defgeneric parse-layout (object)
  (:documentation "Return a layout corresponding to OBJECT.  This is used to
  define a DSL for constructing layouts, unrecognized objects passed
  through by default.")
  (:method (object)
    object))

(defgeneric flatten-into (vector layout object)
  (:documentation "Copy elements of object into vector at offset."))

;;; named layout -- corresponds to a plist

(defstruct+ (named-layout (:constructor named-layout%))
  "Layout where each sub-layout is named by a symbol."
  (keys nil :type vector :read-only t)
  (offsets nil :type simple-fixnum-vector :read-only t)
  (layouts nil :type vector :read-only t))

(defun named-layout (&rest plist)
  "Create a named layout from a plist of layout items."
  (iter
    (with offset := 0)
    (for (key maybe-layout &rest rest) :on plist :by #'cddr)
                                        ; !! check for dangling elements?
    (let ((layout (parse-layout maybe-layout)))
      (collect key :into keys :result-type vector)
      (collect (incf offset (layout-length layout)) :into offsets
               :result-type simple-fixnum-vector)
      (collect layout :into layouts :result-type vector))
    (finally
     (return (named-layout% :keys keys :offsets offsets :layouts layouts)))))

(defmethod layout-length ((named-layout named-layout))
  (vector-last (named-layout-offsets named-layout)))

(defmethod extract (vector (named-layout named-layout) keys &key copy?)
  (declare (optimize debug))
  (let+ (((&named-layout-r/o keys% offsets layouts) named-layout)
         ((first &rest rest) keys)
         (index (aprog1 (position first keys%)
                  (assert it () "key not found")))
         (start (if (zerop index)
                    0
                    (aref offsets (1- index))))
         (end (aref offsets index)))
    (extract (displace-subvector vector start end)
             (aref layouts index) rest :copy? copy?)))

(defmethod extract (vector (named-layout named-layout) (keys null) &key copy?)
  (let+ (((&named-layout-r/o keys offsets layouts) named-layout))
    (iter
      (for key :in-vector keys)
      (for end :in-vector offsets)
      (for start :previous end :initially 0)
      (for layout :in-vector layouts)
      (collect key)
      (collect (extract (displace-subvector vector start end) layout nil
                        :copy? copy?)))))

(defmethod flatten-into (vector (named-layout named-layout) (list list))
  (let+ (((&named-layout-r/o keys offsets layouts) named-layout))
    (iter
      (for (key% object &rest rest) :on list :by #'cddr)
      (for key :in-vector keys)
      (for end :in-vector offsets)
      (for start :previous end :initially 0)
      (for layout :in-vector layouts)
      (assert (eq key key%))
      (flatten-into (displace-subvector vector start end) layout object))))

;;; scalar layout

(defstruct (scalar-layout 
             (:constructor scalar-layout (&optional type)))
  (type t))

(defmethod layout-length ((scalar-layout scalar-layout))
  1)

(defmethod parse-layout ((object (eql 0)))
  (scalar-layout))

(defmethod extract (vector (scalar-layout scalar-layout) keys &key copy?)
  (declare (ignore copy?))
  (assert (null keys) ())
  (aref vector 0))

(defmethod flatten-into (vector (scalar-layout scalar-layout) (number number))
  (assert (= (length vector) 1))
  (setf (aref vector 0) number))

;;; array layout

(defstruct+ (array-layout (:constructor array-layout%))
  "Layout for elements in a row-major order."
  (dimensions nil :type list)
  (row-major-coefficients nil :type vector))

(defun array-layout (dimensions)
  (array-layout% :dimensions dimensions
                 :row-major-coefficients (row-major-coefficients dimensions)))

(defmethod layout-length ((array-layout array-layout))
  (product (array-layout-dimensions array-layout)))

(defmethod parse-layout ((dimensions list))
  (array-layout dimensions))

(defmethod extract (vector (array-layout array-layout) keys
                    &key copy?)
  (let+ (((&array-layout-r/o dimensions row-major-coefficients) array-layout))
    (iter
      (for leftover-dimensions :on dimensions)
      (for key :in keys)
      (for row-major-coefficient :in-vector row-major-coefficients)
      (assert (within? 0 key (car leftover-dimensions)))
      (summing (* key row-major-coefficient) :into sum)
      (finally 
       (return (cl-num-utils::maybe-copy-array
                (displace-array vector leftover-dimensions sum)
                copy?))))))

(defmethod flatten-into (vector (array-layout array-layout) (array array))
  (assert (= (length vector) (array-total-size array)))
  (replace vector (flatten-array array)))

;;; testing

;; (defparameter *a* (make-array* '(3 5) t (let ((i 0))
;;                                           (lambda () (prog1 i (incf i))))))
;; (defparameter *b* 200)

;; (defparameter *l* (named-layout :a '(3 5) :b 0))

;; (defparameter *v* (make-array (layout-length *l*)))

;; (flatten-into *v* *l* (list :a *a* :b *b*))

;; (extract *v* *l* '(:a))


