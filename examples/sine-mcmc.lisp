;;; sine example
(in-package :cl-bayesian)
(asdf:load-system :cl-cairo2-x11)
(asdf:load-system :cl-2d)

(define-mcmc sine ()
    ((x :parameter (atom :updater metropolis) :reader x :initarg :x)))

;; (define-mcmc sinexy ()
;;     ((x :parameter (atom :updater metropolis) :reader x :initarg :x)
;;      (y :parameter (vector :updater gibbs))))


(define-metropolis-updater (sine x) ()
    (flet ((log-likelihood (x)
             (if (< 2 x 100)
                 (log (/ (sin (/ (* 2 pi) x)) (expt x 2)))
                 nil)))
      (let* ((xnext (+ x (rv:draw x-propdist)))
             (lp-ratio (calc-lp-ratio x xnext log-likelihood)))
        (metropolis-step x xnext lp-ratio x-counter))))
 
(defmethod copy-parameters ((mcmc sine))
  (x mcmc))

(defparameter *a* (make-instance 'sine :x 2.5d0 :x-propdist (make-instance 'rv:normal)))

(defparameter *lambda* (run-mcmc *a* 100000))
(defparameter *r* (xcollect (length *lambda*) (rv:generator* 'rv:beta :alpha 3d0 :beta 1d0)))
  
(defparameter *phi1* (xmap '(array :element-type double-float)
                           (lambda (r l)
                             (* 2d0 r (cos (/ (* 2 pi) l))))
                           *r* *lambda*))
(defparameter *phi2* (map 'vector (lambda (r)
                                    (- (expt r 2)))
                          *r*))


(defparameter *frame* (cl-2d:as-frame (cl-cairo2:create-xlib-image-context 800 600)  
                                      :background-color cl-colors:+white+))

(cl-2d:plot-function *frame* (lambda (x) (log (/ (sin (/ (* 2 pi) x)) (expt x 2))))
                     (cl-2d:interval-of 2 10))

(cl-2d:plot-histogram *frame* (cl-numlib:histogram *lambda* :breaks-function 
                                                   (cl-numlib:histogram-evenly-distributed-breaks 100))
                      :x-interval (cl-2d:interval-of 2 8))
(cl-2d:plot-sequence *frame* *lambda*)
(cl-2d:plot-histogram *frame* (cl-numlib:histogram *r*))
(cl-2d:plot-histogram *frame* (cl-numlib:histogram *phi1*))
(cl-2d:plot-histogram *frame* (cl-numlib:histogram *phi2*))

(cl-2d:plot-symbols *frame* *phi1* *phi2*
                    :symbol-drawing-function #'cl-2d:symbol-filled-circle
                    :size-function (constantly 20)
                    :color-function (constantly (cl-colors:add-alpha cl-colors:+blue+ 0.01)))
