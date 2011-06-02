(in-package :cl-timeseries)

(define-mcmc armodel ()
  (;; data
   (x :accessor x :type vector :documentation "time series")
   ;; parameters
   (real-roots :accessor real-roots :type vector
               :documentation "real roots"
               :parameter (vector :updater gibbs))
   (complex-roots :accessor complex-roots :type vector
                  :parameter (vector :updater metropolis)
                  :documentation
                  "complex roots, vector of consecutive pairs, stored
                  as magnitude/period conses.")
   (real-probabilities :accessor real-probabilities :type vector
                       :parameter (atom :updater gibbs)
                       :documentation "pr-1, pr0, and pr1")
   (complex-probabilities :accessor complex-probabilities :type vector
                          :parameter (atom :updater gibbs)
                          :documentation "pc0 and pc1")
   (variance :accessor variance :type double-float
             :parameter (atom :updater gibbs)
             :documentation "variance of innovation")
   ;; dogmatic parameters
   (lambda-upper :accessor lambda-upper :type double-float
                 :documentation "upper bound for periods")))

(defun real-root->poly (real-root)
  "Return polynomial for real root."
  (make-poly (- real-root)))

(defun complex-root->poly (complex-root)
  "Return polynomial for complex root."
  (let+ (((r . lambda) complex-root))
    (make-poly (* -2d0 r (cos (/ (* 2 pi) lambda)))
               (expt r 2))))
    
(defun roots->polynomial (roots complex-or-real &optional remove-index)
  "Return the polynomial for roots, removing root at remove-index if
non-nil."
  (let* ((root->poly (ecase complex-or-real
                       (:real #'real-root->poly)
                       (:complex #'complex-root->poly)))
         (polynomials (iter
                       (for i :from 0)
                       (for r :in-vector roots)
                       (when (and remove-index (= i remove-index))
                         (collecting (funcall root->poly r))))))
    (reduce #'poly* polynomials)))

(defun draw-real-root (mean var probabilities)
  "Draw a real root, with given mean and variance for the likelihood,
and probabilities for mass points (vector of 3 elements)."
  (flet ((calculate-p (x p-index)
           (* (aref probabilities p-index) (exp (/ (expt (- x mean) 2) var -2d0)))))
    (let* ((sd (sqrt var))
           (p-1 (calculate-p -1d0 0))
           (p0 (calculate-p 0d0 1))
           (p1 (calculate-p 1 2))
           (prest (* (- 1 (xsum pi)) 0.5d0 (sqrt (* 2d0 pi var))
                     (- (rv:cdf 'rv:normal 1d0) (rv:cdf 'rv:normal -1d0)))))
      (ecase (rv:draw* 'discrete :probabilities (vector p-1 p0 p1 prest))
        (0 -1d0)
        (1 0d0)
        (2 1d0)
        (3 (rv:draw* 'truncated-normal :mu mean :sigma sd
                  :left -1d0 :right 1d0))))))

(define-updater (armodel real-roots :vector-index i) 
    (real-roots complex-roots real-probabilities variance x)
  "Updater for real roots."
  (let+ ((poly (poly* (roots->polynomial real-roots :real i)
                      (roots->polynomial complex-roots :complex)))
         (w (filter x poly))
         (y (take 'numeric-vector (slice w '(1 -1))))
         (x (take 'numeric-vector (slice w '(0 -2))))
         ((&values beta qr nil) (least-squares y x))
         (var (* variance (xref (least-squares-raw-variance qr) 0 0))))
    (draw-real-root (xref beta 0) var real-probabilities)))

(define-updater (armodel complex-roots :vector-index i)
    (real-roots complex-roots complex-probabilities variance x)
  "Updater for complex roots."
  (let+ ((poly (poly* (roots->polynomial real-roots :real i)
                      (roots->polynomial complex-roots :complex)))
         (w (filter x poly))
         (y 
