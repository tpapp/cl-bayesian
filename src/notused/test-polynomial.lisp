(asdf:oos 'asdf:load-op :lift)
(asdf:oos 'asdf:load-op :cl-timeseries)

(in-package :cl-timeseries)

(use-package :lift)

;;; HELPER FUNCTIONS

;;; TEST SUITE

(deftestsuite cl-timeseries () ()
  :equality-test #'equalp)

;;; polynomials

(addtest (cl-timeseries)
  polynomials
  (let ((a (make-poly 2d0))
        (b (make-poly 3d0))
        (c (make-poly 5d0 7d0))
        (d (make-poly 11d0 13d0 17d0)))
    (ensure-same (poly* a b) #(5.0d0 6.0d0))
    (ensure-same (poly* b c) #(8.0d0 22.0d0 21.0d0))
    (ensure-same (poly* c d) #(16.0d0 75.0d0 159.0d0 176.0d0 119.0d0))))

(addtest (cl-timeseries)
  filter
  (let ((x (make-poly 1d0 2d0 3d0 3d0 5d0 6d0))); abusing make-poly :-)
    (ensure-same (filter x (make-poly))
                 #(1.0d0 2.0d0 3.0d0 3.0d0 5.0d0 6.0d0))
    (ensure-same (filter x (make-poly 0d0))
                 #(2.0d0 3.0d0 3.0d0 5.0d0 6.0d0))
    (ensure-same (filter x (make-poly -1d0))
                 #(1.0d0 1.0d0 0.0d0 2.0d0 1.0d0))
    (ensure-same (filter x (make-poly -0.5d0 -0.5d0))
                 #(1.5d0 0.5d0 2.0d0 2.0d0))))
