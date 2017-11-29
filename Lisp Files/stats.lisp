;;; Copyright (c) 2017 Michael Barbehenn

;;; -*- Mode: LISP; Syntax: Common-Lisp; Base: 10; Package: CL-USER; -*-
(in-package 'user)

;;;;;
;;;;;   FILE:  mouse:main;stats
;;;;;
;;;;;   DEPENDENT ON:  mouse:main;defs
;;;;;                  mouse:studentt;main
;;;;;                  mouse:stdnorm;main
;;;;;
;;;;;   NOTES:
;;;;;

;;;;;
;;;;;  Whoever calls popn-stats, pushes result onto a list.
;;;;;  When it is done collecting repetitions, it processes the list
;;;;;  from '(((N11 A11 D11 P11)...) ... (...(Nrn' Arn' Drn' Prn')))
;;;;;  to '(((N11..Nr1)...) ... (...(P1n'..Prn'))).
;;;;;  Remove the Ahat estimates and feed into statistics below.
;;;;;
;;;;;  In other words, for each repetition and for each square, store
;;;;;  the estimates.  Collect common estimator values and group them by
;;;;;  square.  Some estimates were unattainable and are nil (an atom).
;;;;;


;;;
;;;  dHat is the estimated population density
;;;  pHat is the estimated probability of capture
;;;  nHat is the estimated population size
;;;  aHat is the estimated field area
;;;   hri is the home range index of Calhoun
;;;
;;;  According to Barbehenn (1974) four days of trapping may be
;;;  condensed into two periods of two days each.  These are pd1 and pd2.
;;;
;;;

(defmacro ZERO-NEG-PHAT (phat)
  `(and ,phat (if (minusp ,phat) 0 ,phat)))

(defmacro DROP-NEG-PHAT (phat)
  `(and ,phat (unless (minusp ,phat) ,phat)))


(defun popn.estimators ()
  "
  Assumes *squares* is set up.

  OUTPUT: '((N1 A1 D1 P1k P1d P1z)...(Nn An Dn Pnk Pnd Pnz) (Nn' An' Dn' Pn' Pn'k Pn'd Pn'z))
          where n = *num.squares*, n' = non-cumulative last square
          and k=keep neg phats, d=drop neg phats, and z=zero neg phats
  "
  (if (= *num.days* 4)
      (let (pHat dHat nHat aHat (pd1 0.0) (pd2 0.0) stats)
	;; cumulative incremental squares
	(dotimes (sq *num.squares*)
	  (setq pd1 (+ pd1 (square-toll-day sq 0) (square-toll-day sq 1))
		pd2 (+ pd2 (square-toll-day sq 2) (square-toll-day sq 3))
		aHat (square.area sq))
	  (if (= pd1 pd2)
	      (setq nHat nil dHat nil)
	      (setq nhat (/ (sqr pd1) (- pd1 pd2)) dHat (/ nHat ahat)))
	  (if (= 0 pd1)
	      (setq pHat nil)
	      (setq pHat (- 1.0 (sqrt (/ pd2 pd1)))))
	  (push (list nhat aHat dHat pHat (drop-neg-phat pHat) (zero-neg-phat pHat))
		stats))
	;; outer square in isolation
	(setq pd1 (+ (square-toll-day (1- *num.squares*) 0)
		     (square-toll-day (1- *num.squares*) 1)))
	(setq pd2 (+ (square-toll-day (1- *num.squares*) 2)
		     (square-toll-day (1- *num.squares*) 3)))
	(setq aHat (relative.size (1- *num.squares*)))
	(if (= pd1 pd2)
	    (setq nHat nil dHat nil)
	    (setq nhat (/ (sqr pd1) (- pd1 pd2)) dHat (/ nHat ahat)))
	(if (= 0 pd1)
	    (setq pHat nil)
	    (setq pHat (- 1.0 (sqrt (/ pd2 pd1)))))
	(push (list nHat aHat dHat pHat (drop-neg-phat pHat) (zero-neg-phat pHat))
	      stats)
	(nreverse stats))

      (format t "~%Sorry, no general statistics are currently ~
                   implemented or other than four days of trapping~%")))


;;;;;
;;;;;  Given the data from r repetitions on n groups, get the statistics
;;;;;  per group.  NOTE:  some of the repetitions have nil values.  These
;;;;;  should be removed BEFORE calling statistics:  delete-if #'null
;;;;;

;;;
;;; Takes Globally known, collective output from popn-estimates, and
;;; turns it into the input to statistics.  It returns the output from
;;; statistics on nHat aHat pHat and dHat.
;;;
(defun stats.of.stats ()
  (let ((stats (collect-estimators *stats*))
	 result)
    (dolist (var stats)
      (push (statistics *alpha* var) result))
    (nreverse result)))


;;;
;;;
;;; INPUT:   '((X11..Xr1) ... (X1n..Xrn))  
;;;
;;; OUTPUT:  '((meanX1 std.devX1 kX1 minX1 maxX1) ...)
;;;
(defun statistics (alpha n-by-r)		;n groups, up to r each
  (let (result)
    (dolist (group n-by-r)
      (multiple-value-bind (mean std-dev k min max)
	  (ST.graph.values alpha group)
	(push (list mean std-dev k min max) result)))
    (nreverse result)))


;;;
;;;  Given lists of sets of data points, each list is interpreted as
;;;  random samples (values distributed iid) for related random variables.
;;;  For each sample analyzed, a list is returned containing the population
;;;  mean and standard deviation, the 95% symmetric confidence interval delta,
;;;  and the minimum and maximum values of the sample.
;;;  
;;;  Processes the list
;;;  from '(((N11 A11 D11 P11 ..)...) ... (...(Nrn Arn Drn Prn ..)))
;;;   to  '(((N11..Nr1)...) ... (...(P1n..Prn)...))
;;;
;;;  FROM r sets of n groups of vars TO Var times n groups of up to r
;;;  elements (delete-if #'null)
;;;

;;;
;;;    Where Xij is Variable[rep, group]
;;;
;;;  start  '(((N11 A11 D11 P11 ..)...) ... (...(Nrn Arn Drn Prn ..)))
;;;   then  '(((N11..N1n)..(P11..P1n)..) ... )
;;;   then  '(((N11..N1n)..(Nr1..Nrn)) ... )
;;;   then  '(((N11..Nr1)..(N1n..Nrn)) ...)
;;;
(defun collect-estimators (estimators-list)
  (let (result work1 work2)

    ;; separate variables by repetition
    (dolist (rep estimators-list)
      (push (mapcar* #'list rep) work1))

    ;; collect variables across repetitions
    (setq work2 (mapcar* #'list work1))

    ;; collect variables by group
    (dolist (var work2)
      (push (mapcar* #'list var) result))

    ;; delete nil values
    (dotimes (i (length result))
      (dotimes (j (length (nth i result)))
	(setf (nth j (nth i result))
	      (delete-if #'null (nth j (nth i result))))))

    (nreverse result)))

