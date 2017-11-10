;;; -*- Mode: LISP; Syntax: Common-Lisp; Base: 10; Package: CL-User; -*-
(in-package 'user)
;;;;;
;;;;;   FILE:  mouse:main;defs
;;;;;
;;;;;   DEPENDENT ON:  mouse:lispext
;;;;;
;;;;;

(defparameter *variances*  '(1.0))
(defvar *sigma*)
(defparameter *densities* '(1.0))		; (/ 1.0 (sqr sigma))
(defvar *popn.density*)

(defparameter *trap.shape*  'square)		; ts = {square circle}
(defparameter *forage.areas* '(0.5))		; frequency*patchArea = forageArea
(defvar *forage.area*)
(defparameter *num.squares*     8)		; size of grid = (sqr (* 2 (1- *num.squares*)))
(defparameter *num.days*        4)		; number of trapping days

;;; forage area standard is (/ (sqr sigma) 2.0)

;;;
;;; '(...(#reps frequency-list spacing-list)...) The lists may be atoms.
;;;
(defvar *run.info*   nil)		

(defvar *trap.space*)				; trap.space/sigma
(defvar *frequency*)

;;;
;;;  We model live-trapping by totally random directions at all times.
;;;  We model removal-trapping by moving those mice outside the grid
;;;  only towards the grid (in a rectangular fashion).
;;;
(defparameter *removal.trapping* nil)		; removal trap or live trap
(defvar *d*)					; either distance/trapspace or square.num
(defvar *field.max*)				; (* *trap.space* (- *num.squares* 0.5))
(defvar *field.min*)				; (- *field.max*)

(defvar	*grid*)					; array of tolls indexed on day
(defvar	*squares*)				; vector of tolls indexed on day

(defvar *outfile*)

;;;
;;;   IMPORTANT NOTE:  One refers to the elements of an array (i j)
;;;        as (row col).   We are treating the array *grid* as if it
;;;        were superimposed on the field.  We refer to traps by (x y)
;;;        coordinates with (0 0) at the field center.  It is the col
;;;        index, j, which dictates the x location and the row index,
;;;        i, which dictates the y location.
;;;        Array indices commence at 0.
;;;
;;;   What is stored in the cells of a grid is a vector of tolls:
;;;   #V(day1 day2 ... dayn)  where dayi = number of mice caught on
;;;   day i.  The vector is indexed from 0
;;;

(defmacro trap-tolls (i j)
  `(aref *grid* ,i ,j))

;;;
;;; This sets locations assuming the grid.center is at '(0 . 0)
;;; and that the grid is square 2n x 2n.
;;;
(defmacro trap-location (i j)
  `(cons (* *trap.space* (- ,j *num.squares* -0.5))
	 (* *trap.space* (- *num.squares* ,i 0.5))))

(defmacro trap-x (i j)
  (declare (ignore i))
  `(* *trap.space* (- ,j *num.squares* -0.5)))

(defmacro trap-y (i j)
  (declare (ignore j))
  `(* *trap.space* (- *num.squares* ,i 0.5))) 

(defmacro increment.trap.toll (i j day)
  `(incf (svref (aref *grid* ,i ,j) ,day)))

(defmacro square-tolls (sq)
  `(svref *squares* ,sq))

(defmacro square-toll-day (sq d)
  `(svref (square-tolls ,sq) ,d))

(defmacro increment.square.toll (sq day)
  `(incf (svref (svref *squares* ,sq) ,day)))

;;;
;;;  Want to return the number of the square (ring) the given trap
;;;  (index) belongs.  Squares are numbered from 0 on the inside to
;;;  *num.traps*-1 on the outside.  Indices begin at 0.
;;;
;;;  First we map (i j) to the second quadrant where i,j <= *num.traps*.
;;;  (We fold the grid in half horizontally and again vertically.)  Next
;;;  we upper triangulate so that i <= j <= *num.traps*.  Now i denotes
;;;  the square but it is reversed (increases towards center).  We
;;;  invert the index are done.
;;;
(defmacro square (i j)				; 0 <= i,j < 2*num.squares*
  `(- *num.squares*
      (min (min ,i (- (* 2 *num.squares*) ,i 1))
	   (min ,j (- (* 2 *num.squares*) ,j 1)))
      1 ))

;;;
;;;  Where a point is (x . y) as returned by trap-location.
;;;
(defmacro distance (a b)
  `(sqrt (+ (sqr (- (car ,a) (car ,b)))
	    (sqr (- (cdr ,a) (cdr ,b))))))

(defmacro do-grid (i j &body body)
  `(dotimes (,i (* 2 *num.squares*))
     (dotimes (,j (* 2 *num.squares*))
       ,@body )))


(defmacro trap-is-circle ()
  `(eq *trap.shape* 'circle))

(defmacro trap-is-square ()
  `(eq *trap.shape* 'square))

(defmacro trap-shape-constant ()
  `(if (eq *trap.shape* 'square) 4.0 *pi*))

(defmacro size-of-field ()
  `(+ (* *trap.space* (- (* 2 *num.squares*) 1)) (* 6.0 *sigma*)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                  STATISTICS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *alpha*      0.95)		; confidence coefficient

(defparameter *hri.ratio*  nil)

;;;
;;;  The relative size of the square (not inclusive) 4(tps-1)ss
;;;  where tps is traps-per-side and ss is spacing squared.
;;;
(defmacro relative.size (sq)
  `(* 4 (+ (* 2 ,sq) 1) (sqr *trap.space*)))

;;;
;;;  The sum of the relative areas inclusive = area of the square.
;;;
(defmacro square.area (sq)
  `(sqr (* 2 (+ ,sq 1) *trap.space*)))

(defparameter *histogram*  nil)			; whether to compute or no
(defvar *histogram-width*  20)
(defvar *histogram-height* 40)
(defvar *histogram-max*    4.0)
(defvar *histogram-min*    0.0)

(defvar *distance.list*)
(defmacro clear.distance.list () `(setq *distance.list* nil))
(defmacro add.distance.moved (from to)
  `(push (distance ,from ,to) *distance.list*))


(defvar *stats*)				;'(((X11..Xr1)..(X1n..Xrn))..)
						; Xij ==> random variable X rep i group j
(defmacro clear.stats () `(setq *stats* nil))
(defmacro add.datum.to.stats (datum)
  `(push ,datum *stats*))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Since there is a dangling file, any attempt to terminate
;;;  execution prematurely loses all output accumulated.  This
;;;  flag and function provides and antidote.
;
(defparameter *stop.simulation* nil)
(defmacro stop.simulation-p () '*stop.simulation*)
(defun stop.simulation () (setq *stop.simulation* t) (values))


