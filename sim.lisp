;Mass-Spring Simulator with anisotropic friction, used for modeling of worm loosely based on C.Elegans

; Copyright (C) 2013, Crystalline

;   Permission is hereby granted, free of charge, to any person obtaining a copy of
;this software and associated documentation files (the "Software"), to deal in
;the Software without restriction, including without limitation the rights to
;use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
;of the Software, and to permit persons to whom the Software is furnished to do
;so, subject to the following conditions:

;   The above copyright notice and this permission notice shall be included in all
;copies or substantial portions of the Software.

;   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
;SOFTWARE.

;Controls: arrows for camera rotation, w a s d for movement in XY plane,
;q e for moving in Z axis, F1, F2 for zoom, space for pause

;(declaim (optimize (safety 0) (speed 3) (debug 0)))

;(ql:quickload 'lispbuilder-sdl)
;(ql:quickload 'lispbuilder-sdl-gfx) 

(defpackage :wormsim
   (:use :cl)
   (:export :start-sim :stop-sim :reset-sim :toggle-pause-sim))
(in-package :wormsim)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro def (&rest args) `(defparameter ,@args))

;Worm model tunable parameters

;Surface is XY-plane by default
(def surface-k 5.0)
(def surface-drag 0.28)

(def anisotropic-friction t)
(def surface-drag-tan 0.28)
(def surface-drag-norm 0.01)

(def air-drag 0.01)
(def G 1.0)

;Sim thread is active while this flag is True
(def sim-active t)
(def pause-sim nil)

;Screen parameters
(def height 800)
(def width 800)
(def height/2 (/ height 2))
(def width/2 (/ width 2))
;Camera settings
(def camera-orient (list -0.3 0.5 0.7 0.35))
(def camera-scale 25.0)
(def camera-pos (list -23.0 -29.0 10.0))
(def dphi 0.15)

;Timestep number
(def ts 0)

;Simulation speed
(def timesteps-per-sec 150.0)
(def FPS 30)
(def bgcolor sdl:*black*)

;Controller parameters
(def period 700)
(def start-time 250)
(def freq 2.5)

;Worm parameters

;Middle radius
(def Rmid 1.0)

;End radius
(def Rend 0.59)

(def worm-sections 31)
(def worm-lines 8)
(def worm-section-spacing 1.0)
(def worm-stiffness 30.0)
(def worm-point-mass 0.1)
(def contraction-limit (* 0.70 worm-section-spacing))

;Simulated system is contained in this variable
;elt 0 is list of points
;elt 1 is list of springs
;elt 2 is structured information (muscles etc)
(def sys (list nil nil))

;State variables
(def points nil)

(def indices nil)
(def ground-line nil)
(def left-lines nil)
(def right-lines nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Coloring
(def *orange* (sdl:color :r #xFF :g #xCC :b 66))
(def *pale-blue* (sdl:color :r #xCC :g #xFF :b #xFF))
(def *color-table* nil)

(defun scalar-color (s)
  (when (not *color-table*)
	  (setf *color-table* (make-array 256))
	  (loop for i from 0 to 255 do
	    (setf (aref *color-table* i)
		      (sdl:color :r i :g 0 :b (- 255 i)))))
  (let ((s (if (> s 1.0) 1.0 (if (< s 0.0) 0.0 s))))
	(elt *color-table* (floor (* 255.0 s)))))

(defstruct point
	(mass 1.0 :type float)
	(pos (list 0.0 0.0 0.0))
	(vel (list 0.0 0.0 0.0))
	(force-acc (list 0.0 0.0 0.0))
	(on-ground nil :type boolean)
	(aniso-friction-pair nil)
	(color sdl:*white*))
	
(defstruct spring
	(pa nil)	;Point A
	(pb nil)	;Point B
	(d 1.0 :type float)  ;Equilibrium distance
	(k 1.0 :type float)  ;Stiffness constant
	(color *orange*))

(defun point-get-aniso-normal (p)
  (let ((apair (point-aniso-friction-pair p)))
  (if apair
      (normalize (v- (point-pos apair) (point-pos p)))
	  nil)))
	
(defun apply-force (p F)
   (if p
	(setf (point-force-acc p) (v+ (point-force-acc p) F))))

(defun simple-friction (p)
 (s*v (- surface-drag) (point-vel p)))

;Get aniso friction force
;If there is no aniso-normal (,aniso pair) then falls back to simple friction
(defun aniso-friction (p)
  (let ((anormal (point-get-aniso-normal p))
        (vel (point-vel p)))
    (if anormal
	    (let* ((v-norm (s*v (v.v anormal vel) anormal))
		       (v-tan (v- vel v-norm)))
			(v+ (s*v (- surface-drag-tan) v-tan)
			    (s*v (- surface-drag-norm) v-norm)))
		(simple-friction p))))
	
(defun get-point-force (p)
	(let ((F nil)
		  (Z (.z (point-pos p))))
		;Drag forces
		(if (< Z 0.0)
			(progn
				(setf (point-on-ground p) t)
				(if anisotropic-friction
				  (setf F (aniso-friction p))
				  (setf F (simple-friction p)))
				(setf (elt F 2) (* surface-k (- Z))))
			(progn
			    (setf (point-on-ground p) nil)
				(setf F (s*v (- air-drag) (point-vel p)))))
		;Gravity
		(incf (elt F 2) (* (- G) (point-mass p)))
		F))

;Euler integrator
(defun integrate (p dT)
	(setf (point-vel p)
          (v+ (point-vel p)
              (s*v (/ dT (point-mass p)) (point-force-acc p))))
	(setf (point-pos p) (v+ (point-pos p)
                            (s*v dT (point-vel p))))
	(setf (point-force-acc p) (list 0.0 0.0 0.0)))

(defun timestep (system dT)
  (setf ts (+ ts 1))
  (let ((points (elt system 0))
        (springs (elt system 1)))
	(dolist (s springs)
		(let* ((diff (v- (point-pos (spring-pb s))
			             (point-pos (spring-pa s))))
			   (len (l2-norm diff))
			   (normal (s*v (/ 1.0 len) diff))
			   (spring-force (s*v (* (spring-k s) (- len (spring-d s))) normal)))
			(apply-force (spring-pa s) spring-force)
			(apply-force (spring-pb s) (s*v -1.0 spring-force))))
	(dolist (p points)
		(apply-force p (get-point-force p)))
	(dolist (p points)
		(integrate p dT))))
		
(defun rotate-camera (axis angle)
  (setf camera-orient
        (quaternion* (make-quaternion-rotation axis angle) camera-orient)))

(defun camera-transform (vec)
  (apply-quaternion-rotation camera-orient (s*v camera-scale (v- vec camera-pos))))

(defun draw-line (v1 v2 &key (color sdl:*yellow*))
	(let ((v1 (mapcar #'round (camera-transform v1)))
		  (v2 (mapcar #'round (camera-transform v2))))
	(sdl:draw-line-* (elt v1 0) (- height (elt v1 1))
					 (elt v2 0) (- height (elt v2 1)) :color color :aa nil)))
					 
(defun draw-point (v1 r &key (color sdl:*white*))
	(let ((v1 (mapcar #'round (camera-transform v1))))
		(sdl:draw-filled-circle-* (elt v1 0) (- height (elt v1 1)) r :color color)))
		
(defun draw-grid ()
  (loop for i from -20 to 20 do
	  (draw-line
       (list -20.0 (float i) 0.0) (list 20.0 (float i) 0.0) :color sdl:*cyan*)
	  (draw-line
       (list (float i) -20.0 0.0) (list (float i) 20.0 0.0) :color sdl:*cyan*))
	(draw-line '(0.0 0.0 0.0) '(1.0 0.0 0.0) :color sdl:*red*)
	(draw-line '(0.0 0.0 0.0) '(0.0 1.0 0.0) :color sdl:*green*)
	(draw-line '(0.0 0.0 0.0) '(0.0 0.0 1.0) :color sdl:*blue*))
		
(defun draw-system (system)
	(draw-grid)
	(dolist (s (elt system 1))
		(draw-line (point-pos (spring-pa s))
				   (point-pos (spring-pb s))
				   :color (spring-color s)))
	(dolist (p (elt system 0))
	  (if (point-on-ground p)
		(draw-point (point-pos p) 4 :color sdl:*green*)
		(draw-point (point-pos p) 4 :color (point-color p)))))

(def fading-en t)
		
(defun fade (index)
  (if fading-en
	  (let ((fade-start (floor (* worm-sections 0.75))))
		(if (> index fade-start)
			(+ 1.0
			   (/ (- index fade-start)
				  (- (- worm-sections fade-start))))
			1.0))
	  1.0))
		
(defun control-timestep (dT)
  ;Every 30 timesteps readjust line indices
 ; (when (eq (mod ts 30) 0)
;    (let ((gline (ground-line-index))
	;      (muscle-width (floor (/ worm-lines 4))))
	;  (if (not gline) (setf gline 0))
	;  (setf ground-line (get-line gline))
   ;	  (setf left-line 
  ;t)
  (when (> ts start-time)
	(mapcar (lambda (line) (mapcar (lambda (m index)
							(muscle-contract m
	                          (+ (* (sin (+ (/ (* index freq 2.0 pi) worm-sections)
                                            (/ (* (* 2.0 pi) (mod ts period)) period))) 0.5 (fade index))
								 0.5)))
					   line indices)) left-lines)
    
	(mapcar (lambda (line) (mapcar (lambda (m index)
	                        (muscle-contract m
	                          (+ (* (sin (+ (/ (* index freq 2.0 pi) worm-sections)
                                            (/ (* (* 2.0 pi) (mod ts period)) period)
											pi)) 0.5 (fade index))
								 0.5)))
					   line indices)) right-lines)))

(def time-per-step nil)

(defun sim-thread ()
	(let ((t0 nil)
	      (t1 nil))
	  (loop until (not sim-active) do
			(setf t0 (get-internal-real-time))
			(when (not pause-sim)
				(control-timestep (/ 1.0 timesteps-per-sec))
				(timestep sys (/ 1.0 timesteps-per-sec)))
			(setf t1 (get-internal-real-time))
			(let ((deltaT (- (/ 1.0 timesteps-per-sec) (/ (- t1 t0) 1000.0))))
			    (setf time-per-step (- t1 t0))
				(if (> deltaT 0.0)
					(sleep deltaT))))))

(def delta-ts-save 12)
(def last-save-ts 0)
(def save-screen nil)
(def frame-num 0)

(defun draw-thread ()

(sdl:with-init ()

  (sdl:window	width
				height
				:title-caption "Nematode simulator"
				:FPS (make-instance 'sdl:fps-mixed)
				:FLAGS '(sdl:sdl-doublebuf sdl:sdl-hw-surface))
  
  (setf (sdl:frame-rate) 15)
  
  (sdl:initialise-default-font sdl:*FONT-10X20*)
  
  (draw-system sys)
  
  (sdl:with-events ()
      
	  (:quit-event () (setf sim-active nil) t)
	  
	  ;(:MOUSE-MOTION-EVENT (:x x :y y)
		;(sdl:clear-display bgcolor)
		;(sdl:draw-circle-* x y 20)
		;(sdl:update-display))
	  
	  (:KEY-DOWN-EVENT (:KEY KEY)
      (WHEN (SDL:KEY= KEY :SDL-KEY-ESCAPE)  
        (SDL:PUSH-QUIT-EVENT))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-W)
	    (setf camera-pos (v+ camera-pos (list 1.0 0.0 0.0))))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-S)
	    (setf camera-pos (v+ camera-pos (list -1.0 0.0 0.0))))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-A)
	    (setf camera-pos (v+ camera-pos (list 0.0 1.0 0.0))))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-D)
	    (setf camera-pos (v+ camera-pos (list 0.0 -1.0 0.0))))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-Q)
	    (setf camera-pos (v+ camera-pos (list 0.0 0.0 1.0))))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-E)
	    (setf camera-pos (v+ camera-pos (list 0.0 0.0 -1.0))))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-UP)  
        (rotate-camera '(1.0 0.0 0.0) dphi))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-DOWN)  
        (rotate-camera '(1.0 0.0 0.0) (- dphi)))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-LEFT)  
        (rotate-camera '(0.0 1.0 0.0) dphi))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-RIGHT)  
        (rotate-camera '(0.0 1.0 0.0) (- dphi)))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-F1)  
        (setf camera-scale (* 1.2 camera-scale)))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-F2)  
        (setf camera-scale (* 0.8 camera-scale)))
	  (WHEN (SDL:KEY= KEY :SDL-KEY-SPACE)  
        (setf pause-sim (not pause-sim))))
      
	  (:video-expose-event () (sdl:update-display))
	  
	  (:idle
		
		(sdl:clear-display bgcolor)
			
		;;Draw objects
		(draw-system sys)
		
		(sdl:draw-string-solid-* (format nil "Model Time:~s" ts) 2 2 :color sdl:*white*)
		(sdl:draw-string-solid-* (format nil "CPU Time/Step:~s ms" time-per-step) 2 22 :color sdl:*white*)
		(sdl:draw-string-solid-* (format nil "Cam pos:~s" camera-pos) 2 42 :color sdl:*white*)
		(sdl:draw-string-solid-* (format nil "Cam scale:~s" camera-scale) 2 62 :color sdl:*white*)
		
		(when save-screen
		  (when (> (- ts last-save-ts) delta-ts-save)
		    (setf last-save-ts ts)
			(incf frame-num 1)
			(sdl:save-image SDL:*DEFAULT-DISPLAY* (format nil "frame~4,'0D.bmp" frame-num))))
		    
		
		(sleep 0.01)
		
		;; Flip back/front buffers
		(sdl:update-display)))))


(defun make-ring-z (R N mass k)
	(let* ((points nil)
		   (springs nil)
		   (phi (/ (* 2.0 pi) N))
		   (rvec (list R 0.0 0.0))
		   (axis (list 0.0 0.0 1.0))
		   (dist (* 2.0 R (sin (/ phi 2.0)))))
	  (loop for i from 0 to (- N 1) do
		(push (make-point :mass mass :pos (rotate-axis-angle axis (* i phi) rvec)
					   :vel (list 0.0 0.0 0.0) :force-acc (list 0.0 0.0 0.0))
			   points))
	  (loop for i from 0 to (- N 2) do
		(push (make-spring :pa (elt points i) :pb (elt points (+ i 1))
						:d dist :k k) springs))
							
	  (push (make-spring :pa (elt points (- N 1)) :pb (elt points 0)
					     :d dist :k k) springs)
	  ;Center point
	  (push (make-point :mass mass :pos (list 0.0 0.0 0.0))
			points)
	  ;Radial springs
	  (loop for i from 0 to (- N 1) do
		(push (make-spring :pa (elt points 0) :pb (elt points (+ i 1))
						   :d R :k k) springs))
	  (list points springs)))

(defun rotate-sys (axis phi sys)
   (dolist (p (elt sys 0))
       (setf (point-pos p) (rotate-axis-angle axis phi (point-pos p))))
   sys)

(defun translate-sys (delta sys)
   (dolist (p (elt sys 0))
       (setf (point-pos p) (v+ delta (point-pos p))))
   sys)

(defun get-ring-radius (ring)
  (l2-norm (v- (point-pos (elt (elt ring 0) 0)) (point-pos (elt (elt ring 0) 1)))))
   
;Returns list of springs and muscles
;Sets anisotropic friction pairs
(defun link-rings (A B dist k)
   (if (eq (length (elt A 0)) (length (elt B 0)))
     (let* ((springs nil)	;Bulk list of spring for simulation
	        (muscles nil)	;List of muscles. Muscle is a list of two points and three springs that can be actuated
		    (Ra (get-ring-radius A))
		    (Rb (get-ring-radius B))
	        (diagonal-A (sqrt (+ (expt dist 2) (expt Ra 2))))
			(diagonal-B (sqrt (+ (expt dist 2) (expt Rb 2)))))
	   (loop for i from 0 to (- (length (elt A 0)) 1) do
	     (let ((spr (make-spring :pa (elt (elt A 0) i) :pb (elt (elt B 0) i)
				        	:d dist :k k)))
		   ;Push spring to spring list
		   (push spr springs)
		   ;Add aniso-friction pair
		   (setf (point-aniso-friction-pair (elt (elt A 0) i)) (elt (elt B 0) i))
		   ;Create muscle if not axial fiber
		   (when (>= i 1)
		     (push (list (elt (elt A 0) i) (elt (elt B 0) i) (list spr dist)) muscles))))
	   (loop for i from 1 to (- (length (elt A 0)) 1) do
	     (let ((spA (make-spring :pa (elt (elt A 0) 0) :pb (elt (elt B 0) i)
				        	:d diagonal-A :k k))
			   (spB (make-spring :pa (elt (elt A 0) i) :pb (elt (elt B 0) 0)
				        	:d diagonal-B :k k)))
	       (push spA springs)
	       (push spB springs)
		   (setf (elt muscles (- i 1)) (append (elt muscles (- i 1)) (list (list spA diagonal-A Ra) (list spB diagonal-B Rb))))))
	   (list springs muscles))
	   (print "Error in link-rings~%")))

;Translate point in place
(defun point-translatef (p v)
  (setf (point-pos p) (v+ (point-pos p) v))
  p)
	   
;L is the number of rings, N is a number of sections per ring R is radius of ring, dist is distance between rings	   
(defun make-worm (L N dist k mass proc-R)
	(let ((rings nil)				;Rings of nematoda
	      (axial-springs nil)		;Axial string, goes through whole worm
		  (lines (make-list N)))	;Lines are lists of muscles belonging to one of N outer strings
	  (loop for i from 0 to (- L 1) do
		(let ((R (funcall proc-R (/ i L))))
		(push (translate-sys (list (* i dist) 0.0 1.0)
		                     (rotate-sys '(0.0 1.0 0.0) (/ pi 2.0) (make-ring-z R N mass k)))
			  rings)))
	  (loop for i from 0 to (- L 2) do
	    (let* ((link-res (link-rings (elt rings i) (elt rings (+ i 1)) dist k))
		       (springs (car link-res))
			   (muscles (cadr link-res)))
			(setf axial-springs (append springs axial-springs))
			
			;Push muscles to lines
			(loop for j from 0 to (- N 1) do
			  (push (elt muscles j) (elt lines j)))))
	  (list
		(mapcar (lambda (p) (point-translatef p (list (* (/ L -2.0) dist) 0.0 0.0))) ;Translate worm to center
		        (reduce #'append (mapcar #'car rings)))
		(append (reduce #'append (mapcar #'cadr rings)) axial-springs)
		lines)))

(defun set-muscle-color (M color)
  (dolist (spr (cddr M))
    (setf (spring-color (car spr)) color)))

(defun set-line-color (index color)
  (dolist (muscle (elt (elt sys 2) index))
    (set-muscle-color muscle color)))
	
(defun get-muscle (line-index muscle-index)
  (elt (elt (elt sys 2) line-index) muscle-index))

;Maps [0..1] to worm radius
(defun R-worm (x)
  (let* ((Kl (/ (- Rmid Rend) 0.5))
         (Bl Rend)
		 (Br (+ Rmid (/ Kl 2.0))))
	  (if (< x 0.5)
		  (+ (* Kl x) Bl)
		  (+ (* (- Kl) x) Br))))

;Maps muscle activation to spring equi distance
(defun act-mapping-L2 (s)
  (- worm-section-spacing (* s (- worm-section-spacing contraction-limit))))
  
(defun muscle-contract (muscle s)
  (let* ((s (if (> s 1.0) 1.0 (if (< s 0.0) 0.0 s)))
		 (color (scalar-color s)))
	(setf (spring-color (car (elt muscle 2))) color)
	(setf (point-color (elt muscle 0)) color)
	(setf (point-color (elt muscle 1)) color)
	(setf (spring-d (car (elt muscle 2))) (act-mapping-L2 s))))

(defun muscle-contract-index (line index s)
  (muscle-contract (get-muscle line index) s))

(defun muscle-on-groundp (muscle)
  (and (point-on-ground (elt muscle 0))
       (point-on-ground (elt muscle 1))))
  
(defun get-line (index) (if index (elt (elt sys 2) index) nil))
  
(defun line-contract (line s)
  (mapcar (lambda (m) (muscle-contract m s)) line)
  t)
  
(defun line-contract-index (index s) (line-contract (get-line index) s))

;If half of line's muscle is on the ground then whole line is on the ground
(defun line-on-groundp (line)
  (> (count-if #'muscle-on-groundp line)
    (/ (length line) 2)))

 (defun ground-line-index ()
   (let ((res 0))
     (loop for i from 0 to (length (elt sys 2)) do
         (if (line-on-groundp (get-line i))
			(setf res i)))
	res))

(def sys nil)
(def indices nil)
(def ground-line nil)
(def left-lines nil)
(def right-lines nil)

(defun create-initial-conditions ()
  (setf sys (make-worm worm-sections worm-lines
                       worm-section-spacing
                       worm-stiffness worm-point-mass
                       #'R-worm))
;Muscle indices in lines
  (setf indices (build-list worm-sections (lambda (x) x)))
  (setf ground-line (get-line 0))
  (setf left-lines (list (get-line 1) (get-line 2)))
  (setf right-lines (list (get-line 5) (get-line 6)))
  t)

;Running simulation and graphics threads
(def sim-thread nil)
(def gr-thread nil)

(defun start-sim ()
  (when (not sys)
    (create-initial-conditions))
  (setf sim-active t)
  (setf sim-thread (sb-thread:make-thread  #'sim-thread))
  (setf gr-thread (sb-thread:make-thread  #'draw-thread)))

(defun stop-sim ()
  (setf sim-active nil)
  (when (sb-thread:thread-alive-p gr-thread)
  	(sb-thread:terminate-thread gr-thread))
  (when (sb-thread:thread-alive-p sim-thread)
	(sb-thread:terminate-thread sim-thread)))

(defun reset-sim ()
  (create-initial-conditions))

(defun toggle-pause-sim ()
  (setf pause-sim (not pause-sim)))
