;(defvar pi 3.14159265)

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

(in-package :wormsim)

(defun @ (a i) (elt a i))

(defun I (x) x)

(defun build-list-helper (size proc)
  (if (eq size -1)
      nil
	  (cons (funcall proc size) (build-list-helper (- size 1) proc))))

(defun build-list (size proc)
  (let ((result (make-list size)))
    (loop for i from 0 to (- size 1) do
	  (setf (elt result i) (funcall proc i)))
	  result))

(defun .x (v) (elt v 0))
(defun .y (v) (elt v 1))
(defun .z (v) (elt v 2))

;Vectors
(defun l2-norm (lst)
  (sqrt (reduce #'+ (mapcar (lambda (x) (* x x)) lst))))

(defun m*v (mat lst)
  (mapcar (lambda (row) (reduce #'+ (mapcar * row lst)))
       mat))
	   
(defun s*v (s v)
  (mapcar (lambda (x) (* s x)) v))

(defun normalize (v)
  (let ((norm (l2-norm v)))
    (if (eq norm 0.0)
	    '(0.0 0.0 0.0)
		 (s*v (/ 1.0 norm) v))))

(defun v+v (A B)
  (mapcar #'+ A B))

(defun v+ (&rest args)
  (reduce #'v+v args))

(defun v- (A B)
      (mapcar #'- A B))

(defun v.v (A B)
  (reduce #'+ (mapcar #'* A B)))

(defun cross* (A B)
  (let ((a1 (elt A 0)) (a2 (elt A 1)) (a3 (elt A 2))
        (b1 (elt B 0)) (b2 (elt B 1)) (b3 (elt B 2)))
    (list (- (* a2 b3) (* a3 b2))
          (- (* a3 b1) (* a1 b3))
          (- (* a1 b2) (* a2 b1)))))
		  
(defun project (target vect)
  (s*v (v.v target vect) (normalize target)))

(defun normal-part (normal vect) (project normal vect))

(defun tangential-part (normal vect) (v- vect (project normal vect)))

;Matrices
(defun take-row (A i)
  (@ A i))

(defun take-col (A j)
  (mapcar (lambda (row) (elt row j)) A))

(defun row*col (A i B j)
  (v.v (take-row A i) (take-col B j)))

(defun build-mat (M N proc)
  (build-list M (lambda (i) (build-list N (lambda (j) (funcall proc i j))))))

(defun n-rows (mat)
  (length mat))

(defun n-cols (mat)
  (length (car mat)))

(defun m*m (A B)
  (build-mat (n-rows A)
             (n-cols B)
             (lambda (i j) (row*col A i B j))))

;(defun (print-v3 lst)
;  (printf "[X:~s Y:~s Z:~s]\n" (@ lst 0) (@ lst 1) (@ lst 2)))

(defun ang->rad (ang)
  (* 2 pi (/ ang 360)))

;Rotation matrices 
(defun rot-x-m3 (phi)
  (list
   (list 1.0 0.0 0.0)
   (list 0.0 (cos phi) (- (sin phi)))
   (list 0.0 (sin phi)    (cos phi))))

(defun rot-y-m3 (phi)
  (list
   (list (cos phi) 0.0 (sin phi))
   (list 0.0 1.0 0.0)
   (list (- (sin phi)) 0.0 (cos phi))))

(defun rot-z-m3 (phi)
  (list
   (list (cos phi) (- (sin phi)) 0.0)
   (list (sin phi)    (cos phi) 0.0)
   (list 0.0 0.0 1.0)))

(defun rot-matrix (rot-v)
  (m*m (rot-z-m3 (@ rot-v 2))
       (m*m (rot-y-m3 (@ rot-v 1)) (rot-x-m3 (@ rot-v 0)))))

(defun rotate-angles (angles v)
  (m*v (rot-matrix angles) v))

;Rodriguez formula
(defun rotate-axis-angle (unit-axis phi v)
  (let ((k unit-axis))
    (v+ (s*v (cos phi) v)
        (s*v (sin phi) (cross* k v))
        (s*v (* (v.v k v) (- 1.0 (cos phi))) k))))

;Quaternions
(defun quaternion+ (A B)
  (mapcar #'+ A B))
  
(defun quaternion* (A B)
  (let ((s1 (car A)) (s2 (car B))
        (v1 (cdr A)) (v2 (cdr B)))
    (cons (- (* s1 s2) (v.v v1 v2))
          (v+ (s*v s1 v2)
              (s*v s2 v1)
              (cross* v1 v2)))))

(defun quaternion-conjugate (A)
  (cons (car A) (mapcar #'- (cdr A))))

(defun quaternion-inverse (A)
  (s*v (/ 1.0 (l2-norm A)) (quaternion-conjugate A)))

(defun make-quaternion-rotation (unit-axis phi)
  (cons (cos (/ phi 2)) (s*v (sin (/ phi 2)) unit-axis)))

(defun apply-quaternion-rotation (q v)
  (cdr (quaternion* (quaternion* q (cons 0.0 v))
                    (quaternion-conjugate q))))

(defun S (x) (* x x))
					
(defun quaternion->rot-matrix (q)
  (let* ((s (car q)) (v (cdr q))
                     (vx (.x v))
                     (vy (.y v))
                     (vz (.z v)))
    (list (list (- 1.0 (* 2.0 (+ (S vy) (S vz)))) (* 2.0 (- (* vx vy) (* s vz))) (* 2.0 (+ (* vx vz) (* s vy))))
          (list (* 2.0 (+ (* vx vy) (* s vz))) (- 1.0 (* 2.0 (+ (S vx) (S vz)))) (* 2.0 (- (* vy vz) (* s vx))))
          (list (* 2.0 (- (* vx vz) (* s vy))) (* 2.0 (+ (* vy vz) (* s vx))) (- 1.0 (* 2.0 (+ (S vx) (S vy))))))))

(defun quaternion->axis/angle (q)
  (let ((phi (* 2.0 (acos (first q)))))
  (cons (normalize (cdr q))
        phi)))
		
;(def N (expt 10 4))

;(def rot (make-quaternion-rotation (normalize '(-1.0 10.0 1.0)) 0.6))

;(defun random-float () (float (/ (random (expt 10 9))
;                                      (expt 10 9))))

; (def v (build-list N (lambda (j)
                                ; (normalize
                                  ; (build-list 3 (lambda (i) (random-float)))))))

;(time (mapcar (lambda (vec) (apply-quaternion-rotation rot vec)) v))
