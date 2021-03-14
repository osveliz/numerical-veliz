;;;; Computing pi using Machin-like formula
;;; Example code for the Oscar Veliz YouTube channel
;;; @author Oscar Veliz

(defconstant +digits+ 333)				; 2^-333 ~= 10^-100
(defconstant +itr-cap+ 75)				; only need ~75 itr for 100 places
(SETF (EXT:LONG-FLOAT-DIGITS) +digits+)	; changes accuracy
(defun arctan(x) "arctangent of x defined using Taylor Series"
	(let ((sum 0))
		(loop for i from 0 to +itr-cap+
			for k = (+ (* 2 i) 1)		; k=2i+1
			do (incf sum (* (/ (expt -1 i) k) (expt x k))) ; sum+=(-1^i/k)*x^k
		)
		sum
	)
)
(defun gregory-leibnitz() "Gregory-Leibnitz Method of computing pi"
	(* 4 (arctan 1))					; 4*arctan(1)
)
(defun euler() "Euler's Method of computing pi"
	(* 4 (+ (arctan (/ 1 2)) (arctan (/ 1 3)))) ; 4(arctan(1/2)+arctan(1/3))
)
(defun euler2() "Euler's advanced Method of computing pi"
	(+ (* 20 (arctan(/ 1 7))) (* 8 (arctan(/ 3 79)))) ; 20arctan(1/7) + 8arctan(3/79)
)
(defun machin() "Machin's Method of computing pi"
	(* 4 (- (* 4 (arctan (/ 1 5))) (arctan (/ 1 239)))) ; 4(4arctan(1/5)-arctan(1/239))
)
(defun printer(message number) "Prints a message with coerced number"
	(format t (concatenate 'string message " = ~f~%") (coerce number 'long-float))
)
;;; main
(printer "22/7     " (/ 22 7))
(printer "355/133  " (/ 355 113))
(printer "Greg-Leib" (gregory-leibnitz))
(printer "Euler  Pi" (euler))
(printer "Machin Pi" (machin))
(printer "Euler2 Pi" (euler2))
(printer "Buitin Pi" pi)				; built in pi from CLISP