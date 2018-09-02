(* Power Method and Inverse Power Method
 * to find eigenvalues and eigenvectors
 * (for simplicity arrays represent vectors)
 * @author Oscar Veliz
 *)
open Printf

(*print an array using a loop*)
let println a = for i = 0 to (Array.length (a)) do
					if i < Array.length (a) then
						printf "%f " a.(i)
					else printf "\n"
				done ;;

(*print a matrix*)
let printMatrix m = Array.iter (println) m;;

(*couldn't use create_float for some reson*)
let createFloatArray n = let x = Array.make_matrix 1 n 0.0 in x.(0);;

(*create a nxn identity matrix*)
let identity n =	let m = Array.make_matrix n n 0.0 in
					for i = 0 to (n -1) do
						m.(i).(i) <- 1.0
					done;
					m;;

(*multiply matrix m by vector b store result vector in x*)
let matrixByVector m b x = for i = 0 to (Array.length m) -1 do
								x.(i) <- 0.0;
								for j = 0 to (Array.length x) -1 do
									x.(i) <- x.(i) +. m.(i).(j) *. b.(j)
								done
							done;;
(*same as matrixByVector except returns result*)
let matrixByVector2 m b =	let n = ref (Array.length b)in
							let x = createFloatArray !n in
							for i = 0 to (Array.length m) -1 do
								for j = 0 to (Array.length x) -1 do
									x.(i) <- x.(i) +. m.(i).(j) *. b.(j)
								done
							done;
							x;;

(*vector transpose times vector returns number 1xn by nx1 = 1x1*)
let vtByV vt v = 	let n = ref (Array.length v) in
					let sum = ref 0.0 in
					for i = 0 to !n - 1 do
						sum := !sum +. vt.(i) *. v.(i);
					done;
					!sum;;

(*overwrite array with random values*)
let randomize x = for i = 0 to (Array.length x) -1 do 
					x.(i) <- Random.float 10.0
					done;;

(*find maximum value in vector x*)
let max x = let m = ref (x.(0)+.0.0) in
					for i = 1 to (Array.length x) -1 do
						if !m < x.(i) then
							m := x.(i)
					done;
					!m;;
(*create a new matrix where every value is multiplied by s*)
let scaler s x =	let n = ref (Array.length x) in
					let m = Array.make_matrix !n !n 0.0 in
					for i = 0 to (!n -1) do
						for j = 0 to (!n -1) do
							m.(i).(j) <- s *. x.(i).(j)
						done
					done;
					m;;

(*alter array to multiply every value by s*)
let multScalar x s =	for i = 0 to (Array.length x) - 1 do
							x.(i) <- x.(i) *. s
						done;;

(*alter array to divide every value by s*)
let divideScalar x s = let d = 1.0 /. s in multScalar x d;;

(*return result of two matrix a minus matrix b*)
let sub a b = 	let n = ref (Array.length a) in
				let x = Array.make_matrix (!n) (!n) 0.0 in
				for i = 0 to !n -1 do
					for j = 0 to !n -1 do
						x.(i).(j) <- a.(i).(j) -. b.(i).(j)
					done;
				done;
				x;;

(*store the values from vector b into vector a*)
let copyTo a b = for i = 0 to Array.length a -1 do a.(i) <- b.(i) done;;

(*determinant of 2x2 matrix*)
let det2x2 m = m.(0).(0) *. m.(1).(1) -. m.(0).(1) *. m.(1).(0);;

(*invert 2x2 matrix a into new matrix*)
let inv2x2 a = 	let d = det2x2 a in
				let b = Array.make_matrix 2 2 0.0 in
				b.(0).(0) <- a.(1).(1) /. d;
				b.(0).(1) <- 0.0 -. a.(0).(1) /. d;
				b.(1).(0) <- 0.0 -. a.(1).(0) /. d;
				b.(1).(1) <- a.(0).(0) /. d;
				b;;

(*compute Rayleigh quotient*)
let rayleigh a v = vtByV v (matrixByVector2 a v) /. (vtByV v v);;

(*power method*)
let power a b = let m = ref 0.0 in
				let n = ref (Array.length a)in
				let x = createFloatArray !n in
				printf "\nApply Power to\n";
				printMatrix a;
				for i = 0 to 10 do
					printf "b%d " i;
					println b;
					matrixByVector a b x;
					printf "x ";
					println x;
					m := max x;
					printf "max %f\n" !m;
					divideScalar x !m;
					copyTo b x;
				done;;

(*inverse power method*)
let invPower a b =	let matrix = inv2x2 a in
					printf "\nApply Inverse Power to A\n";
					printMatrix a;
					printf "Inverse\n";
					printMatrix matrix;
					power matrix b;
					printf "lambda %f\n" (rayleigh a b);;

(*inverse power method with shift*)
let invPowerShift a b mu =  let n = ref (Array.length a)in
							let ident = identity !n in
							let muI = scaler mu ident in
							let tempM = sub a muI in
							let matrix = inv2x2 tempM in
							printf "\nApply Inverse Power to A\n";
							printMatrix a;
							printf "Identity\n";
							printMatrix ident;
							printf "muI\n";
							printMatrix muI;
							printf "A - muI\n";
							printMatrix tempM;
							printf "Inverse\n";
							printMatrix matrix;
							power matrix b;
							printf "lambda %f\n" (rayleigh a b);;

let a = [|[|0.0;1.0|];
		  [|1.0;1.0|]|];;			  
let b = [| 1.0; 1.0|];;
power a b;;
(*use if power method doesn't quite converge*)
printf "lambda %f\n" (rayleigh a b);;
let b2 = [| 1.0; 1.0|];; 
(*invPower a b2;;*)
invPowerShift a b2 (-0.5);;
