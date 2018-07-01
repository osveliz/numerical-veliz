(* 
 * Example of Steffensen's Method
 * and Aitken's Delta-squared Method
 * created for the YouTube channel
 * https://youtube.com/OscarVeliz
 * @author Oscar Veliz
 *)

(* Print more of lists in SML/NJ *)
Control.Print.printLength := 100;

(* epsilon *)
val eps = 10E~12;

(* Example Function f(x) *)
fun f x = x*x - x - 1.0;

(* Example g(x) for Fixed Point Iteration *)
fun g x = 1.0 + 1.0 / x;

(* Derivative Approximation (f prime) *)
fun fp x = f(x + f(x))/f(x) - 1.0;

(* Fixed Point Iteration *)
fun fixed x =
	if abs(f(x)) < eps then x
	else fixed(g(x));

(* Fixed Point iteration using a list *)
fun fixedlist x: real list =
	if abs(f(List.last(x))) < eps then x
	else fixedlist(x @ [g(List.last(x))]);

(* Aitken's Delta-squared Method *)
fun aitken(a, b, c) = a - (b-a)*(b-a)/(a-2.0*b+c);

(* Aitken's Delta-squared Method using a list of 3 numbers *)
fun aitkenlist x: real list = [aitken(List.nth(x,0),List.nth(x,1),List.nth(x,2))];

(* Steffensen's Method using Aitken's *)
fun steffaitken x =
  if abs(f(x)) < eps then x
  else steffaitken(aitken(x, g(x), g(g(x))));

(* Steffensen's Method using Aitken's showing each step *)
fun showsteffaitken x =
  let
    val p : real ref = ref x 
    val a : real ref = ref x
    val b : real ref = ref x
    val c : real ref = ref x
    val n : int ref = ref 1
  in
	print("n\ta\t\tb\t\tc\t\tp\n");
    while abs(f(!p))> eps andalso !n < 1000 do(
	  a := !p;
	  b := g(!a);
	  c := g(!b);
	  p := aitken(!a,!b,!c);
      print((Int.toString(!n)) ^ "\t" ^ (Real.toString(!a)) ^ "\t" ^ (Real.toString(!b)) ^ "\t" ^ (Real.toString(!c)) ^ "\t" ^ (Real.toString(!p)) ^ "\n");
      n := !n + 1
    );
    !p
end;

(* Steffensen's Method using Aitken's with a list *)
fun steffaitkenlist x =
  if abs(f(List.last(x))) < eps then x
  else steffaitkenlist(x @ aitkenlist([List.last(x), g(List.last(x)), g(g(List.last(x)))]));

(* Steffensen's Method with derivative *)
fun steff x =
	if abs(f(x)) < eps then x
	else steff(x - f(x)/fp(x));

(* Steffensen's Method with derivative using a list *)
fun stefflist x: real list =
	if abs(f(List.last(x))) < eps then x
	else stefflist(x @ [(List.last(x) - f(List.last(x))/fp(List.last(x)))]);

(* main *)
val start = 2.0;

print("\nFixed Point Iteration\n");
fixedlist([start]);

print("\nSteffensen's Method using Fixed Point Iteration and Aitken's Delta-squared Method\n");
showsteffaitken(start);
steffaitkenlist([start]);

print("\nSteffensen's Method using Derivative Approximation\n");
stefflist([start]);
