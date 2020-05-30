(**************************
 * Ternary, Dichotomous,
 * Fibonacci, and Golden-section
 * Search for minimizing f(x)
 * @author Oscar Veliz
 **************************)
program TernDichFibGold;
const
    eps = 10E-7;
    phi = (1+sqrt(5))/2.0;
    max = 75;
(* A function f(x) to evaluate *)
function f(x: real):real;
begin
    f:=x*x*x/3-x*x/2-x-1.0;
end;
(**
 * Ternary Search
 * @param a - left side of interval
 * @param b - right side of interval
 * @return midpoint of shrunken [a,b]
 *)
function ternary(a,b: real):real;
var
    c,d : real;
    itr : integer;
begin
    itr := 0;
    writeln('n,a,c,d,b');
    while (abs(b - a) > eps) and (itr < max) do
    begin
        c := (2.0*a+b)/3.0;
        d := (a+2.0*b)/3.0;
        writeln(itr,',',a,',',c,',',d,',',b);
        if f(c) < f(d)  then
            b := d
        else
            a := c;
        inc(itr);
    end;
    writeln('iterations ',itr-1);
    ternary := (a+b)/2;
end;

(**
 * Dichotomous Search
 * @param a - left side of interval
 * @param b - right side of interval
 * @return midpoint of shrunken [a,b]
 *)
function dichotomous(a,b: real):real;
var
    c,d,m,halfeps : real;
    itr : integer;
begin
	halfeps := eps / 2.0;
    itr := 0;
    writeln('n,a,c,d,b');
    while ( abs(b - a)/2.0 > eps) and (itr < max) do
    begin
		m := (a+b)/2.0;
        c := m - halfeps;
        d := m + halfeps;
        writeln(itr,',',a,',',c,',',d,',',b);
        if f(c) < f(d)  then
            b := d
        else
            a := c;
        inc(itr);
    end;
    writeln('iterations ',itr-1);
    dichotomous := m;
end;

(* Coming soon Fibonacci and Golden-section Search *)

(* Main *)
begin
  writeln('ternary ', ternary(1,2));
  writeln('dichotomous ', dichotomous(1,2));
end.
