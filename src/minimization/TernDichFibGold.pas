(**************************
 * Ternary, Dichotomous,
 * Fibonacci, and Golden-section
 * Search for minimizing f(x)
 * @author Oscar Veliz
 **************************)
program TernDichFibGold;
uses math;
const
    eps = 10E-7; (*10^-6*)
    phi = (1+sqrt(5))/2.0;
    max = 75;

(* A function f(x) to evaluate *)
function f(x: real):real;
begin
    f := power(x,3)/3 - power(x,2)/2 - x - 1;
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
        if f(c) < f(d) then
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
        if f(c) < f(d) then
            b := d
        else
            a := c;
        inc(itr);
    end;
    writeln('iterations ',itr-1);
    dichotomous := m;
end;

(**
 * Fibonacci Search
 * @param a - left side of interval
 * @param b - right side of interval
 * @return midpoint of shrunken [a,b]
 *)
function fibonacci(a,b: real):real;
var
    c,d,fc,fd : real;
    itr, n : integer;
    fib: array[0..max] of qword; (*unsigned 64-bit int*)
begin
	fib[0] := 1;
	fib[1] := 1;
	for n := 2 to max do
	begin
		fib[n] := fib[n-1] + fib[n-2];
		if fib[n] >((b-a)/eps) then
			break;
	end;
	writeln('n=',n,' fibn=',fib[n]);
    c := a + fib[n-2]/fib[n]*(b-a);
    d := a + fib[n-1]/fib[n]*(b-a);    
    fc := f(c);
    fd := f(d);    
    writeln('n,a,c,d,b');
    writeln('0,',a,',',c,',',d,',',b);
    itr := 0;
    for n := n-1 downto 2 do
    begin
        if(fc < fd) then
        begin
            b := d;
            d := c;
            fd := fc;
            c := a + (fib[n-2]/fib[n])*(b-a);
            fc := f(c);
        end
        else
        begin
            a := c;
            c := d;
            fc := fd;
            d := a + (fib[n-1]/fib[n])*(b-a);
            fd := f(d);
        end;
        inc(itr);
        writeln(itr,',',a,',',c,',',d,',',b);
    end;
    writeln('iterations ',itr);
    fibonacci := (c+d)/2;
end;

(**
 * Golden-section Search
 * @param a - left side of interval
 * @param b - right side of interval
 * @return midpoint of shrunken [a,b]
 *)
function goldensection(a,b: real):real;
var
    c,d,fc,fd : real;
    itr : integer;
begin
    c := b + (a-b)/phi;
    d := a + (b-a)/phi;
    fc := f(c);
    fd := f(d);
    itr := 0;
    writeln('n,a,c,d,b');
    writeln(itr,',',a,',',c,',',d,',',b);
    while (abs(b - a) > eps) and (itr < max) do
    begin
        if(fc < fd) then
        begin
            b := d;
            d := c;
            fd := fc;
            c := b + (a-b)/phi;
            fc := f(c)
        end
        else
        begin
            a := c;
            c := d;
            fc := fd;
            d := a + (b-a)/phi;
            fd := f(d);
        end;
        inc(itr);
        writeln(itr,',',a,',',c,',',d,',',b);
    end;
    writeln('iterations ',itr);
    goldensection := (c+d)/2;
end;

(* Main *)
begin
  writeln('ternary ', ternary(1,2));
  writeln('dichotomous ', dichotomous(1,2));
  writeln('fibonacci ', fibonacci(1,2));
  writeln('goldensect ', goldensection(1,2));
end.
