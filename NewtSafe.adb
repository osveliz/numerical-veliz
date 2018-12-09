-----------------------------------------
-- Example Newton-Bisection Hybrid
-- Root-Finding Method created by
-- Oscar Veliz based off of NewtSafe from 
-- https://www.ldeo.columbia.edu/~mspieg/e4300/BlankPDFs/Lecture06_blank.pdf
-- and rtsafe from Numerical Recipies in C p.366-367 and
-- Numerical Methods That Work by Acton p.51-52
-- www.youtube.com/OscarVeliz
-----------------------------------------
with Ada.Text_IO; use  Ada.Text_IO;
with Ada.Integer_Text_IO; use  Ada.Integer_Text_IO;
with Ada.Float_Text_IO; use  Ada.Float_Text_IO;
with Ada.Task_Identification;  use Ada.Task_Identification;
with Ada.Numerics.Generic_Elementary_Functions;

procedure NewtSafe is
	---------------------
	-- f(x) = x^3 - x^2 - x - 1
	---------------------
	function f (x: Float) return Float is
	begin
		return x**3 - x**2 - x - 1.0;
	end;
	-------------------
	-- f'(x) = 3x^2 - 2x - 1
	-------------------
	function fp (x: Float) return Float is
	begin
		return 3.0*x*x - 2.0*x - 1.0;
	end;
	----------
	-- Regular Newton's Method (unsafe)
	----------
	function newt (xinit: Float) return Float is
		x: Float := xinit;
		eps: Float := 10.0**(-6);
	begin
		while abs(f(x)) > eps loop
			put(x,2,6,0);put_line("");
			x := x - f(x) / fp(x);
		end loop;
		return x;
	end;
	-- header so forceBis can call hybrid
	function hybrid (inita, initb: Float) return Float;
	----------
	-- Bisection between a and b then restart Hybrid
	----------
	function forceBis(a, b: Float) return Float is
		c: Float := (a + b) / 2.0;
	begin
		--put_line("forced bisection");
		if f(a) * f(c) < 0.0 then -- fa & fc have different signs
			return hybrid(a, c);
		else
			return hybrid(c, b);
		end if;
	end;
	----------
	-- Newton-Bisection Hybrid, inita and initb specify interval
	----------
	function hybrid (inita, initb: Float) return Float is
		x, fa, fb, fx, fpx, oldx, lastStep: Float;
		a: Float := inita;
		b: Float := initb;
		eps: Float := 10.0**(-6);
	begin
		if(a > b) then --ensure a is on left, swap
			a := initb;
			b := inita;
		end if;
		if a = b then
			return a;
		end if;
		-- check interval points
		fa := f(a);
		fb := f(b);
		if fa = 0.0 then
			return a;
		end if;
		if fb = 0.0 then
			return b;
		end if;
		if fa * fb > 0.0 then -- invalid interval
			put_line("error: interval must contain root.");
			Abort_Task (Current_Task);
		end if;
		lastStep := b - a;
		oldx := b; -- intitial oldx is on boundary
		x := (a + b) / 2.0; -- start from midpoint
		fx := f(x);
		while abs(fx) > eps loop
			put(x,2,6,0);put_line("");
			fpx := fp(x);
			if fpx = 0.0 then -- avoid divide by zero
				return forceBis(a,b);
			end if;
			lastStep := abs(x - oldx);
			if abs(fx * 2.0) > abs(lastStep * fpx) then -- decreasing slowly
				--put_line("too slow");
				return forceBis(a,b);
			end if;
			oldx := x;
			x := x - fx / fpx;
			if x = oldx then -- converged
				return x;
			end if;
			if x <= a or else x >= b then -- avoid going outside interval
				--put_line("outside interval");
				return forceBis(a,b);
			end if;
			--shrink interval
			fx := f(x);			
			if fx * fa < 0.0 then
				b:= x;
				fb:= fx;
			else
				a:= x;
				fa:= fx;
			end if;
		end loop;
		return x;
	end;
	begin --main--
		put_line("Newton's Method");
		put(newt(1.5),2,6,0); put_line("");
		put_line("Newt-Safe");
		put(hybrid(1.0,2.0),2,6,0); put_line("");
end NewtSafe;
