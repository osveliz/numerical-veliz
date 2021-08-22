%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized Newton's Method
% Solving Systems of Nonlinear Eq
% If GNU Octave install symbolic package
% youtube.com/OscarVeliz
% @author Oscar Veliz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load symbolic
syms x y F(x,y)
F(x,y) = [x^2-y-1;x-y^2+1]
J = jacobian(F,[x,y])
JX = matlabFunction(J)
invJ = inv(J)
iJF = invJ*F
H = horzcat(diff(F,x,x),diff(F,y,y))
HX = matlabFunction(H)
eps = 10^-6;
Xn = sym([1;2])
while (norm(double(F(Xn(1),Xn(2)))) > eps)
  invJ = inv(JX(Xn(1),Xn(2)));
  a = - invJ * F(Xn(1),Xn(2));
  a2 = a .* a;
  b = 0.5 * invJ * HX(Xn(1),Xn(2)) * a2;
  Xn = Xn + a2 ./ (a + b); #halley
  double(Xn)
end
