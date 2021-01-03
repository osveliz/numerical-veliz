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
invJ = inv(J)
iJF = invJ*F
eps = 10^-7;
Xn = sym([1;1])
while (norm(double(F(Xn(1),Xn(2)))) > eps)
  Xn = Xn - iJF(Xn(1),Xn(2))
end
