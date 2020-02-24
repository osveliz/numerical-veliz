#################################
# GNU Octave program plotting
# Taylor and MacLaurin Series
# created for the YouTube channel
# https://youtube.com/OscarVeliz
# @author Oscar Veliz
#################################
clear all;
clf;
pkg load symbolic

################################
# Plot Taylor Series
# @param f function string
# @param a point approximating around
# @param n number of terms in the series
#################################
function plotTaylor(f,a,n)
  syms x;
  fx = function_handle(sym(f));
  limits = [a-2,a+2,fx(a)-2,fx(a)+2];
  fplot(fx,limits,1000,'k');
  ylabel("y");
  xlabel("x");
  title(strcat("Taylor Series n= ",num2str(n)," a=",num2str(a)));
  derivs = cell(1,n);
  derivs(1) = char(simplify(diff(sym(f))));
  for i = 2:n
    derivs(i) = char(simplify(diff(sym(derivs(i-1)))));
  endfor
  derivs
  tn = num2str(fx(a));
  for i = 1:n
    fp = function_handle(sym(derivs(i)));
    tn = strcat(tn,"+",num2str(fp(a)),"/(",num2str(i),"!)*(x-",num2str(a),")**",num2str(i));
  endfor
  tn = char(simplify(sym(tn)))
  tnx = function_handle(sym(tn));
  hold on;
  fplot(tnx,limits, '--r')
  legend(f,tn, "location", "southoutside");
  print -dsvg "taylorplot.svg";
endfunction

################################
# Plot MacLaruin Series
# @param f function string
# @param n number of terms in the series
#################################
function plotMacLaurin(f,n)
  plotTaylor(f,0,n)
endfunction

###Main
func = "atan(x)"
plotTaylor(func,0,5);
