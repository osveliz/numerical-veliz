/**
 * Demonstration of False Position Method
 * as part of online learning video series on
 * YouTube created by Oscar Veliz.
 * https://www.youtube.com/user/oscarsveliz
 * @author Oscar Veliz
 */

#include <stdio.h>
#include <math.h>
//headers
double f(double);
int sign(double);

int main(){
	double startA = 0, startB = 2;//interval
	double a = startA, b = startB;
	double c = 0, fa = 0, fb = 0, fc = 1000;
	const double eps = 0.0000001;//10^-7
	int i = 0;
	char hbar[90] = {[0 ... 88] = '-'};
	printf("Bisection Method\n");
	printf("a\t\tb\t\tc\t\tf(a)\t\tf(b)\t\tf(c)\n");
	printf("%s\n", hbar);
	while(fabs(fc) > eps){
		c = (a + b) / 2.0;//midpoint
		fa = f(a);
		fb = f(b);
		fc = f(c);
		printf("%lf\t%lf\t%lf\t", a, b, c);
		printf("%lf\t%lf\t%lf\n", fa, fb, fc);
		if(sign(fc) == sign(fa))
			a = c;
		else
			b = c;
		i++;
	}
	printf("%d iterations\n\n", i);
	//reinitialize
	a = startA,	b = startB;
	c = 0, fa = 0, fb = 0, fc = 1000;
	i = 0;
	printf("False Position Method\n");
	printf("a\t\tb\t\tc\t\tf(a)\t\tf(b)\t\tf(c)\n");
	printf("%s\n", hbar);
	while(fabs(fc) > eps){	
		fa = f(a);
		fb = f(b);
		c = b - fb * (b - a) / (fb - fa);//x-intercept
		fc = f(c);
		printf("%lf\t%lf\t%lf\t", a, b, c);
		printf("%lf\t%lf\t%lf\n", fa, fb, fc);
		if(sign(fc) == sign(fa))
			a = c;
		else
			b = c;
		i++;
	}
	printf("%d iterations\n", i);
	return 0;
}

double f(double x){
	return x * x - x - 1;//x^2 - x - 1
}

int sign(double x){
	if(x >= 0)
		return 1;
	return 0;
}
