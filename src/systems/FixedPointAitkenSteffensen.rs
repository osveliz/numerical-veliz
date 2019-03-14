/**
 * Fixed Point Iteration for Systems of Equations
 * with Generalized Aitken-Steffensen Method
 * Created by Oscar Veliz
 * youtube.com/OscarVeliz
 */
const EPS : f64 = 10e-7;//terminal epsilon
const MAX : usize = 100;//maximum number iterations
const VARS: usize = 2;  //number of unknown variables
const ZERO: [f64;VARS] = [0.0;VARS];//zeros

/**
 * Call fpi with starting values and  false for normal Fixed Point Iteration
 * Use true for applying Generalized Aitken-Steffensen
 */
fn main() {
	let start: [f64;VARS] = [1.0, 2.0];
    fpi(start,false);//fixed point iteration
    fpi(start,true);//aitken-steffensen
}

/**
 * G(X)
 * @param s the metric space
 * @return iteration on space
 */
fn g(s:[f64;VARS]) -> [f64;VARS]{
	let mut xn: [f64;VARS] = ZERO;
	//the following variables make equations easier to read
	let x = s[0];
	let y = s[1];
	xn[0] = (y + 1.0).sqrt(); //x = sqrt(y + 1)
	xn[1] = (x + 1.0).sqrt(); //y = sqrt(x + 1)
	
	//xn[0] = -(y + 1.0).sqrt();//x = -sqrt(y + 1)
	//xn[1] = -(x + 1.0).sqrt();//y = -sqrt(x + 1)
	
	//xn[0] = (y + 1.0) / x; //x = (y + 1)/x
	//xn[1] = (x + 1.0) / y; //y = (x + 1)/y
	
	//these equations can not be combined with the previous ones
	//xn[0] = y*y - 1.0; //x = y^2 - 1
	//xn[1] = x*x - 1.0; //y = x^2 - 1
	return xn;
}

/**
 * F(X) = 0
 * @param s metric space
 * @return F applied to metric space
 */
fn f(s:[f64;VARS]) -> [f64;VARS]{
	let mut xn: [f64;VARS] = ZERO;
	//the following variables make the equations easier to read
	let x = s[0];
	let y = s[1];
	xn[0] = x - y*y + 1.0; //x - y^2 = -1
	xn[1] = x*x - y - 1.0; //x^2 - y = 1
	return xn;
}

/**
 * Fixed Point Iteration or Steffensen's Method when useaitken is true
 * @param s starting values
 * @param useaitken false = normal fixed point, true = steffensen's method
 * @return X* solution
 */
fn fpi(s:[f64;VARS], useaitken:bool) -> [f64;VARS]{
	if useaitken{ println!("Aitken-Steffensen"); }
	else{ println!("Fixed Point Iteration"); }
	let mut x:[[f64;VARS];MAX] = [ZERO;MAX];
	x[0] = s;
	let mut i:usize = 0;
	let mut done = false;
	println!("X0 = {:?}",x[0]);
	while !done{
		if useaitken{ x[i+1] = aitken(x[i]); }
		else{ x[i+1] = g(x[i]);	}
		i += 1;
		println!("X{} = {:?}", i,x[i]);
		if norm(f(x[i])) < EPS || d(x[i], x[i-1]) < EPS || i >= MAX-1 {
			done = true;
		}
	}
	println!("X*= {:?}" ,x[i]);
	println!("final error = {:?}",norm(f(x[i])));
	println!("{} iterations",i);
	println!("order = {}", (norm(f(x[i]))/norm(f(x[i-1]))).ln()/(norm(f(x[i-1]))/norm(f(x[i-2]))).ln());
	return x[i];
}

/**
 * Vector subtraction A - B
 * @param a (n element vector)
 * @param b (n element vector)
 * @return a - b
 */
fn minus(a:[f64;VARS], b:[f64;VARS]) -> [f64;VARS]{
	let mut s: [f64;VARS] = ZERO;
	for i in 0..VARS{
		s[i] = a[i] - b[i];
	}
	return s;
}
/**
 * Matrix subtraction A - B
 * @param a (nxn matrix)
 * @param b (nxn matrix)
 * @return a - b
 */
fn minusnxn(a:[[f64;VARS];VARS], b:[[f64;VARS];VARS]) -> [[f64;VARS];VARS]{
	let mut s:[[f64;VARS];VARS] = [ZERO;VARS];
	for i in 0..VARS{
		for j in 0..VARS{
			s[i][j] = a[i][j] - b[i][j];
		}
	}
	return s;
}	
/**
 * Distance between two spaces ||A - B||_2
 * @param a space1
 * @param b space2
 * @return distance
 */
fn d(a:[f64;VARS], b:[f64;VARS]) -> f64{
	return norm(minus(a,b));
}

/**
 * l^2 norm of a space = sqrt(sum(s[i]^2))
 * @param s the input space
 * @return the l^2 norm
 */
fn norm(s:[f64;VARS]) -> f64{
	let mut sum: f64 = 0.0;
	for i in 0..VARS{
		sum += (s[i]).powi(2);//s[i]^2
	}
	return sum.sqrt();
}

/**
 * Combines two vectors into side-by-side matrix
 * @param a vector 1
 * @param b vector 2
 * @return [a,b]
 */
fn lateral(a:[f64;VARS], b:[f64;VARS]) -> [[f64;VARS];VARS]{
	let mut x:[[f64;VARS];VARS] = [ZERO;VARS];
	for i in 0..VARS{
		x[i][0] = a[i]
	}
	for i in 0..VARS{
		x[i][1] = b[i]
	}
	return x;
}
/**
 * Perform one iteration of generalized Aitken's Delta squared
 * @param s starting vector
 * @return X-hat
 */
fn aitken(s:[f64;VARS]) -> [f64;VARS]{
	let a:[f64;VARS] = s;
	let b:[f64;VARS] = g(a);
	let c:[f64;VARS] = g(b);
	let d:[f64;VARS] = g(c);
	let bminusa:[f64;VARS] = minus(b,a);
	let cminusb:[f64;VARS] = minus(c,b);
	let dminusc:[f64;VARS] = minus(d,c);
	let dx1:[[f64;VARS];VARS] = lateral(bminusa,cminusb);//delta x1
	let dx2:[[f64;VARS];VARS] = lateral(cminusb,dminusc);//delta x2
	let d2x:[[f64;VARS];VARS] = minusnxn(dx2,dx1);//delta squared x
	let det:f64 = 1.0 / (d2x[0][0]*d2x[1][1] - d2x[0][1]*d2x[1][0]);//determinant
	let inverse:[[f64;VARS];VARS] = [[det*d2x[1][1],-det*d2x[0][1]],[-det*d2x[1][0],det*d2x[0][0]]];//inverse delta squared x
	let mut nxn:[[f64;VARS];VARS] = [ZERO;VARS];//deltax1 * inverse
	for i in 0..VARS{
		for j in 0..VARS{
			for k in 0..VARS{
				nxn[i][j] += dx1[i][k] * inverse[k][j]
			}
		}
	}
	let mut temp:[f64;VARS] = ZERO;//(deltax1 * inverse) * (b-a)
	for i in 0..VARS{
		for j in 0..VARS{
			temp[i] += nxn[i][j]*bminusa[j];
		}
	}
	return minus(a,temp); // X-hat = X_n - DeltaX_n*(DeltaSquaredX_n^-1)*(X_(n+1) - X_n)
}
