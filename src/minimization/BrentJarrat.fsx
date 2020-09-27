(*
    Brent's Minimization Method & Jarrat's Method
    @author Oscar Veliz
*)
[<Literal>]
let Max = 100
let (eps : float) = 10.0**(-7.0)
let invphi = (1.0 + sqrt(5.0)) / 2.0 - 1.0
let ratio = (3.0 - sqrt(5.0)) / 2.0
let lookup = new System.Collections.Generic.Dictionary<_,_>()

(*
    Function f(x) that saves results to avoid recomputing
    @param x
    @return f(x)
*)
let f (x : float) : float =
    let key , value = lookup.TryGetValue(x)
    if key then 
        value
    else
        let value = x**3.0 / 3.0 - x**2.0 / 2.0 - x - 1.0
        lookup.Add(x, value)
        value

(*
    Single iteration of Jarrat's Method
    @parameters x_(n-2), x_(n-1), x_n
    @return x_(n+1)
*)
let jarrat (xnm2 : float) (xnm1 : float) (xn : float) : float = 
    let fn = f(xn)
    let fnm1 = f(xnm1)
    let fnm2 = f(xnm2)
    xn + 0.5 * ((xnm1 - xn)**2.0 * (fn - fnm2) + (xnm2 - xn)**2.0 * (fnm1 - fn)) / ((xnm1 - xn) * (fn - fnm2) + (xnm2 - xn) * (fnm1 - fn))

(*
    One iteration of Jarrat's Method using a list
    @param [x_(n-2); x_(n-1); x_n]
    @return [x_(n-2); x_(n-1);  x_n; x_(n+1)]
*)
let jarratList (x : List<float>) : List<float> =
    let n = (x.Length)
    x @ [(jarrat (x.Item(n - 3)) (x.Item(n - 2)) (x.Item(n - 1)))]

(*
    Jarrat's Method iterating until xmin is found
    @param [x_0; x_1; x_2]
    @return [x_0; x_1; x_2; ...; x_n]
*)
let rec polynomialInterpolation (x : List<float>) : List<float> =
    let n = (x.Length)
    if abs(x.Item(n - 1) - x.Item(n - 2)) < eps || n > Max then x
    else polynomialInterpolation (jarratList x)

(*
    Recursive loop of Golden-section Search
    @param a left interval
    @param c first test point
    @param d second test point
    @param b right interval
    @return minimum
*)
let rec golden (a : float) (c : float) (d : float) (b : float)  n: float =
    if b - a < eps then
        (a + b) / 2.0
    elif f(c) < f(d) then
        let t = (d + invphi * (a - d))
        //printfn "%d %.10f %.10f %.10f %.10f" n a t c d // Uncomment to see itr
        golden a t c d (n + 1)
    else
        let t = (c + invphi * (b - c))
        //printfn "%d %.10f %.10f %.10f %.10f" n c d t b // Uncomment to see itr
        golden c d t b (n + 1)
(*
    Golden-section Search Initializer
    @param a left interval
    @param b right interval
    @return minimum
*)
let goldenSection (a : float) (b : float) : float =
    let c = b + invphi * (a - b)
    let d = a + invphi * (b - a)
    //printfn "%d %.10f %.10f %.10f %.10f" 0 a c d b // Uncomment to see itr
    golden a c d b 1

(*
    Brent's Minimization Method
    Coming Soon
*)

// Main
printfn "Jarrat's Method \n%A" (polynomialInterpolation [2.0; 1.0; 1.5])
printfn "Golden-section Search \n%.10f" (goldenSection 1.0 2.0)