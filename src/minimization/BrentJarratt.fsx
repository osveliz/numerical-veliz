(*
    Brent's Minimization Method with
    Jarratt's Method and Golden-section Search
    @author Oscar Veliz
*)
[<Literal>]
let Max = 100
let (eps : float) = 10.0**(-6.0)
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
    Single iteration of Jarratt's Method
    @parameters x_(n-2), x_(n-1), x_n
    @return x_(n+1)
*)
let jarratt (xnm2 : float) (xnm1 : float) (xn : float) : float = 
    let fn = f(xn)
    let fnm1 = f(xnm1)
    let fnm2 = f(xnm2)
    xn + 0.5 * ((xnm1 - xn)**2.0 * (fn - fnm2) + (xnm2 - xn)**2.0 * (fnm1 - fn)) / ((xnm1 - xn) * (fn - fnm2) + (xnm2 - xn) * (fnm1 - fn))

(*
    One iteration of Jarratt's Method using a list
    @param [x_(n-2); x_(n-1); x_n]
    @return [x_(n-2); x_(n-1);  x_n; x_(n+1)]
*)
let jarrattList (x : List<float>) : List<float> =
    let n = (x.Length)
    x @ [(jarratt (x.Item(n - 3)) (x.Item(n - 2)) (x.Item(n - 1)))]

(*
    Jarratt's Method iterating until xmin is found
    @param [x_0; x_1; x_2]
    @return [x_0; x_1; x_2; ...; x_n]
*)
let rec polynomialInterpolation (x : List<float>) : List<float> =
    let n = (x.Length)
    if abs(x.Item(n - 1) - x.Item(n - 2)) < eps || n > Max then x
    else polynomialInterpolation (jarrattList x)

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
    Recursive loop of Brent's Minimization Method
    @param a left interval
    @param b right interval
    @param v previous iterate value
    @param w previous iterate value
    @param x previous iterate value
    @param dold last delta step
    @param eold last golden interval size
    @param i iteration counter
    @return minimum
*)
let rec brent (a : float) (b : float) (v : float) (w : float) (x : float) (dold : float) (eold : float) (i : int) : float =
    let fv = f(v)
    let fw = f(w)
    let fx = f(x)
    let newi = i + 1
    let m = 0.5 * (a + b)
    if b - a <= eps || i > Max then
        m
    else
        let r = (x - w) * (fx - fv)
        let tq = (x - v) * (fx - fw)
        let tp = (x - v) * tq - (x - w) * r
        let tq2 = 2.0 * (tq - r)
        let p = if tq2 > 0.0 then -tp else tp
        let q = if tq2 > 0.0 then tq2 else -tq2
        let safe = q <> 0.0
        let deltax = if safe then p / q else 0.0
        let parabolic = safe && a < x + deltax && x + deltax < b && abs(deltax) < 0.5 * abs(eold)
        let e = if parabolic then dold elif x < m then b - x else a - x
        let d = if parabolic then deltax else ratio * e
        let u = x + d
        let fu = f(u)
        //printfn "i=%d\ta=%f\tb=%f\tx=%f\tw=%f\tv=%f\tu=%f" i a b x w v u //uncomment to see each variable
        //printfn "i=%d\tu=%f\tfu=%f\t%s" newi u fu (if parabolic then "parabolic" else "golden") //uncomment to see itr
        if fu <= fx then
            let newa = if u < x then a else x
            let newb = if u < x then x else b
            brent newa newb w x u d e newi
        else
            let newa = if u < x then u else a
            let newb = if u < x then b else u
            if fu <= fw || w = x then brent newa newb w u x d e newi
            elif fu <= fv || v = x || v = w then brent newa newb u w x d e newi
            else brent newa newb v w x d e newi

(*
    Brent's Minimization Method Initializer
    @param a left interval
    @param b right interval
    @return minimum
*)
let brentMethod (a: float) (b : float) : float =
    let x = b + invphi * (a - b)
    brent a b x x x 0.0 0.0 0

// Main
printfn "Jarratt's Method \n%A" (polynomialInterpolation [2.0; 1.0; 1.5])
printfn "Golden-section Search \n%.10f" (goldenSection 1.0 2.0)
printfn "Brent's Minimizatin Method \n%.10f" (brentMethod 1.0 2.0)