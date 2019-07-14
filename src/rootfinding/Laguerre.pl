##################################
# Laguerre's Method for solving polynomials
# Created for the YouTube channel
# https://youtube.com/OscarVeliz
# @author Oscar Veliz
##################################
use Math::Complex;

##################################
# Find the derivative of P(x)
# @param a array (-1,0,1) = x^2-1
# @return coefficients for P'(x)
##################################
sub powerRule {
    my ($a) = @_;
    my @coef = @{$a};
    my @deriv;
    for (my $i = 1; $i < @coef; $i++){
        push(@deriv,$coef[$i] * $i);
    }
    return \@deriv;
}

##################################
# Horner's Method computing P(x)
# @param a array (-1,0,1) = x^2-1
# @param x the value for input
# @return P(x)
##################################
sub horner {
    my ($a,$x) = @_;
    my @coef = @{$a};
    my $result = $coef[@coef-1];
    for(my $i = @coef - 2; $i >= 0; $i--){
        $result = $coef[$i] + $x * $result;
    }
    return $result;
}

##################################
# Maximum in absolute value
# @param a a number
# @param b another number
# @return a or b depending on abs
##################################
sub maxabs {
    my ($a,$b) = @_;
    if(abs($a) > abs($b)){
        return $a;
    }
    else{
        return $b;
    }
}

##################################
# Laguerre's Method
# @param a array (-1,0,1) = x^2-1
# @param x starting value
# @return a root of P(x)
##################################
sub laguerre {
    my ($a,$x) = @_;
    my ($G,$R,$H,$sr);
    my @coef = @{$a};
    my @dcoef = @{powerRule(\@coef)};
    my @ddcoef = @{powerRule(\@dcoef)};
    my $eps = 10e-12;
    my $degree = @coef-1;
    my $k = 0;
    my $max = 100;
    my $p = horner(\@coef,$x);
    while(abs($p) > $eps && $k < $max){
        print "x${k} = ${x}\n";
        $G = horner(\@dcoef,$x) / $p;
        $H = $G**2 - horner(\@ddcoef,$x) / $p;
        $sr = sqrt(($degree - 1)*($degree * $H - $G**2));
        $x = $x - $degree / maxabs($G - $sr, $G + $sr);
        $p = horner(\@coef,$x);
        $k++;
    }
    print "x${k} = ${x}\n";
    if($k == $max){#never expect to enter this if
        print "maximum iteration count reached\n"
    }
    return $x;
}

### main ###
@coefficents = (-1, -1, -1/2, 1/3);
$start = 0;
$root = laguerre(\@coefficents,$start);
