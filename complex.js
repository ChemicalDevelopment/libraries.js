/*

Complex functions

*/

var pi = 3.141592653589793238;

function re(a) {
  return [a[0], 0];
}

function im(a) {
  return [0, a[1]];
}


function neg(a) {
  return [-a[0], -a[1]];
}


// 1 / a
function rec(a) {
  return scale(conj(a), 1 / (a[0] * a[0] + a[1] * a[1]));
}

//Conjugate
function conj(a) {
  return [a[0], - a[1]];
}


//Angle
function ang(a) {
  return Math.atan2(a[1], a[0]);
}


//Norm
function norm(a) {
  return Math.sqrt(norm_sqr(a));
}

//Norm * Norm
function norm_sqr(a) {
  return a[0] * a[0] + a[1] * a[1];
}

//scales a by b
function scale(a, b) {
  return [a[0] * b, a[1] * b];
}

//Adds a and b
function add(a, b) {
  return [a[0] + b[0], a[1] + b[1]];
}

//a - b
function sub(a, b) {
  return [a[0] - b[0], a[1] - b[1]];
}

//a * b
function mul(a, b) {
  return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
}

//a * a
function sqr(a) {
  return [a[0] * a[0] - a[1] * a[1], 2 * a[0] * a[1]];
}

// a/b
function div(a, b) {
  return scale(mul(a, conj(b)), (1 / norm_sqr(b)));
}

//e^a
function exp(a) {
  return scale([Math.cos(a[1]), Math.sin(a[1])], Math.exp(a[0]));
}

//ln(a)
function log(a) {
  var real = Math.log(norm(a));
  var imag = Math.atan2(a[1], a[0]);
  return [real, imag];
}

//a^b
function pow(a, b) {
  var lga = log(a);
  var _e = mul(lga, b);
  return exp(_e);
}
//a^b, b is an integer
function pow_i(a, b) {
    if (b == 0) {
      return 1;
    }
    if (b == 1) {
      return a;
    }
    if (b % 2 == 0) {
      return pow_i(sqr(a), Math.floor(b / 2));
    }
    if (b % 2 == 1) {
      return mul(a, pow_i(sqr(a), Math.floor((b - 1) / 2)));
    }
}


/*

Functions class. Implements basic functions as well as some lesser known ones

*/

function dragon_geometric(x) {
  var sum = [0, 0];
  var sign = 1;
  for (var i = 1; i < 35; ++i) {
      sum = add(sum, scale(pow([i, 0], sub(x, [i, 0])), sign));
      sign = sign * -1;
  }
  return sum;
}


/*

Riemann zeta function. Sum of the recipricols to some power x

*/
/*var cheb_term = [0.1, 20.1, 680.1, 9128.1, 64040.1, 269044.89999999, 734964.899999999999, 1390324.9, 1947380.9, 2209524.9, 2261953.699999];
var cheb_n = 10;*/

var cheb_term = [0.05, 40.05, 5360.05, 286256.05, 8131280.050000001, 142019689.65, 1663478889.65, 13835152489.65, 85039443049.65, 397779856489.6499, 1447929244777.65, 
4175589993577.6494, 9690208463977.648, 18377853561961.65, 28996086459497.65, 38955256625462.445, 45982896863542.445, 49590669392182.445, 50861979711798.445, 51136857618742.445, 51164345409436.84];
var cheb_n = 20;

function init_cheb() {

}


function zeta(x) {
  if (x[0] < 0) {
    var ref_x = sub([1, 0], x);
    var ze_ref = zeta(ref_x);
    var c_xp = sin(scale(x, pi / 2));
    var g_ref = gamma(ref_x);
    var tpi_tx = pow([2 * pi , 0], x);
    var a = mul(ze_ref, g_ref);
    var b = mul(c_xp, tpi_tx);
    return scale(mul(a, b), 25.0 / 78.0);
  } 
  var sum = [0, 0];
  var sign = 1;
  var _c = sub([1, 0], pow([2, 0], sub([1,0], x)));
  for (var i = 0; i <= cheb_n - 1; ++i) {
    var tmp = pow([i + 1, 0], x);
    var scl = sign * (cheb_term[i] - cheb_term[cheb_n]);
    var ot = div([1, 0], tmp);
    sum = add(sum, scale(ot, scl));
    sign *= -1;
  }
  sum = div(sum, _c);
  sum = scale(sum, - 1 / cheb_term[cheb_n]);
  return sum;
}


function gamma(x) {
  if (x[0] < 4) {
    var xp1 = add(x, [1, 0]);
    return div(gamma(xp1), x);
  }
  var lgx = log(x);
  var x_n2 = sqr(div([1, 0], x));
  var x_k_n1 = div([1, 0], x);
  var _ex = mul(x, lgx);
  _ex = sub(_ex, x);
  _ex = sub(_ex, scale(log(scale(x, 1 / (2 * 3.141592653589793238))), 1 / 2));
  _ex = add(_ex, scale(x_k_n1, 1 / 12));
  x_k_n1 = mul(x_k_n1, x_n2);
  _ex = sub(_ex, scale(x_k_n1, 1 / 360));
  x_k_n1 = mul(x_k_n1, x_n2);
  _ex = add(_ex, scale(x_k_n1, 1 / 1260));
  return exp(_ex);
}


function sin(x) {
  var e_ix = exp(mul(x, [0, 1]));
  return div(sub(e_ix, rec(e_ix)), [0, 2]);
}

function cos(x) {
  return sin(add(x, [pi / 2, 0]));
}

function tan(x) {
  return div(sin(x), cos(x));
}

function sinh(x) {
  var e_x = exp(x);
  return scale(sub(e_x, rec(e_x)), 1 / 2);
}

function sinh(x) {
  var e_x = exp(x);
  return scale(sub(e_x, rec(e_x)), 1 / 2);
}

function cosh(x) {
  var e_x = exp(x);
  return scale(add(e_x, rec(e_x)), 1 / 2);
}

function tanh(x) {
  return div(sinh(x), cosh(x));
}

/*

Product function.

*/
function prod(x) {
  return mul(prod_p1(x), prod_p2(x));
}

function prod_p1(x) {
  var p = [1, 0];
  var x_i = x;
  for (var i = 1; i < 200; ++i) {
    if (norm_sqr(x_i) >= 1000000000) {
      break;
    }
    p = mul(p, div(x_i, sub(x_i, [1, 0])));
    x_i = mul(x_i, x);
  }
  return p;
}

function prod_p2(x) {
  var p = [1, 0];
  var x_i = x;
  for (var i = 1; i < 200; ++i) {
    if (norm_sqr(x_i) >= 1000000000) {
      break;
    }
    p = mul(p, div(x_i, add(x_i, [1, 0])));
    x_i = mul(x_i, x);
  }
  return p;
}


/*

Other functions

*/

//Pure guassian function, with guass(x + iy)
function guass(x, a, b, c) {
    if (!a) a = [1, 0]; 
    if (!b) b = [0, 0]; 
    if (!c) c = [1, 0]; 
    return mul(exp(div(scale(sqr(sub(x, b)), -1), scale(sqr(c), 2))), a);
}

//Guassian function in 2 variables, not imaginary. all zero indicies are grouped together
function guass_re(V, a, b, c) {
    if (!a) a = [1, 1]; 
    if (!b) b = [0, 0]; 
    if (!c) c = [1, 1]; 
    var end_sc = a[0] * a[1];
    var exp_0 = scale(sqr(sub(V[0], [b[0], 0])), -1.0 / (2.0 * c[0] * c[0]));
    var exp_1 = scale(sqr(sub(V[1], [b[1], 0])), -1.0 / (2.0 * c[1] * c[1]));
    var exp = exp(add(exp_0, exp_1))
    return scale(exp, end_sc);
}