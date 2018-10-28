/*
 * polyEquationsSolver.h
 *
 *  Created on: Oct 28, 2018
 *      Author: guo
 */

#ifndef POLYEQUATIONSSOLVER_H_
#define POLYEQUATIONSSOLVER_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct {
  double r, i;
} Cdouble;

static const Cdouble _C0 = { 0.0, 0.0 };

static inline Cdouble
_Cdouble(double r, double i)
{
  Cdouble res;

  res.r = r;
  res.i = i;
  return res;
}

static inline Cdouble
_Cdouble1(double r)
{
  return _Cdouble(r, 0.0);
}

static inline int
cis0(Cdouble x)
{
  return x.r == 0.0 && x.i == 0.0;
}

static inline Cdouble
cneg(Cdouble x)
{
  x.r = -x.r;
  x.i = -x.i;
  return x;
}

static inline Cdouble
cadd(Cdouble a, Cdouble b)
{
  a.r += b.r;
  a.i += b.i;
  return a;
}

static inline Cdouble
csub(Cdouble a, Cdouble b)
{
  a.r -= b.r;
  a.i -= b.i;
  return a;
}

static inline Cdouble
cmul(Cdouble a, Cdouble b)
{
  double ar = a.r;
  double ai = a.i;
  double br = b.r;
  double bi = b.i;

  a.r = ar * br - ai * bi;
  a.i = ar * bi + ai * br;
  return a;
}

static inline Cdouble
cdiv(Cdouble a, Cdouble b)
{
  double ar = a.r;
  double ai = a.i;
  double br = b.r;
  double bi = b.i;
  double bb = br * br + bi * bi;

  if (bb == 0.0) {
    fprintf(stderr, "cdiv((%e %e), (%e, %e)) Error: divided by 0.\n",
	    ar, ai, br, bi);
    exit(1);
  }
  a.r = (ar * br + ai * bi) / bb;
  a.i = (ai * br - ar * bi) / bb;
  return a;
}

static inline double
Cabs(Cdouble a)
{
  double r = a.r;
  double i = a.i;
  return sqrt(r * r + i * i);
}

Cdouble
Csqrt(Cdouble x)
{
  double r = x.r;
  double i = x.i;
  double len;

  if (i == 0.0) {
    if (r >= 0.0)
      return _Cdouble1(sqrt(r));
    else
      return _Cdouble(0.0, sqrt(-r));
  }
  len = sqrt(r * r + i * i);
  x.r = sqrt((len + r) * 0.5);
  i = sqrt((len - r) * 0.5);
  if (i < 0.0) i = -i;
  x.i = i;
  return x;
}

Cdouble
ccbrt(Cdouble x)
{
  double r = x.r;
  double i = x.i;
  double clen, ctheta;

  if (i == 0.0)
    return _Cdouble1(cbrt(r));

  clen = cbrt(sqrt(r * r + i * i));
  ctheta = atan2(i, r) / 3.0;
  x.r = clen * cos(ctheta);
  x.i = clen * sin(ctheta);
  return x;
}

typedef struct {
  int n;
  Cdouble ans[4];
} Cans;

Cans
equation1(double a, double b)
{
  Cans ans;

  ans.ans[0] = ans.ans[1] = ans.ans[2] = ans.ans[3] = _C0;
  if (a == 0.0) {
    if (b == 0.0)
      ans.n = -1;
    else
      ans.n = 0;
  }
  else {
    ans.n = 1;
    ans.ans[0] = _Cdouble1(-b/a);
  }
  return ans;
}

Cans
cequation1(Cdouble a, Cdouble b)
{
  Cans ans;

  ans.ans[0] = ans.ans[1] = ans.ans[2] = ans.ans[3] = _C0;

  if (cis0(a)) {
    if (cis0(b))
      ans.n = -1;
    else
      ans.n = 0;
  }
  else {
    ans.n = 1;
    ans.ans[0] = cdiv(cneg(b), a);
  }
  return ans;
}

Cans
equation2(double a, double b, double c)
{
  Cans ans;

  ans.ans[0] = ans.ans[1] = ans.ans[2] = ans.ans[3] = _C0;

  if (a == 0.0)
    return equation1(b, c);

  b = -0.5 * b / a;
  c /= a;
  a = b * b - c;
  if (a >= 0.0) {
    a = sqrt(a);
    ans.ans[0] = _Cdouble1(b + a);
    ans.ans[1] = _Cdouble1(b - a);
  }
  else {
    a = sqrt(-a);
    ans.ans[0] = _Cdouble(b, a);
    ans.ans[1] = _Cdouble(b, -a);
  }
  ans.n = 2;
  return ans;
}

Cans
cequation2(Cdouble a, Cdouble b, Cdouble c)
{
  Cans ans;

  ans.ans[0] = ans.ans[1] = ans.ans[2] = ans.ans[3] = _C0;

  if (cis0(a))
    return cequation1(b, c);

  b = cdiv(cmul(_Cdouble1(-0.5), b), a);
  c = cdiv(c, a);
  a = Csqrt(csub(cmul(b, b), c));
  ans.ans[0] = cadd(b, a);
  ans.ans[1] = csub(b, a);
  ans.n = 2;
  return ans;
}

Cans
equation3(double a, double b, double c, double d)
{
  Cans ans;
  Cdouble w, ww, cA3;
  double A, B, C, A3;
  double p, q, p3, q2;

  if (a == 0.0)
    return equation2(b, c, d);
  A = b / a;
  B = c / a;
  C = d / a;

  w = _Cdouble(-0.5, 0.5 * sqrt(3.0));
  ww = _Cdouble(-0.5, -0.5 * sqrt(3.0));
  A3 = A / 3.0;
  cA3 = _Cdouble1(A3);
  p = B - A * A3;
  q = 2.0 * A3 * A3 * A3 - A3 * B + C;
  p3 = p / 3.0;
  q2 = -0.5 * q;

  d = q2 * q2 + p3 * p3 * p3;
  if (d >= 0.0) {
    d = sqrt(d);
    double m = cbrt(q2 + d);
    double n = cbrt(q2 - d);
    ans.ans[0] = _Cdouble1(m + n - A3);
    ans.ans[1] = csub(cadd(cmul(_Cdouble1(m), w), cmul(_Cdouble1(n), ww)), cA3);
    ans.ans[2] = csub(cadd(cmul(_Cdouble1(m), ww), cmul(_Cdouble1(n), w)), cA3);
  }
  else {
    d = sqrt(-d);
    Cdouble cm = ccbrt(_Cdouble(q2, d));
    Cdouble cn = ccbrt(_Cdouble(q2, -d));
    ans.ans[0] = csub(cadd(cm, cn), cA3);
    ans.ans[1] = csub(cadd(cmul(cm, w), cmul(cn, ww)), cA3);
    ans.ans[2] = csub(cadd(cmul(cm, ww), cmul(cn, w)), cA3);
  }
  ans.n = 3;
  return ans;
}

static double
Find1RealAns(Cans res)
{
  int n = res.n;
  double ii;
  int id, i;

  if (n <= 0) return 0.0;
  ii = fabs(res.ans[0].i);
  id = 0;
  for (i = 1; i < n; ++i) {
    double ti = fabs(res.ans[i].i);
    if (ti < ii) {
      ii = ti;
      id = i;
    }
  }
  return res.ans[id].r;
}

Cans
equation4f(double A, double B, double C)
{
  if (B == 0.0) {
    Cans res = equation2(1.0, A, C);
    Cdouble x1 = Csqrt(res.ans[0]);
    Cdouble x2 = Csqrt(res.ans[1]);
    res.ans[0] = x1;
    res.ans[1] = cneg(x1);
    res.ans[2] = x2;
    res.ans[3] = cneg(x2);
    res.n = 4;
    return res;
  }
  else {
    /* y^4 + Ay^2 + By + C =
     * (y^2+L)^2 + (A-2L)(y^2 + B/(A-2L)y + (C-L^2)/(A-2L) =
     * (y^2+L)^2 + (A-2L)[(y+B/(2(A-2L)))^2+(C-L^2)/(A-2L)-B^2/(4(A-2L)^2)) =
     * (y^2+L)^2 + (A-2L)[(y+B/(2(A-2L)))^2,
     * if (C-L^2)-B^2/(4(A-2L)) = 0.
     * (y^2+L)^2 + M(y + N)^2 = 0
     */
    Cdouble X, b1, c1, b2, c2;
    double L, M, N;
    Cans res1, res2;

    res1 = equation3(1.0, -0.5 * A, -C, 0.5 *  A * C - 0.125 * B * B);
    L = Find1RealAns(res1);
    M = A - L * 2.0;
    N = B / (2.0 * M);
    X = Csqrt(_Cdouble1(-M));
    b1 = cneg(X);
    c1 = csub(_Cdouble1(L), cmul(X, _Cdouble1(N)));
    b2 = X;
    c2 = cadd(_Cdouble1(L), cmul(X, _Cdouble1(N)));

    res1 = cequation2(_Cdouble(1.0, 0.0), b1, c1);
    res2 = cequation2(_Cdouble(1.0, 0.0), b2, c2);
    res1.ans[2] = res2.ans[0];
    res1.ans[3] = res2.ans[1];
    res1.n = 4;
    return res1;
  }
}

Cans
equation4sub(double a, double b, double c, double d)
{
  double A, B, C;
  double aa = a * a;
  Cdouble a4 = _Cdouble1(a / 4.0);
  Cans res;

  /* x^4 + ax^3 + bx^2 + cx + d = y^4 + Ay^2 + By + C,
   *  where x = y - a/4
   */

  A = b - (3.0 / 8.0) * aa;
  B = c + a * aa / 8.0 - a * b / 2.0;
  C = d - (3.0 / 256.0) * aa * aa + aa * b / 16.0 - a * c / 4.0;

  res = equation4f(A, B, C);
  res.ans[0] = csub(res.ans[0], a4);
  res.ans[1] = csub(res.ans[1], a4);
  res.ans[2] = csub(res.ans[2], a4);
  res.ans[3] = csub(res.ans[3], a4);
  return res;
}

Cans
equation4(double a, double b, double c, double d, double e)
{
  if (a == 0.0)
    return equation3(b, c, d, e);

  return equation4sub(b / a, c / a, d / a, e / a);
}



Cdouble
poly4(double a, double b, double c, double d, double e, Cdouble x)
{
  Cdouble xx = cmul(x, x);
  Cdouble xxx = cmul(xx, x);
  Cdouble xxxx = cmul(xx, xx);
  Cdouble ar = cmul(_Cdouble1(a), xxxx);
  Cdouble br = cmul(_Cdouble1(b), xxx);
  Cdouble cr = cmul(_Cdouble1(c), xx);
  Cdouble dr = cmul(_Cdouble1(d), x);
  Cdouble er = _Cdouble1(e);
  return cadd(cadd(cadd(ar, br), cadd(cr, dr)), er);
}

double getadouble(void)
{
  double x;

  while (scanf("%lf", &x) != 1) {
    fprintf(stderr, "Please type a number.\n");
  }
  return x;
}

#endif /* POLYEQUATIONSSOLVER_H_ */
