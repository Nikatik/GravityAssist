#include "lib.h"

#include <stdio.h>

void RKq (__float128 t, __float128* x0, __float128* y0, __float128* z0,
          __float128* r0, __float128 h, unsigned int s, __float128** k,
          __float128** cab,
          __float128   f (__float128, __float128, __float128, __float128,
                          __float128, __float128),
          __float128   g (__float128, __float128, __float128, __float128,
                          __float128, __float128),
          __float128   q (__float128, __float128, __float128, __float128,
                          __float128, __float128),
          __float128   e (__float128, __float128, __float128, __float128,
                          __float128, __float128),
          __float128 M, bool check)
{
    __float128 tempx;
    __float128 tempy;
    __float128 tempz;
    __float128 tempr;

    for (unsigned int i = 0; i < s; i++)
    {
        tempx = 0;
        tempy = 0;
        tempz = 0;
        tempr = 0;
        for (unsigned int j = 0; j < i; j++)
        {
            tempx += cab[i + 1][j] * k[0][j];
            tempy += cab[i + 1][j] * k[1][j];
            tempz += cab[i + 1][j] * k[2][j];
            tempr += cab[i + 1][j] * k[3][j];
        }
        k[0][i] = f (t + h * cab[0][i], *x0 + h * tempx, *y0 + h * tempy,
                     *z0 + h * tempz, *r0 + h * tempr, M);
        k[1][i] = g (t + h * cab[0][i], *x0 + h * tempx, *y0 + h * tempy,
                     *z0 + h * tempz, *r0 + h * tempr, M);
        k[2][i] = q (t + h * cab[0][i], *x0 + h * tempx, *y0 + h * tempy,
                     *z0 + h * tempz, *r0 + h * tempr, M);
        k[3][i] = e (t + h * cab[0][i], *x0 + h * tempx, *y0 + h * tempy,
                     *z0 + h * tempz, *r0 + h * tempr, M);
    }

    tempx = 0;
    tempy = 0;
    tempz = 0;
    tempr = 0;
    if (check)
    {
        for (unsigned int i = 0; i < s; i++)
        {
            tempx += cab[s + 2][i] * k[0][i];
            tempy += cab[s + 2][i] * k[1][i];
            tempz += cab[s + 2][i] * k[2][i];
            tempr += cab[s + 2][i] * k[3][i];
        }
    }
    else
    {
        for (unsigned int i = 0; i < s; i++)
        {
            tempx += cab[s + 1][i] * k[0][i];
            tempy += cab[s + 1][i] * k[1][i];
            tempz += cab[s + 1][i] * k[2][i];
            tempr += cab[s + 1][i] * k[3][i];
        }
    }

    *x0 += h * tempx;
    *y0 += h * tempy;
    *z0 += h * tempz;
    *r0 += h * tempr;
}

__float128 astep (__float128 T, __float128* x, __float128* y, __float128* z,
                  __float128* r, long long unsigned* i, long long unsigned* j,
                  unsigned int p, unsigned int s, __float128** k,
                  __float128** cab, __float128 tol,
                  __float128 f (__float128, __float128, __float128, __float128,
                                __float128, __float128),
                  __float128 g (__float128, __float128, __float128, __float128,
                                __float128, __float128),
                  __float128 q (__float128, __float128, __float128, __float128,
                                __float128, __float128),
                  __float128 e (__float128, __float128, __float128, __float128,
                                __float128, __float128),
                  __float128 M, bool first)
{
    __float128 tempx, tempy, tempz, tempr, x_, y_, z_, r_, dist, h, temp, fac,
        err;
    x_ = tempx = *x;
    y_ = tempy = *y;
    z_ = tempz = *z;
    r_ = tempr = *r;
    dist       = 0;
    h          = 0.01Q;
    fac        = 1.7Q;
    err        = 0;
    FILE* outf = NULL;
    if (first) outf = fopen ("trajectory.txt", "w");
    else outf = fopen ("trajectoryb.txt", "w");
    if (outf == NULL) { printf ("File for trajectories can`t be opened!\n"); }
    for (*i = *j = 0; T - dist > EPS;)
    {
        if (dist + h > T) h = T - dist;

        RKq (dist, &tempx, &tempy, &tempz, &tempr, h, s, k, cab, f, g, q, e, M,
             false);
        RKq (dist, &x_, &y_, &z_, &r_, h, s, k, cab, f, g, q, e, M, true);

        temp = h;
        h *= fminq (fac, fmaxq (0.7Q, powq (0.98Q * tol /
                                                sqrtq (powq (tempx - x_, 2) +
                                                       powq (tempy - y_, 2) +
                                                       powq (tempz - z_, 2) +
                                                       powq (tempr - r_, 2)),
                                            1.Q / (p + 1))));

        if (h < (__float128) pow (10, -8) && fabsq (T - dist - h) > EPS)
        {
            printf ("\nSomething goes wrong...\n");
            return 0;
        }
        // if (h > 0.01Q) h = 0.01Q;
        if (sqrtq (powq (tempx - x_, 2) + powq (tempy - y_, 2) +
                   powq (tempz - z_, 2) + powq (tempr - r_, 2)) > tol)
        {
            x_ = tempx = *x;
            y_ = tempy = *y;
            z_ = tempz = *z;
            r_ = tempr = *r;
            fac        = 1;
            *j += 1;
            continue;
        }

        err += sqrtq (powq (tempx - x_, 2) + powq (tempy - y_, 2) +
                      powq (tempz - z_, 2) + powq (tempr - r_, 2));
        dist += temp;
        *x  = tempx;
        *y  = tempy;
        *z  = tempz;
        *r  = tempr;
        fac = 1.7Q;
        *i += 1;
        fprintf (
            outf,
            "%16.8Lf  |  %50.35Lf  |  %50.35Lf  |  %50.35Lf  |  %50.35Lf  |  "
            "%10.3Lf\n",
            (long double) dist, (long double) *x, (long double) *y,
            (long double) *z, (long double) *r, (long double) M);
        // printf("%.10Lf\t", (long double)
        // (powq(q(dist,*x,*y,*z,*r,M),2)+powq(e(dist,*x,*y,*z,*r,M),2)));

        if (powq (q (dist, *x, *y, *z, *r, M), 2) +
                    powq (e (dist, *x, *y, *z, *r, M), 2) <
                powq (10, -7) &&
            dist / T > 0.3Q)
            break;
    }
    fclose (outf);
    return err;
}
