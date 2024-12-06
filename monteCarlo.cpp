#pragma once
#include "gaussianQuadratures.h"
#include "sampleFunctions.h"

typedef double Doub;
template <class T> 
struct NRf3 
{
    Doub xsav, ysav; 
    T* func3d; // Must be a pointer
    Doub operator()(const Doub z) // The integrand f(x, y, z) evaluated at fixed x and y
    {
        return (*func3d)(xsav, ysav, z); // Must dereference explicitly
    }
};
template <class T, class Z1, class Z2>
struct NRf2 {
    NRf3<T> f3; 

    Z1& z1;
    Z2& z2;
    NRf2(Z1& zz1, Z2& zz2) : z1(zz1), z2(zz2) {}

    Doub operator()(const Doub y) // This is G of eq. (4.8.4)
    {
        f3.ysav = y;
        return qgaus(f3, z1(f3.xsav, y), z2(f3.xsav, y));
    }
};


template <class T, class Y1, class Y2, class Z1, class Z2>
struct NRf1 {
    Y1& y1;
    Y2& y2;

    NRf2<T, Z1, Z2> f2;
 
    NRf1(Y1& yy1, Y2& yy2, Z1& z1, Z2& z2) : y1(yy1), y2(yy2), f2(z1, z2) {}

    Doub operator()(const Doub x) // This is H of eq. (4.8.5).
    {
        f2.f3.xsav = x;
        return qgaus(f2, y1(x), y2(x));
    }
};

template <class T, class Y1, class Y2, class Z1, class Z2>
Doub quad3d(T& func, const Doub x1, const Doub x2, Y1& y1, Y2& y2, Z1& z1, Z2& z2)
// Returns the integral of a user-supplied function func over a three-dimensional region specified
// by the limits x1, x2, and by the user-supplied functions y1, y2, z1, and z2, as defined in (4.8.2).
// Integration is performed by calling qgaus recursively.
{
    NRf1<T, Y1, Y2, Z1, Z2> f1(y1, y2, z1, z2);
    f1.f2.f3.func3d = &func; // Must set as a pointer
    return qgaus(f1, x1, x2);
}


