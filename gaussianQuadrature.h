#pragma once
// DISCLAIMER: This is just an example of Gaussian quadrature routine that contains the tabulated
// abscissas and weights for the case W(x) = 1 and N = 10
typedef double Doub;
template <class T>
Doub qgaus(T& func, const Doub a, const Doub b)
// Returns the integral of the function or functor func between a and b, by ten-point Gauss
// Legendre integration: the function is evaluated exactly ten times at interior points in the range
// of integration.
{
    //Here are the abscissas and weights:
    static const Doub x[] = { 0.1488743389816312, 0.4333953941292472,
    0.6794095682990244, 0.8650633666889845, 0.9739065285171717 };

    static const Doub w[] = { 0.2955242247147529, 0.2692667193099963,
    0.2190863625159821, 0.1494513491505806, 0.0666713443086881 };

    Doub xm = 0.5 * (b + a);
    Doub xr = 0.5 * (b - a);
    Doub s = 0;

    for (int j = 0; j < 5; j++) {
        Doub dx = xr * x[j];
        // Will be twice the average value of the function, since the
        // ten weights (five numbers above each used twice) sum to 2.
        s += w[j] * (func(xm + dx) + func(xm - dx));
    }

    return s *= xr; // Scale the answer to the range of integration
}


