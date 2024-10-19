#include <bits/stdc++.h>
using namespace std;

// Library to calculate intersections between geometries
#include "intersections.hpp"
// Simple 3D vector class
#include "vector.hpp"
// Library to implement geometries
#include "geometry.hpp"

// Computes reflectivity of surface interface according to Fresnel equations
double reflectivity(double n1, double n2, double theta1)
{
    double n1n2sintheta = n1 / n2 * sin(theta1);
    double sqrtn1n2sintheta = sqrt(1 - n1n2sintheta * n1n2sintheta);

    double R_s = (n1 * cos(theta1) - n2 * sqrtn1n2sintheta) / (n1 * cos(theta1) + n2 * sqrtn1n2sintheta);
    double R_p = (n1 * sqrtn1n2sintheta - n2 * cos(theta1)) / (n1 * sqrtn1n2sintheta + n2 * cos(theta1));
    double R = (R_s * R_s + R_p * R_p) / 2;
    return R;
}

