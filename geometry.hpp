#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <bits/stdc++.h>
using namespace std;
#include "vector.hpp"

template <typename T>
struct Ray
{
    Vector3D<T> p;
    Vector3D<T> v;
    Ray(Vector3D<T> p, Vector3D<T> v) : p(p), v(v) {}
};

template <typename T>
struct Sphere
{
    T r;
    Vector3D<T> c;
    Sphere(T r, Vector3D<T> c) : r(r), c(c) {}
    Vector3D<T> normal(Vector3D<T> p) const
    {
        return normalize(p - c);
    }
};

template <typename T>
struct Cylinder
{
    T r, w;
    Vector3D<T> c, v;
    Cylinder(T r, T w, Vector3D<T> c, Vector3D<T> v) : r(r), w(w), c(c), v(v) {}
    Vector3D<T> normal(Vector3D<T> p) const
    {
        Vector3D<T> n = p - c;
        n = n - (n * v) * v;
        return normalize(n);
    }
};

template <typename T>
struct SphereSector
{
    T r, theta;
    Vector3D<T> c, v;
    SphereSector(T r, T theta, Vector3D<T> c, Vector3D<T> v) : r(r), theta(theta), c(c), v(v) {}
    Vector3D<T> normal(Vector3D<T> p) const
    {
        return normalize(p - c);
    }
};

template <typename T>
struct Plane
{
    Vector3D<T> n, p;
    Plane(Vector3D<T> p, Vector3D<T> n) : n(n), p(p) {}
    Vector3D<T> normal(Vector3D<T> p) const
    {
        return n;
    }
};

template <typename T>
struct Drop
{
    Cylinder<T> cyl;
    SphereSector<T> l, r;
    Drop(Cylinder<T> cyl, SphereSector<T> l, SphereSector<T> r) : cyl(cyl), l(l), r(r) {}
    Vector3D<T> normal(Vector3D<T> p) const
    {
        if (abs((p - cyl.c) * cyl.v) <= cyl.w / 2)
            return cyl.normal(p);
        if ((p - cyl.c) * cyl.v > 0)
            return r.normal(p);
        return l.normal(p);
    }
};

// this is just to make the priority queue work
template <typename T>
bool operator<(const Ray<T> &r1, const Ray<T> &r2)
{
    return norm(r1.p) < norm(r2.p);
}

// CHECK IF POINT IS INSIDE DETECTOR SHAPE (CIRCLE)
double in_shape(Vector3D<double> const &p, Vector3D<double> const &det_corn, double const DETECTOR_SIZE)
{
    return hypot(det_corn.y + DETECTOR_SIZE / 2 - p.y, det_corn.z + DETECTOR_SIZE / 2 - p.z) <= DETECTOR_SIZE / 2;
}

#endif