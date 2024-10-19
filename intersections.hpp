#ifndef INTERSECTIONS_HPP
#define INTERSECTIONS_HPP

#include <bits/stdc++.h>
using namespace std;
#include "geometry.hpp"

template <typename T>
bool has_solution(T a, T b, T c)
{
    return b * b - 4 * a * c >= 0;
}

template <typename T>
bool is_tangent(T a, T b, T c)
{
    return b * b - 4 * a * c == 0;
}

template <typename T>
pair<T, T> solve(T a, T b, T c)
{
    T delta = b * b - 4 * a * c;
    if (delta < 0)
        throw runtime_error("imaginary roots");

    T t1 = (-b - sqrt(delta)) / (2 * a);
    T t2 = (-b + sqrt(delta)) / (2 * a);
    if (t1 < t2)
        return {t1, t2};
    return {t2, t1};
}

template <typename T>
T distance(const Ray<T> &r, const Vector3D<T> &p)
{
    return ((p - r.p) ^ r.v) / norm(r.v);
}

template <typename T>
bool operator<(const Ray<T> &r, const Sphere<T> &s)
{
    Vector3D<T> oc = r.p - s.c;
    return has_solution(r.v * r.v, 2 * (oc * r.v), oc * oc - s.r * s.r);
}

template <typename T>
bool operator<(const Ray<T> &r, const Cylinder<T> &cyl)
{
    Vector3D<T> oc = r.p - cyl.c;
    Vector3D<T> v = normalize(r.v);
    Vector3D<T> d = normalize(cyl.v);
    Vector3D<T> d_parallel = d * (v * d);
    Vector3D<T> d_perpendicular = v - d_parallel;

    Vector3D<T> dp_parallel = d * (oc * d);
    Vector3D<T> dp_perpendicular = oc - dp_parallel;

    T a = d_perpendicular * d_perpendicular;
    T b = 2. * d_perpendicular * dp_perpendicular;
    T c = dp_perpendicular * dp_perpendicular - cyl.r * cyl.r;

    if (!has_solution(a, b, c))
        return false;
    pair<T, T> t = solve(a, b, c);

    Vector3D<T> p1 = r.p + v * t.first;
    Vector3D<T> p2 = r.p + v * t.second;

    if (max(t.first, t.second) <= 0)
        return false;
    if (t.first > 0 && abs((p1 - cyl.c) * cyl.v) < cyl.w / 2)
        return true;
    return abs((p2 - cyl.c) * cyl.v) < cyl.w / 2;
}

template <typename T>
bool operator<(const Ray<T> &r, const SphereSector<T> &s)
{
    Vector3D<T> oc = r.p - s.c;
    Vector3D<T> v = normalize(r.v);
    T a = v * v;
    T b = 2.0 * oc * v;
    T c = oc * oc - s.r * s.r;
    if (!has_solution(a, b, c))
        return false;
    pair<T, T> t = solve(a, b, c);
    Vector3D<T> p1 = r.p + r.v * t.first;
    Vector3D<T> p2 = r.p + r.v * t.second;
    if (max(t.first, t.second) <= 0)
        return false;
    if ((normalize(p2 - s.c) * s.v) >= cos(s.theta))
        return true;
    return (normalize(p1 - s.c) * s.v) >= cos(s.theta);
}

template <typename T>
bool operator<(const Ray<T> &r, const Plane<T> &p)
{
    return r.v * (p.p - r.p) > 0;
}

template <typename T>
bool operator<(const Ray<T> &r, const Drop<T> &d)
{
    return (r < d.cyl || r < d.l || r < d.r);
}

template <typename T>
Vector3D<T> operator&(const Ray<T> &r, const Sphere<T> &s)
{
    Vector3D<T> oc = r.p - s.c;
    T a = r.v * r.v;
    T b = 2 * (oc * r.v);
    T c = oc * oc - s.r * s.r;
    T delta = b * b - 4 * a * c;
    T t1 = (-b - sqrt(delta)) / (2 * a);
    T t2 = (-b + sqrt(delta)) / (2 * a);
    Vector3D<T> p1 = r.p + r.v * t1;
    Vector3D<T> p2 = r.p + r.v * t2;
    if (t1 < t2)
        return p1;
    return p2;
}

template <typename T>
Vector3D<T> operator&(const Ray<T> &r, const Cylinder<T> &cyl)
{
    Vector3D<T> oc = r.p - cyl.c;
    Vector3D<T> v = normalize(r.v);
    Vector3D<T> d = normalize(cyl.v);
    Vector3D<T> d_parallel = d * (v * d);
    Vector3D<T> d_perpendicular = v - d_parallel;

    Vector3D<T> dp_parallel = d * (oc * d);
    Vector3D<T> dp_perpendicular = oc - dp_parallel;

    T a = d_perpendicular * d_perpendicular;
    T b = 2. * d_perpendicular * dp_perpendicular;
    T c = dp_perpendicular * dp_perpendicular - cyl.r * cyl.r;

    if (!has_solution(a, b, c))
        return Vector3D<T>(0, 0, 0);
    pair<T, T> t = solve(a, b, c);

    Vector3D<T> p1 = r.p + v * t.first;
    Vector3D<T> p2 = r.p + v * t.second;
    if (max(t.first, t.second) <= 0)
        return Vector3D<T>(0, 0, 0);
    if (t.first > 0 && abs((p1 - cyl.c) * cyl.v) < cyl.w / 2)
        return p1;
    return p2;
}

template <typename T>
Vector3D<T> operator&(const Ray<T> &r, const SphereSector<T> &s)
{
    Vector3D<T> oc = r.p - s.c;
    Vector3D<T> v = normalize(r.v);
    Vector3D<T> d = normalize(s.v);

    T a = v * v - (v * d) * (v * d);
    T b = 2 * ((oc * v) - (oc * d) * (v * d));
    T c = (oc * oc) - (oc * d) * (oc * d) - s.r * s.r;
    pair<T, T> t = solve(a, b, c);
    Vector3D<T> p1 = r.p + r.v * t.first;
    Vector3D<T> p2 = r.p + r.v * t.second;
    if (max(t.first, t.second) <= 0)
        return Vector3D<T>(0, 0, 0);
    if ((normalize(p1 - s.c) * s.v) >= cos(s.theta) && t.first > 0)
        return p1;
    return p2;
}

template <typename T>
Vector3D<T> operator&(const Ray<T> &r, const Plane<T> &p)
{
    Vector3D<T> n = p.n;
    Vector3D<T> p0 = p.p;
    Vector3D<T> v = normalize(r.v);
    Vector3D<T> l0 = r.p;
    T d = (n * (p0 - l0)) / (n * v);
    return l0 + v * d;
}

template <typename T>
Vector3D<T> operator&(const Ray<T> &r, const Drop<T> &d)
{
    Vector3D<T> p = Vector3D<T>(0, 0, 1) * d.cyl.w * 200.;
    if (r < d.cyl)
    {
        auto tmp = r & d.cyl;
        if (norm(tmp - r.p) < norm(p - r.p))
            p = tmp;
    }
    if (r < d.l)
    {
        auto tmp = r & d.l;
        if (norm(tmp - r.p) < norm(p - r.p))
            p = tmp;
    }
    if (r < d.r)
    {
        auto tmp = r & d.r;
        if (norm(tmp - r.p) < norm(p - r.p))
            p = tmp;
    }
    return p;
}

#endif