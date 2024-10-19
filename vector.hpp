#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <bits/stdc++.h>
using namespace std;

template <typename T>
struct Vector3D
{
    T x, y, z;
    Vector3D<T>(T x, T y, T z) : x(x), y(y), z(z) {}
    Vector3D<T>() : x(0), y(0), z(0) {}
};

template <typename T>
Vector3D<T> operator*(const Vector3D<T> &v, T scalar)
{
    return Vector3D<T>(v.x * scalar, v.y * scalar, v.z * scalar);
}

template <typename T>
Vector3D<T> operator*(T scalar, const Vector3D<T> &v)
{
    return Vector3D<T>(v.x * scalar, v.y * scalar, v.z * scalar);
}

template <typename T>
Vector3D<T> operator/(const Vector3D<T> &v, T scalar)
{
    return Vector3D<T>(v.x / scalar, v.y / scalar, v.z / scalar);
}

template <typename T>
Vector3D<T> operator+(const Vector3D<T> &v1, const Vector3D<T> &v2)
{
    return Vector3D<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

template <typename T>
Vector3D<T> operator-(const Vector3D<T> &v1, const Vector3D<T> &v2)
{
    return Vector3D<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

// cross product
template <typename T>
Vector3D<T> operator^(const Vector3D<T> &v1, const Vector3D<T> &v2)
{
    return Vector3D<T>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// dot product
template <typename T>
T operator*(const Vector3D<T> &v1, const Vector3D<T> &v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
T norm(const Vector3D<T> &v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

template <typename T>
Vector3D<T> normalize(const Vector3D<T> &v)
{
    T n = norm(v);
    return Vector3D<T>(v.x / n, v.y / n, v.z / n);
}

template <typename T>
Vector3D<T> reflect(const Vector3D<T> &v, const Vector3D<T> &normal)
{
    return v - 2 * (v * normal) * normal;
}

template <typename T>
Vector3D<T> refract(const Vector3D<T> &v, const Vector3D<T> &normal, T n1, T n2)
{
    Vector3D<T> nrm = normal * v < 0 ? normal : normal * -1.;
    T ior_ratio = n1 / n2;
    T cos_theta = -(v * nrm) / norm(v);
    T sin_theta = sqrt(1 - cos_theta * cos_theta);
    T sin_phi = ior_ratio * sin_theta;
    T cos_phi = sqrt(1 - sin_phi * sin_phi);
    return ior_ratio * v + (ior_ratio * cos_theta - cos_phi) * nrm * norm(v);
}

template <typename T>
ostream &operator<<(ostream &os, const Vector3D<T> &v)
{
    os << "(" << v.x << ",\t" << v.y << ",\t" << v.z << ")";
    return os;
}

#endif