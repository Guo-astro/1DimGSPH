#pragma once

#include<iostream>
#include<iomanip>
template<class T>
class Vector1 {
public:
	T x;
	Vector1() :
			x(T(0)) {
	}
	Vector1(const T _x) :
			x(_x) {
	}
	Vector1(const Vector1 & src) :
			x(src.x) {
	}

	typedef T DataType;
	static const int DIM = 1;

	const Vector1 & operator =(const Vector1 & rhs) {
		x = rhs.x;
		return (*this);
	}

	const Vector1 & operator =(const T s) {
		x = s;
		return (*this);
	}

	Vector1 operator +(const Vector1 & rhs) const {
		return Vector1(x + rhs.x);
	}
	const Vector1 & operator +=(const Vector1 & rhs) {
		(*this) = (*this) + rhs;
		return (*this);
	}
	Vector1 operator -(const Vector1 & rhs) const {
		return Vector1(x - rhs.x);
	}

	const Vector1 & operator -=(const Vector1 & rhs) {
		(*this) = (*this) - rhs;
		return (*this);
	}

	const Vector1 & operator *=(const T s) {
		(*this) = (*this) * s;
		return (*this);
	}
	Vector1 operator *(const T s) const {
		return Vector1(x * s);
	}
	friend Vector1 operator *(const T s, const Vector1 & v) {
		return (v * s);
	}
	Vector1 operator /(const T s) const {
		return Vector1(x / s);
	}

	Vector1 inv() const {
		return Vector1(1.0 / x);
	}
	const Vector1 & operator /=(const T s) {
		(*this) = (*this) / s;
		return (*this);
	}

	const Vector1 & operator +() const {
		return (*this);
	}

	const Vector1 operator -() const {
		return Vector1(-x);
	}

	T operator *(const Vector1 & rhs) const {
		return (x * rhs.x);
	}

	Vector1 operator ^(const Vector1 & rhs) const {
		return Vector1(0);
	}
	template<typename U>
	operator Vector1<U>() const {
		return Vector1<U>(static_cast<U>(x));
	}

	T getMax() const {
		return x;
	}

	T getMin() const {
		return x;
	}

	template<class F>
	Vector1 applyEach(F f) const {
		return Vector1(f(x));
	}

	template<class F>
	friend Vector1 ApplyEach(F f, const Vector1 & arg1, const Vector1 & arg2) {
		return Vector1(f(arg1.x, arg2.x));
	}

	friend std::ostream & operator <<(std::ostream & c, const Vector1 & u) {
		c << u.x;
		return c;
	}

	friend std::istream & operator >>(std::istream & c, Vector1 & u) {
		c >> u.x;
		return c;
	}

	const T & operator[](const int i) const {
		return (&x)[i];
	}

	T & operator[](const int i) {
		return (&x)[i];
	}

	T getDistanceSQ(const Vector1 & u) const {
		T dx = x - u.x;
		return dx * dx;
	}
	bool operator ==(const Vector1 & u) const {
		return ((x == u.x));
	}
	bool operator !=(const Vector1 & u) const {
		return ((x != u.x));
	}
};

template<>
inline Vector1<float> Vector1<float>::operator /(const float s) const {
	const float inv_s = 1.0f / s;
	return Vector1(x * inv_s);
}
template<>
inline Vector1<double> Vector1<double>::operator /(const double s) const {
	const double inv_s = 1.0 / s;
	return Vector1(x * inv_s);
}
