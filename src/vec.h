#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <ostream>

#include "scalar.h"

namespace MeshCut
{
    template <class Type>
    class Vector
    {
    private:
        Type value_[3];
        static constexpr ComponentIndex nComponents_ = 3; // - Number of components

    public:
        //Constructors
        // for temporary vector, and initialized to zero
        inline Vector();

        //for given size and initialized to given value
        inline Vector(const Type &initValue);

        //as copy
        inline Vector(const Vector<Type> &v);

        //from (vx,vy,vz) tuple
        inline Vector(const Type &vx, const Type &vy, const Type &vz);

        ~Vector(){};

        static constexpr unsigned char nComponents() { return nComponents_; }

        inline Type &component(const ComponentIndex d)
        {
            return this->value_[d];
        }

        inline Type &operator[](const ComponentIndex i)
        {
            return this->value_[i];
        }

        inline const Type &operator[](const ComponentIndex i) const
        {
            return this->value_[i];
        }

        inline const Type &x() const
        {
            return this->value_[0];
        };
        inline const Type &y() const
        {
            return this->value_[1];
        };
        inline const Type &z() const
        {
            return this->value_[2];
        };

        inline Type &x()
        {
            return this->value_[0];
        };
        inline Type &y()
        {
            return this->value_[1];
        };
        inline Type &z()
        {
            return this->value_[2];
        };

        Vector<Type> &normalize()
        {
            Type norm = this->norm();
            value_[0] /= norm;
            value_[1] /= norm;
            value_[2] /= norm;
            return *this;
        };

        inline Type norm()
        {
            return sqrt(value_[0] * value_[0] + value_[1] * value_[1] + value_[2] * value_[2]);
        }


        Vector<Type> &operator=(const Vector<Type> &right)
        {
            if (this == &right)
                return *this;
            value_[0] = right.value_[0];
            value_[1] = right.value_[1];
            value_[2] = right.value_[2];
            return *this;
        }

        Vector<Type> &operator+=(const Vector<Type> &v)
        {
            value_[0] += v.value_[0];
            value_[1] += v.value_[1];
            value_[2] += v.value_[2];
            return *this;
        }
        Vector<Type> &operator*=(double s)
        {
            value_[0] *= s;
            value_[1] *= s;
            value_[2] *= s;
            return *this;
        }
    };

    //Global Operators
    template <class Type>
    inline Vector<Type> operator+(const Vector<Type> &left, const Vector<Type> &right)
    {
        return Vector<Type>(
            left.x() + right.x(),
            left.y() + right.y(),
            left.z() + right.z());
    }

    template <class Type>
    inline Vector<Type> operator-(const Vector<Type> &left, const Vector<Type> &right)
    {
        return Vector<Type>(
            left.x() - right.x(),
            left.y() - right.y(),
            left.z() - right.z());
    }

    // 向量叉乘
    template <class Type>
    inline Vector<Type> operator*(const Vector<Type> &left, const Vector<Type> &right)
    {
        return Vector<Type>(
            left.y() * right.z() - left.z() * right.y(),
            left.z() * right.x() - left.x() * right.z(),
            left.x() * right.y() - left.y() * right.x());
    }

    // 向量数乘
    template <class Type>
    inline Vector<Type> operator*(double s, const Vector<Type> &v)
    {
        return Vector<Type>(v.x() * s, v.y() * s, v.z() * s);
    }

    // 向量数乘
    template <class Type>
    inline Vector<Type> operator*(const Vector<Type> &v, double s)
    {
        return s * v;
    }

    //the inner product for vectors
    template <class Type>
    inline Type operator&(const Vector<Type> &left, const Vector<Type> &right)
    {
        return Type(
            left.x() * right.x() +
            left.y() * right.y() +
            left.z() * right.z());
    }

    template <class Type>
    inline Vector<Type> operator/(const Vector<Type> &v, double s)
    {
        return Vector<Type>(v.x() / s, v.y() / s, v.z() / s);
    }

    template <typename T>
    std::ostream &operator<<(std::ostream &os, const Vector<T> &rhs)
    {
        os << "[" << rhs.x() << ", " << rhs.y() << ", " << rhs.z() << "]";
        return os;
    }

    typedef Vector<Scalar> ScalarVector; //Type:double shuold be configured in compile-time

    /**
 * @brief Traits<Type> 的模板特化
 * 
 * @tparam Type = Vector
 */
    template <>
    class Traits<ScalarVector>
    {
        static constexpr ComponentIndex nComponents_ = 3;

    public:
        static constexpr ComponentIndex nComponents() { return nComponents_; }

        static Scalar &component(ScalarVector &v, const ComponentIndex d)
        {
            return v.component(d);
        }
    };

} // namespace MeshCut

#include "vec.cpp"

#endif
