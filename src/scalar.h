#ifndef SCALAR_H
#define SCALAR_H


namespace MeshCut
{

typedef unsigned char ComponentIndex;

#ifdef SINGLE_PRECISE
typedef float Scalar;
#else
typedef double Scalar;
#endif

template <class Type>
class Traits : public Type
{
public:
    explicit Traits(const Type &t) : Type(t)
    {
    }

    static constexpr unsigned char nComponents() { return Type::nComponents(); }
};


/**
 * @brief Traits<Type> 的模板特化
 * 
 * @tparam Type = Scalar
 */
template <>
class Traits<Scalar>
{
    Scalar p_;
    static constexpr ComponentIndex nComponents_ = 1; // - Number of components

public:
    static constexpr ComponentIndex nComponents() { return nComponents_; }

    explicit Traits(const Scalar &val) : p_(val)
    {
    }

    operator Scalar() const
    {
        return p_;
    }

    operator Scalar &()
    {
        return p_;
    }

    static Scalar &component(Scalar &s, const ComponentIndex)
    {
        return s;
    }
};

} // namespace MeshCut

#endif
