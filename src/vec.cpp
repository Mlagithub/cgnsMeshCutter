template <class Type>
MeshCut::Vector<Type>::Vector()
{
    value_[0] = value_[1] = value_[2] = 0;
}

template <class Type>
MeshCut::Vector<Type>::Vector(const Type &initValue)
{
    value_[0] = value_[1] = value_[2] = initValue;
}

template <class Type>
MeshCut::Vector<Type>::Vector(const Vector<Type> &v)
{
    this->value_[0] = v.x();
    this->value_[1] = v.y();
    this->value_[2] = v.z();
}

template <class Type>
MeshCut::Vector<Type>::Vector(const Type &vx, const Type &vy, const Type &vz)
{
    value_[0] = vx;
    value_[1] = vy;
    value_[2] = vz;
}