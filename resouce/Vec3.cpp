//
// Created by Tanakit Sachati on 2019-07-11.
//

#include "../header/Vec3.h"
#include "../header/Resources.h"
#include "../header/SphereNode.h"


#include <math.h>


namespace pd
{

// ============================================================================
// ==== De-/Constructor(s) ====================================================
// ============================================================================

/**
 * @brief Default constructor, initialize all components with zero.
 */
    template <class T>
    Vec3<T>::Vec3()
            : x(0), y(0), z(0)
    {
    }

/**
 * @brief Initialize each component with the respective parameter value.
 * @param x The first component, x.
 * @param y The second component, y.
 * @param z The third component, z.
 */
    template <class T>
    Vec3<T>::Vec3( T x, T y, T z )
            : x(x), y(y), z(z)
    {
    }

/**
 * @brief Initialize each component with the respective parameter value.
 * @param values Pointer to the three vector components (x,y and z).
 */
    template <class T>
    Vec3<T>::Vec3( T values[3] )
            : x(values[0]), y(values[1]), z(values[2])
    {
    }

// ============================================================================
// ==== Static ================================================================
// ============================================================================

/**
 * @brief Find intersections of a ray with a sphere.
 * @param rayOrigin 3D point from which the ray emits.
 * @param rayDir Normalized 3D vector that is the direction the ray is pointing.
 * @param sphere The 3D sphere that is to be intersected with the ray.
 * @param t1 The resulting distance of the first ray intersection.
 * @param t2 The resulting distance of the second intersection.
 *
 * @see https://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm
 *      on a more in depth explanation of this procedure.
 */
    template <class T>
    bool Vec3<T>::calcIntersection( const Vec3<T> &rayOrigin, const Vec3<T> &rayDir, SphereNode<T> *sphere, T &t1, T &t2 )
    {
        // Catch zero-normals, otherwise division by zero later on
        if (rayDir.isZero())
            return false;

        Vec3<T> w = rayOrigin - sphere->m_center;
        T A = rayDir.dot( rayDir );
        T B = 2 * w.dot(rayDir);
        T C = w.dot(w) - (sphere->m_radius * sphere->m_radius);
        T D = B*B - 4*A*C;

        if (D >= 0)
        {
            t1 = (-sqrtf(D)-B) / (2*A);
            t2 = (+sqrtf(D)-B) / (2*A);
            return true;
        }

        // D<0 => point not in sphere & ray does _NOT_ intersect sphere.
        return false;
    }

/**
 * @brief Find intersections of a ray with a sphere.
 * @param rayOrigin 3D point from which the ray emits.
 * @param rayDir Normalized 3D vector that is the direction the ray is pointing.
 * @param sphere The 3D sphere that is to be intersected with the ray.
 * @param t1 The resulting distance of the first ray intersection.
 *
 * @see https://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm
 *      on a more in depth explanation of this procedure.
 */
    template <class T>
    bool Vec3<T>::calcIntersection( const Vec3<T> &rayOrigin, const Vec3<T> &rayDir, SphereNode<T> *sphere, T &t1 )
    {
        // Catch zero-normals, otherwise division by zero later on
        if (rayDir.isZero())
            return false;

        Vec3<T> w = rayOrigin - sphere->m_center;
        T A = rayDir.dot( rayDir );
        T B = 2 * w.dot(rayDir);
        T C = w.dot(w) - (sphere->m_radius * sphere->m_radius);
        T D = B*B - 4*A*C;

        if (D >= 0)
        {
            t1 = (-sqrtf(D)-B) / (2*A);
            return true;
        }

        // D<0 => point not in sphere & ray does _NOT_ intersect sphere.
        return false;
    }

/**
 * @brief Check if the sphere's center lies behind the plane.
 * @param planePoint A point that lies in the plane.
 * @param planeNormal Normal that defines the same plane.
 * @param sphere The 3D sphere that is to be checked against the plane.
 * @return Whether the sphere center lies behind the plane.
 * @retval true   Iff sphere center lies behind plane.
 * @retval false  Iff sphere center lies in front of plane.
 *
 * @todo If it is desirable to check whether the sphere is completely behind
 *       the plane the signed_distance needs to be >= sphere->m_radius.
 */
    template <class T>
    bool Vec3<T>::isValidPath( const Vec3<T> &planePoint, const Vec3<T> &planeNormal, SphereNode<T> *sphere )
    {
        // Convert to HNF
        T D = -planeNormal.dot( planePoint );
        T signed_distance = planeNormal.dot(sphere->m_center) + D;

        if (signed_distance >= 0)
            return false;
        else
            return true;
    }

// ============================================================================
// ==== Queries ===============================================================
// ============================================================================

/**
 * @brief Return copy that has all components negative of own values.
 */
    template <class T>
            Vec3<T> Vec3<T>::operator -() const
    {
        return Vec3<T>( -x, -y, -z );
    }

/**
 * @brief Check if this vector has a NaN component.
 * @return Whether this vector has a NaN component.
 * @retval true   If one of the components is NaN.
 * @retval false  Iff no component is NaN.
 */
    template <class T>
    bool Vec3<T>::isNaN() const
    {
        return x != x && y != y && z != z;
    }

/**
 * @brief Checks if this vector is the zero vector.
 * @return Whether this vector is the zero vector.
 * @retval true   Iff this vector has zero for all components.
 * @retval false  If any component is not zero.
 */
    template <class T>
    bool Vec3<T>::isZero() const
    {
        return x == 0 && y == 0 && z == 0;
    }

/**
 * @brief Checks if this vector is <i>nearly</i> the zero vector.
 * @return Whether this vector is <i>nearly</i> the zero vector.
 * @retval true   Iff this vector has a value <i>near</i> zero for all components.
 * @retval false  If any component is not <i>near</i> zero.
 *
 * @note \f$x\f$ is <i>nearly</i> zero means \f$-\varepsilon<x<\varepsilon\f$,
 *       with \f$\varepsilon\f$ being a very small value like 1e-10f.
 */
    template <class T>
    bool Vec3<T>::isNearZero() const
    {
        return (-__PD__EPSILON < x && x < __PD__EPSILON) &&
               (-__PD__EPSILON < y && y < __PD__EPSILON) &&
               (-__PD__EPSILON < z && z < __PD__EPSILON);
    }

/**
 * @brief Calculate the squared length of this vector (for performance).
 * @return The squared length of this vector.
 */
    template <class T>
            T Vec3<T>::calcSqLength() const
    {
        return x*x + y*y + z*z;
    }

/**
 * @brief Calculate the length of this vector.
 * @return The length of the vector.
 *
 * @todo Maybe call calcSqLength inside sqrtf()?
 */
    template <class T>
            T Vec3<T>::calcLength() const
    {
        return sqrtf( x*x + y*y + z*z );
    }

/**
 * @brief Return copy that is the normalized version of this vector.
 * @see Vec3<T>::normalize() for details.
 */
    template <class T>
            Vec3<T> Vec3<T>::normalizedCopy() const
    {
        Vec3<T> result( x, y, z );
        result.normalize();
        return result;
    }

// ============================================================================
// ==== Operations ============================================================
// ============================================================================

/**
 * @brief Equality means all components match.
 */
    template <class T>
    bool Vec3<T>::operator ==( const Vec3<T> &operand ) const
    {
        return x == operand.x && y == operand.y && z == operand.z;
    }

/**
 * @brief Inequality means equality is not present.
 */
    template <class T>
    bool Vec3<T>::operator !=( const Vec3<T> &operand ) const
    {
        return !( *this == operand );
    }

/**
 * @brief True if magnitude smaller than operand.
 */
    template <class T>
    bool Vec3<T>::operator <( const Vec3<T> &operand ) const
    {
        return this->calcLength() < operand.calcLength();
    }

/**
 * @brief True if magnitude larger than operand.
 */
    template <class T>
    bool Vec3<T>::operator >( const Vec3<T> &operand ) const
    {
        return this->calcLength() > operand.calcLength();
    }

/**
 * @brief Return copy that has difference of operand's values and own values.
 */
    template <class T>
            Vec3<T> Vec3<T>::operator -( const Vec3<T> &operand ) const
    {
        return Vec3<T>( x-operand.x, y-operand.y, z-operand.z );
    }

/**
 * @brief Return copy that is the sum of operand's values and the own values.
 */
    template <class T>
            Vec3<T> Vec3<T>::operator +( const Vec3<T> &operand ) const
    {
        return Vec3<T>( x+operand.x, y+operand.y, z+operand.z );
    }

    template <class T>
            Vec3<T> Vec3<T>::operator *( const Vec3<T> &operand ) const
    {
        return Vec3<T>( x*operand.x, y*operand.y, z*operand.z );
    }

    template <class T>
            Vec3<T> Vec3<T>::operator /( const Vec3<T> &operand ) const
    {
        T nx=0,ny=0,nz=0;
        if (operand.x!=0) nx = x/operand.x;
        if (operand.y!=0) ny = y/operand.y;
        if (operand.z!=0) nz = z/operand.z;
        return Vec3<T>( nx, ny, nz );
    }

/**
 * @brief TODO
 */
    template <class T>
            Vec3<T> Vec3<T>::operator -( const T &operand ) const
    {
        return Vec3<T>( x-operand, y-operand, z-operand );
    }

/**
 * @brief TODO
 */
    template <class T>
            Vec3<T> Vec3<T>::operator +( const T &operand ) const
    {
        return Vec3<T>( x+operand, y+operand, z+operand );
    }

/**
 * @brief Return copy that is the product of this vector and the operand.
 */
    template <class T>
            Vec3<T> Vec3<T>::operator *( const T &operand ) const
    {
        return Vec3<T>( x*operand, y*operand, z*operand );
    }

/**
 * @brief Return copy that is the quotient of this vector and the operand.
 */
    template <class T>
            Vec3<T> Vec3<T>::operator /( const T &operand ) const
    {
        return Vec3<T>( x/operand, y/operand, z/operand );
    }

/**
 * @brief Calculates the minimum length from this point to the sphere surface.
 * @param node The 3D sphere to calculate the distance to.
 * @return The minimal distance from the point to the sphere surface.
 */
    template <class T>
            T Vec3<T>::calcMinDistance( SphereNode<T> *node ) const
{
    return (node->m_center-*this).calcLength() - node->m_radius;
}

/**
 * @brief Calculates the maximum length from this point to the sphere surface.
 * @param node The 3D sphere to calculate the distance to.
 * @return The maximal distance from the point to the sphere surface.
 */
template <class T>
        T Vec3<T>::calcMaxDistance( SphereNode<T> *node ) const
{
return (node->m_center-*this).calcLength() + node->m_radius;
}

/**
 * @brief Return dot product of this and the other vector.
 * @param other The vector that will be used to calculate the dot product with.
 * @return The dot product of this and the other vector.
 */
template <class T>
        T Vec3<T>::dot( const Vec3<T> &other ) const
{
    return x*other.x + y*other.y + z*other.z;
}

/**
 * @brief Return the result of the cross product of this and the other vector.
 * @param other The vector that will be used to calculate the cross product.
 * @return The result of the cross product of this and the other vector.
 */
template <class T>
        Vec3<T> Vec3<T>::crossProduct( const Vec3<T> &other ) const
{
    return Vec3<T>( y*other.z - z*other.y,
                    z*other.x - x*other.z,
                    x*other.y - y*other.x );
}

/**
 * @brief Multiply this vector like it's a row and operand like it's a column.
 */
template <class T>
        T Vec3<T>::timesColumn( const Vec3<T> &operand ) const
{
    return x*operand.x + y*operand.y + z*operand.z;
}

///**
// * @brief Multiply this vector like it's a column and operand like it's a row.
// */
//template <class T>
//Matrix3f Vec3<T>::timesRow( const Vec3<T> &operand ) const
//{
//    return Matrix3f( x*operand.x, x*operand.y, x*operand.z,
//                     y*operand.x, y*operand.y, y*operand.z,
//                     z*operand.x, z*operand.y, z*operand.z );
//}

// ============================================================================
// ==== Modification ==========================================================
// ============================================================================

/**
 * @brief Make this vector the difference of this vector and the operand.
 */
template <class T>
        Vec3<T> &Vec3<T>::operator -=( const Vec3<T> &operand )
{
    Vec3<T> result = *this - operand;
    *this = result;
    return *this;
}

/**
 * @brief Make this vector the sum of this vector and the operand.
 */
template <class T>
        Vec3<T> &Vec3<T>::operator +=( const Vec3<T> &operand )
{
    Vec3<T> result = *this + operand;
    *this = result;
    return *this;
}

template <class T>
        Vec3<T> &Vec3<T>::operator *=( const Vec3<T> &operand )
{
    Vec3<T> result = *this * operand;
    *this = result;
    return *this;
}

/**
 * @brief Make this vector the product of this vector and the operand.
 */
template <class T>
        Vec3<T> &Vec3<T>::operator *=( T operand )
{
    Vec3<T> result = *this * operand;
    *this = result;
    return *this;
}

/**
 * @brief Make this vector the quotient of this vector and the operand.
 */
template <class T>
        Vec3<T> &Vec3<T>::operator /=( T operand )
{
    Vec3<T> result = *this / operand;
    *this = result;
    return *this;
}

/**
 * @brief Make the length of this vector unit while keeping the direction.
 */
template <class T>
void Vec3<T>::normalize()
{
    T norm = calcLength();

    // Check to not divide by zero
    if (norm != 0) {
        x /= norm;
        y /= norm;
        z /= norm;
    }
}

// ============================================================================
// ==== Misc. =================================================================
// ============================================================================

/**
 * @brief Print a Vec3<T> by printing each of the components.
 */
template <class T>
std::ostream &operator<<( std::ostream &o, const Vec3<T> &v )
{
    o << "Vec3<T> x=" << v.x << " y=" << v.y << " z=" << v.z;
    return o;
}

/**
 * @brief TODO
 */
template <class T>
        Vec3<T> operator *( const T &lhs, Vec3<T> rhs )
{
    return rhs*lhs;
}

template std::ostream &operator<<( std::ostream &o, const Vec3<float> &v );
template std::ostream &operator<<( std::ostream &o, const Vec3<double> &v );

template Vec3<float> operator *( const float &lhs, Vec3<float> rhs );
template Vec3<double> operator *( const double &lhs, Vec3<double> rhs );

template class Vec3<int>;
template class Vec3<float>;
template class Vec3<double>;

}