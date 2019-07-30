//
// Created by Tanakit Sachati on 2019-07-11.
//

#ifndef UNTITLED4_VEC3_H
#define UNTITLED4_VEC3_H

#include <iostream>

namespace pd
{

    template <class T> class SphereNode;
    template <class T> class Matrix3;

    template <class T>
    class Vec3
    {
    public:
        /*************************************************************************/
        /*************************** public Functions ****************************/
        /*************************************************************************/

        // ---- De-/Constructor(s) ------------------------------------------------

        Vec3();
        Vec3( T x, T y, T z );
        Vec3( T values[3] );

        // ---- Static ------------------------------------------------------------

        static bool calcIntersection( const Vec3<T> &rayOrigin, const Vec3<T> &rayDir, SphereNode<T> *sphere, T &t1, T &t2 );
        static bool calcIntersection( const Vec3<T> &rayOrigin, const Vec3<T> &rayDir, SphereNode<T> *sphere, T &t1 );
        static bool isValidPath( const Vec3<T> &rayOrigin, const Vec3<T> &rayDir, SphereNode<T> *sphere );
        static Vec3<T> ZERO() { return Vec3<T>(0,0,0); };
        static Vec3<T> ONE() { return Vec3<T>(1,1,1); };
        static Vec3<T> UNIT_X() { return Vec3<T>(1,0,0); };
        static Vec3<T> UNIT_Y() { return Vec3<T>(0,1,0); };
        static Vec3<T> UNIT_Z() { return Vec3<T>(0,0,1); };

        // ---- Queries -----------------------------------------------------------

        Vec3<T> operator -()     const;
        bool    isNaN()          const;
        bool    isZero()         const;
        bool    isNearZero()     const;
        T       calcSqLength()   const;
        T       calcLength()     const;
        Vec3<T> normalizedCopy() const;

        // ---- Operations --------------------------------------------------------

        bool  operator ==( const Vec3<T> &operand ) const;
        bool  operator !=( const Vec3<T> &operand ) const;
        bool  operator < ( const Vec3<T> &operand ) const;
        bool  operator > ( const Vec3<T> &operand ) const;
        Vec3<T> operator - ( const Vec3<T> &operand ) const;
        Vec3<T> operator + ( const Vec3<T> &operand ) const;
        Vec3<T> operator * ( const Vec3<T> &operand ) const;
        Vec3<T> operator / ( const Vec3<T> &operand ) const;
        Vec3<T> operator - ( const T &operand ) const;
        Vec3<T> operator + ( const T &operand ) const;
        Vec3<T> operator * ( const T &operand ) const;
        Vec3<T> operator / ( const T &operand ) const;
        T calcMinDistance( SphereNode<T> *node ) const;
        T calcMaxDistance( SphereNode<T> *node ) const;
        T dot(          const Vec3<T> &other  ) const;
        Vec3<T> crossProduct( const Vec3<T> &other  ) const;
        T timesColumn( const Vec3<T> &operand ) const;
        Matrix3<T> timesRow( const Vec3<T> &operand ) const;

        // ---- Modification ------------------------------------------------------

        Vec3<T> &operator -=( const Vec3<T> &operand );
        Vec3<T> &operator +=( const Vec3<T> &operand );
        Vec3<T> &operator *=( const Vec3<T> &operand );
        Vec3<T> &operator *=( T              operand );
        Vec3<T> &operator /=( T              operand );
        void     normalize();

        /*************************************************************************/
        /************************** public Attributes ****************************/
        /*************************************************************************/

        /** The x component of the 3D-vector. */
        T x;

        /** The y component of the 3D-vector. */
        T y;

        /** The z component of the 3D-vector. */
        T z;
    };

    template <class T>
    std::ostream &operator<<( std::ostream &o, const Vec3<T> &v );
    template <class T>
    Vec3<T> operator *( const T &lhs, Vec3<T> rhs );

}


#endif //UNTITLED4_VEC3_H

