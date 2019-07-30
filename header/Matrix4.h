//
// Created by Tanakit Sachati on 2019-07-11.
//

#ifndef UNTITLED4_MATRIX4_H
#define UNTITLED4_MATRIX4_H

#include "Resources.h"
#include "Matrix4.h"
#include "Vec3.h"

namespace pd
{
    template <class T>
    class Matrix4
    {
    public:
        /*************************************************************************/
        /*************************** public Functions ****************************/
        /*************************************************************************/

        // ---- De-/Constructor ---------------------------------------------------

        Matrix4<T>();
        Matrix4<T>( T value );
        Matrix4<T>( T values[4][4] );
        Matrix4<T>( T values[16] );

        Matrix4<T>( T a00, T a01, T a02, T a03,
                    T a10, T a11, T a12, T a13,
                    T a20, T a21, T a22, T a23,
                    T a30, T a31, T a32, T a33 );

        Matrix4<T>(Vec3<T> axis, T angle );

        // ---- Static ------------------------------------------------------------

        static Matrix4<T> ZERO() { return Matrix4<T>(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0); };
        static Matrix4<T> IDENTITY() { return Matrix4<T>(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1); };

        // ---- Operators ---------------------------------------------------------

        bool operator==( const Matrix4<T> &M ) const;
        bool operator!=( const Matrix4<T> &M ) const;

        T &operator ()( unsigned int row, unsigned int column );
        T  operator ()( unsigned int row, unsigned int column ) const;

        const Matrix4<T> operator ~() const;

        const Matrix4<T> operator +( const T &rhs ) const;
        const Matrix4<T> operator +( const Matrix4<T> &rhs ) const;

        const Matrix4<T> operator *( const T &rhs ) const;
        const Matrix4<T> operator *( const Matrix4<T> &rhs ) const;

        // ---- Modification ------------------------------------------------------

        void identity();
        void orthonormalize();
        bool invert();
        void copyVals( Matrix4<T> M );

        // ---- Analysis ----------------------------------------------------------

        void print() const;

        T det() const;
        T det3() const;
        static T det3( T a00, T a01, T a02,
        T a10, T a11, T a12,
                T a20, T a21, T a22 );

        bool toEulerAnglesZYX( T &ex, T &ey, T &ez ) const;
        bool toEulerAnglesZXY( T &ex, T &ey, T &ez ) const;
        bool toEulerAnglesYZX( T &ex, T &ey, T &ez ) const;

        Vec3<T> getColumn( int col ) const;
        Vec3<T> getRow( int row ) const;

        Vec3<T> extractPosition() const;
        Vec3<T> extractEulerAngles() const;
        Vec3<T> extractScale() const;

        Matrix4<T> reduceToTranslation() const;
        Matrix4<T> reduceToRotation() const;
        Matrix4<T> reduceToScale() const;

        // ---- Multiplication ----------------------------------------------------

        void multFullMatrixPnt( Vec3<T> src, Vec3<T> &dst ) const;
        void multRotMatrixPnt(  Vec3<T> src, Vec3<T> &dst ) const;

        /*************************************************************************/
        /************************** public Attributes ****************************/
        /*************************************************************************/

        /** The underlying memory where all 16 elements are stored. */
        T data[4][4];
    };

/** Matrix addition with scalar, with matrix as rhs. */
    template <class T>
    const Matrix4<T> operator +( const T &lhs, const Matrix4<T> &rhs );

/** Matrix multiplication with scalar, with matrix as rhs. */
    template <class T>
    const Matrix4<T> operator *( const T &lhs, const Matrix4<T> &rhs );

}

#endif //UNTITLED4_MATRIX4_H
