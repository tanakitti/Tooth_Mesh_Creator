//
// Created by Tanakit Sachati on 2019-07-11.
//

//#include <cmath>

#include "../header/Matrix4.h"
#include "../header/Vec3.h"
#include <assert.h>
#include <math.h>

namespace pd
{

// ============================================================================
// ==== De-/Constructor =======================================================
// ============================================================================

/**
 * @brief Default constructor.
 */
    template <class T>
    Matrix4<T>::Matrix4()
    {
    }

/**
 * @brief Initialize all elements of this matrix with value.
 */
    template <class T>
    Matrix4<T>::Matrix4( T value )
    {
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) = value;
    }

/**
 * @brief Initialize all elements of this matrix with value.
 */
    template <class T>
    Matrix4<T>::Matrix4( T values[4][4] )
    {
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) = values[j][i];
    }

/**
 * @brief Initialize all elements of this matrix with value.
 */
    template <class T>
    Matrix4<T>::Matrix4( T values[16] )
    {
        int k=0;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) = values[k++];
    }

/**
 * @brief Initialize this matrix elements with a set of values.
 * @param a00 The value to be stored at M(0,0).
 * @param a01 The value to be stored at M(0,1).
 * @param a02 The value to be stored at M(0,2).
 * @param a03 The value to be stored at M(0,3).
 * @param a10 The value to be stored at M(1,0).
 * @param a11 The value to be stored at M(1,1).
 * @param a12 The value to be stored at M(1,2).
 * @param a13 The value to be stored at M(1,3).
 * @param a20 The value to be stored at M(2,0).
 * @param a21 The value to be stored at M(2,1).
 * @param a22 The value to be stored at M(2,2).
 * @param a23 The value to be stored at M(2,3).
 * @param a30 The value to be stored at M(3,0).
 * @param a31 The value to be stored at M(3,1).
 * @param a32 The value to be stored at M(3,2).
 * @param a33 The value to be stored at M(3,3).
 */
    template <class T>
    Matrix4<T>::Matrix4( T a00, T a01, T a02, T a03,
                         T a10, T a11, T a12, T a13,
                         T a20, T a21, T a22, T a23,
                         T a30, T a31, T a32, T a33 )
    {
        (*this)(0,0) = a00; (*this)(0,1) = a01; (*this)(0,2) = a02; (*this)(0,3) = a03;
        (*this)(1,0) = a10; (*this)(1,1) = a11; (*this)(1,2) = a12; (*this)(1,3) = a13;
        (*this)(2,0) = a20; (*this)(2,1) = a21; (*this)(2,2) = a22; (*this)(2,3) = a23;
        (*this)(3,0) = a30; (*this)(3,1) = a31; (*this)(3,2) = a32; (*this)(3,3) = a33;
    }

/** TODO */
    template <class T>
    Matrix4<T>::Matrix4( Vec3<T> axis, T angle )
    {
        T c = cosf(angle);
        T s = sinf(angle);
        T t = 1.f - c;
        (*this)(0,0) = c + axis.x*axis.x*t;
        (*this)(1,1) = c + axis.y*axis.y*t;
        (*this)(2,2) = c + axis.z*axis.z*t;

        T tmp1 = axis.x*axis.y*t;
        T tmp2 = axis.z*s;
        (*this)(1,0) = tmp1 + tmp2;
        (*this)(0,1) = tmp1 - tmp2;
        tmp1 = axis.x*axis.z*t;
        tmp2 = axis.y*s;
        (*this)(2,0) = tmp1 - tmp2;
        (*this)(0,2) = tmp1 + tmp2;
        tmp1 = axis.y*axis.z*t;
        tmp2 = axis.x*s;
        (*this)(2,1) = tmp1 + tmp2;
        (*this)(1,2) = tmp1 - tmp2;
        (*this)(3,3) = 1;
    }

// ============================================================================
// ==== Operators =============================================================
// ============================================================================

/**
 * @brief Check if all elements are the same.
 * @param M The other matrix that will be tested for equality.
 * @return If this and the other matrix have element-wise equality.
 * @retval true   If all elements of both matrices are the same.
 * @retval false  If any elements of the two matrices differ.
 */
    template <class T>
    bool Matrix4<T>::operator ==( const Matrix4<T> &M ) const
    {
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                if ((*this)(i,j) != M(i,j))
                    return false;
        return true;
    }

/**
 * @brief Check if the two matrices are unequal.
 * @param M The other matrix that will be tested for inequality.
 * @return If this and the other matrix are unequal element-wise.
 * @retval true   If any of the elements of both matrices differ.
 * @retval false  If all elements of both matrices are the same.
 */
    template <class T>
    bool Matrix4<T>::operator !=( const Matrix4<T> &M ) const
    {
        return !(*this == M);
    }

/**
 * @brief Access a certain element of this matrix.
 * @return The element of this matrix at row row & column column as a reference.
 *
 * @note This was only implemented to have a easier time switching between
 *       column=major & row=major order, however in kernels it is not used for
 *       performance reasons.
 */
    template <class T>
            T &Matrix4<T>::operator ()( unsigned int row, unsigned int column )
    {
        // Row-Major-Order
        //return data[row][column];

        // Column-Major-Order
        return data[column][row];
    }

/**
 * @brief Access a certain element of this matrix.
 * @return The value of the element of this matrix at row row & column column.
 *
 * @note This was only implemented to have a easier time switching between
 *       column=major & row=major order, however in kernels it is not used for
 *       performance reasons.
 */
    template <class T>
            T Matrix4<T>::operator ()( unsigned int row, unsigned int column ) const
    {
        // Row-Major-Order
        //return data[row][column];

        // Column-Major-Order
        return data[column][row];
    }

/**
 * @brief Return the transpose of this matrix.
 * @return The transpose of this matrix.
 */
    template <class T>
    const Matrix4<T> Matrix4<T>::operator ~() const
    {
        Matrix4<T> result;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                result(i,j) = (*this)(j,i);
        return result;
    }

/**
 * @brief Get the sum of this matrix plus the scalar parameter rhs.
 * @param rhs The scalar that will be added from the right side to the matrix.
 * @return The sum of this matrix and scalar rhs.
 */
    template <class T>
    const Matrix4<T> Matrix4<T>::operator +( const T &rhs ) const
    {
        Matrix4<T> result;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                result(i,j) = (*this)(i,j) + rhs;
        return result;
    }

/**
 * @brief Call the Matrix4<T>::operator+(const T&) const. '+' is commutative.
 */
    template <class T>
    const Matrix4<T> operator +( const T &lhs, const Matrix4<T> &rhs )
    {
        return rhs + lhs;
    }

    template const Matrix4<float> operator +( const float &lhs, const Matrix4<float> &rhs );
    template const Matrix4<double> operator +( const double &lhs, const Matrix4<double> &rhs );

/**
 * @brief Get the product of this matrix and the scalar parameter rhs.
 * @param rhs The scalar that will be multiplied from the right side.
 * @return The product of this matrix and scalar rhs.
 */
    template <class T>
    const Matrix4<T> Matrix4<T>::operator *( const T &rhs ) const
    {
        Matrix4<T> result;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                result(i,j) = (*this)(i,j) * rhs;
        return result;
    }

/**
 * @brief Call the Matrix4<T>::operator*(const T&) const. '*' is commutative.
 */
    template <class T>
    const Matrix4<T> operator *( const T &lhs, const Matrix4<T> &rhs )
    {
        return rhs * lhs;
    }

    template const Matrix4<float> operator *( const float &lhs, const Matrix4<float> &rhs );
    template const Matrix4<double> operator *( const double &lhs, const Matrix4<double> &rhs );

/**
 * @brief Return the sum of this matrix and the matrix parameter rhs.
 * @param rhs The matrix that will be added from the right side to this matrix.
 * @return The sum of this matrix and the matrix rhs.
 */
    template <class T>
    const Matrix4<T> Matrix4<T>::operator +( const Matrix4<T> &rhs ) const
    {
        Matrix4<T> result;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                result(i,j) = (*this)(i,j) + rhs(i,j);
        return result;
    }

/**
 * @brief Get the product of this matrix and the matrix parameter rhs.
 * @param rhs The matrix that will be multiplied from the right hand side.
 * @return The product of this matrix and the other matrix rhs.
 */
    template <class T>
    const Matrix4<T> Matrix4<T>::operator *( const Matrix4<T> &rhs ) const
    {
        Matrix4<T> result;
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<4; j++)
            {
                T sum = 0;
                for (int k=0; k<4; k++)
                    sum += (*this)(i,k) * rhs(k,j);
                result(i,j) = sum;
            }
        }
        return result;
    }

// ============================================================================
// ==== Modification ==========================================================
// ============================================================================

/**
 * @brief Overwrite the matrix with the 4x4 identity matrix.
 */
    template <class T>
    void Matrix4<T>::identity()
    {
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                i == j ? (*this)(i,j) = 1.0f : (*this)(i,j) = 0.0f;
    }

/**
 * @brief Ortho-normalize this matrix.
 *
 * This is done by taking two of the axis and generating a third one by cross=
 * product, so it is guaranteed to be perpendicular to the first one. Now we
 * replace the second axis by calculating the crossproduct of the first and
 * third axis, so now we have 3 axis that are guaranteed to be perpendicular to
 * one another.
 *
 * @warning This assumes that this matrix is a rotation matrix.
 * @warning This is used only for the virtual tool simulation of torque, so it
 *          was not thoroughly tested. It can be tested by checking if the
 *          determinant of the matrix is very close to 1 afterwards.
 */
    template <class T>
    void Matrix4<T>::orthonormalize()
    {
        Vec3<T> x( (*this)(0,0), (*this)(1,0), (*this)(2,0) );
        Vec3<T> y( (*this)(0,1), (*this)(1,1), (*this)(2,1) );
        Vec3<T> z;

        x.normalize();
        z = x.crossProduct( y );
        z.normalize();
        y = z.crossProduct( x );
        y.normalize();

        (*this)(0,0) = x.x; (*this)(0,1) = y.x; (*this)(0,2) = z.x;
        (*this)(1,0) = x.y; (*this)(1,1) = y.y; (*this)(1,2) = z.y;
        (*this)(2,0) = x.z; (*this)(2,1) = y.z; (*this)(2,2) = z.z;
    }

/**
 * @brief Invert this matrix.
 *
 * First checks if this matrix is invertible by checking if it's determinant is
 * not zero. If it is invertible, the inverse matrix is calculated and all
 * values of this matrix are replaced by the inverse.
 *
 * @return false if matrix not invertible, otherwise true.
 */
    template <class T>
    bool Matrix4<T>::invert()
    {
        Matrix4<T> result;
        T rDet;
        T a00, a01, a02, a03,
                a10, a11, a12, a13,
                a20, a21, a22, a23,
                a30, a31, a32, a33;

        a00=(*this)(0,0); a01=(*this)(0,1); a02=(*this)(0,2); a03=(*this)(0,3);
        a10=(*this)(1,0); a11=(*this)(1,1); a12=(*this)(1,2); a13=(*this)(1,3);
        a20=(*this)(2,0); a21=(*this)(2,1); a22=(*this)(2,2); a23=(*this)(2,3);
        a30=(*this)(3,0); a31=(*this)(3,1); a32=(*this)(3,2); a33=(*this)(3,3);

        rDet = det();

        if (rDet == 0)
        {
            this->identity();
            return false;
        }

        rDet = 1/rDet;

        result(0,0) =  det3( a11, a12, a13, a21, a22, a23, a31, a32, a33 ) * rDet;
        result(0,1) = -det3( a01, a02, a03, a21, a22, a23, a31, a32, a33 ) * rDet;
        result(0,2) =  det3( a01, a02, a03, a11, a12, a13, a31, a32, a33 ) * rDet;
        result(0,3) = -det3( a01, a02, a03, a11, a12, a13, a21, a22, a23 ) * rDet;

        result(1,0) = -det3( a10, a12, a13, a20, a22, a23, a30, a32, a33 ) * rDet;
        result(1,1) =  det3( a00, a02, a03, a20, a22, a23, a30, a32, a33 ) * rDet;
        result(1,2) = -det3( a00, a02, a03, a10, a12, a13, a30, a32, a33 ) * rDet;
        result(1,3) =  det3( a00, a02, a03, a10, a12, a13, a20, a22, a23 ) * rDet;

        result(2,0) =  det3( a10, a11, a13, a20, a21, a23, a30, a31, a33 ) * rDet;
        result(2,1) = -det3( a00, a01, a03, a20, a21, a23, a30, a31, a33 ) * rDet;
        result(2,2) =  det3( a00, a01, a03, a10, a11, a13, a30, a31, a33 ) * rDet;
        result(2,3) = -det3( a00, a01, a03, a10, a11, a13, a20, a21, a23 ) * rDet;

        result(3,0) = -det3( a10, a11, a12, a20, a21, a22, a30, a31, a32 ) * rDet;
        result(3,1) =  det3( a00, a01, a02, a20, a21, a22, a30, a31, a32 ) * rDet;
        result(3,2) = -det3( a00, a01, a02, a10, a11, a12, a30, a31, a32 ) * rDet;
        result(3,3) =  det3( a00, a01, a02, a10, a11, a12, a20, a21, a22 ) * rDet;

        *this = result;

        return true;
    }

/**
 * @brief Overwrite the elements of this matrix with those of matrix M.
 * @param M The matrix to take all values from.
 */
    template <class T>
    void Matrix4<T>::copyVals( Matrix4<T> M )
{
    for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    (*this)(i,j) = M(i,j);
}

// ============================================================================
// ==== Analysis ==============================================================
// ============================================================================

/**
 * @brief Calculates the determinant of the 3x3 sub-matrix of this matrix.
 * @return The resulting determinant.
 */
template <class T>
        T Matrix4<T>::det3() const
{
    return (*this)(0,0) * (*this)(1,1) * (*this)(2,2)
           + (*this)(0,1) * (*this)(1,2) * (*this)(2,0)
           + (*this)(0,2) * (*this)(1,0) * (*this)(2,1)
           - (*this)(0,2) * (*this)(1,1) * (*this)(2,0)
           - (*this)(0,1) * (*this)(1,0) * (*this)(2,2)
           - (*this)(0,0) * (*this)(1,2) * (*this)(2,1);
}

/**
 * @brief Calculates the determinant according to given parameters.
 *
 * Treats the given parameters as the values of a 3x3 matrix and calculates
 * the determinant accordingly and returns it.
 *
 * @param a11 The value of the element at first row and first column
 * @param a12 The value of the element at first row and second column
 * @param a13 The value of the element at first row and third column
 * @param a21 The value of the element at second row and first column
 * @param a22 The value of the element at second row and second column
 * @param a23 The value of the element at second row and third column
 * @param a31 The value of the element at third row and first column
 * @param a32 The value of the element at third row and second column
 * @param a33 The value of the element at third row and third column
 * @return The resulting determinant
 */
template <class T>
        T Matrix4<T>::det3( T a11, T a12, T a13,
T a21, T a22, T a23,
        T a31, T a32, T a33 )
{
return a11 * a22 * a33
+ a12 * a23 * a31
+ a13 * a21 * a32
- a13 * a22 * a31
- a12 * a21 * a33
- a11 * a23 * a32;
}

/**
 * @brief Calculates the determinant of this matrix.
 * @return The resulting determinant.
 */
template <class T>
        T Matrix4<T>::det() const
{
    T a1, a2, a3, a4,
            b1, b2, b3, b4,
            c1, c2, c3, c4,
            d1, d2, d3, d4;

    a1 = (*this)(0,0); a2 = (*this)(0,1); a3 = (*this)(0,2); a4 = (*this)(0,3);
    b1 = (*this)(1,0); b2 = (*this)(1,1); b3 = (*this)(1,2); b4 = (*this)(1,3);
    c1 = (*this)(2,0); c2 = (*this)(2,1); c3 = (*this)(2,2); c4 = (*this)(2,3);
    d1 = (*this)(3,0); d2 = (*this)(3,1); d3 = (*this)(3,2); d4 = (*this)(3,3);

    return a1 * det3( b2, b3, b4, c2, c3, c4, d2, d3, d4 )
           - b1 * det3( a2, a3, a4, c2, c3, c4, d2, d3, d4 )
           + c1 * det3( a2, a3, a4, b2, b3, b4, d2, d3, d4 )
           - d1 * det3( a2, a3, a4, b2, b3, b4, c2, c3, c4 );
}

/**
 * @brief Prints the content of the matrix.
 */
template <class T>
void Matrix4<T>::print() const
{
    printf( "[Matrix4<T>] %p\n", static_cast<const void*>(this) );
    for (int i=0; i<4; i++)
    {
        printf( "[ " );
        for (int j=0; j<4; j++)
            printf( "%f ", (*this)(i,j) );
        printf( "]\n" );
    }
    printf( "\n" );
}

/**
 * @brief TODO
 */
template <class T>
bool Matrix4<T>::toEulerAnglesZYX (T &ex, T &ey, T &ez) const
{
    // rot =  cy*cz           cz*sx*sy-cx*sz  cx*cz*sy+sx*sz
    //        cy*sz           cx*cz+sx*sy*sz -cz*sx+cx*sy*sz
    //       -sy              cy*sx           cx*cy

    ey = asinf(-data[2][0]);
    if (ey < __PD__PIHALF)
    {
        if (ey > -__PD__PIHALF)
        {
            ex = atan2f(data[1][0],data[0][0]);
            ez = atan2f(data[2][1],data[2][2]);
            return true;
        }
        else
        {
            // WARNING.  Not a unique solution.
            T fRmY = atan2f(-data[0][1],data[0][2]);
            ez = 0.f;  // any angle works
            ex = ez - fRmY;
            return false;
        }
    }
    else
    {
        // WARNING.  Not a unique solution.
        T fRpY = atan2f(-data[0][1],data[0][2]);
        ez = 0.f;  // any angle works
        ex = fRpY - ez;
        return false;
    }
}

/**
 * @brief TODO
 */
template <class T>
bool Matrix4<T>::toEulerAnglesZXY(T &ex, T &ey, T &ez) const
{
    // rot =  cy*cz-sx*sy*sz -cx*sz           cz*sy+cy*sx*sz
    //        cz*sx*sy+cy*sz  cx*cz          -cy*cz*sx+sy*sz
    //       -cx*sy           sx              cx*cy

    ey = asinf(-data[2][1]);
    if (ey < __PD__PIHALF)
    {
        if (ey > -__PD__PIHALF)
        {
            ex = atan2f(-data[0][1],data[1][1]);
            ez = atan2f(-data[2][0],data[2][2]);
            return true;
        }
        else
        {
            // WARNING.  Not a unique solution.
            T fRmY = atan2f(data[0][2],data[0][0]);
            ez = 0.f;  // any angle works
            ex = ez - fRmY;
            return false;
        }
    }
    else
    {
        // WARNING.  Not a unique solution.
        T fRpY = atan2f(data[0][2],data[0][0]);
        ez = 0.f;  // any angle works
        ex = fRpY - ez;
        return false;
    }
}

/**
 * @brief TODO
 */
template <class T>
bool Matrix4<T>::toEulerAnglesYZX( T &ex, T &ey, T &ez ) const
{
    ey = asinf(data[1][0]);
    if (ey < __PD__PIHALF)
    {
        if (ey > -__PD__PIHALF)
        {
            ex = atan2f(-data[2][0],data[0][0]);
            ez = atan2f(-data[1][2],data[1][1]);
            return true;
        }
        else
        {
            // WARNING.  Not a unique solution.
            T fRmY = atan2f(data[2][1],data[2][2]);
            ez = 0.f; // any angle works
            ex = ez - fRmY;
            return false;
        }
    }
    else
    {
        // WARNING.  Not a unique solution.
        T fRpY = atan2f(data[2][1],data[2][2]);
        ez = 0.f; // any angle works
        ex = fRpY - ez;
        return false;
    }
}

/**
 * @brief Return column vector at column col.
 */
template <class T>
        Vec3<T> Matrix4<T>::getColumn( int col ) const
{
    return Vec3<T>( (*this)(0,col), (*this)(1,col), (*this)(2,col) );
}

/**
 * @brief Return row vector at row row.
 */
template <class T>
        Vec3<T> Matrix4<T>::getRow( int row ) const
{
    return Vec3<T>( (*this)(row,0), (*this)(row,1), (*this)(row,2) );
}

/**
 * @brief Return translation vector of this matrix.
 */
template <class T>
Vec3<T> Matrix4<T>::extractPosition() const
{
    return Vec3<T>( (*this)(0,3), (*this)(1,3), (*this)(2,3) );
}

/**
 * @brief Return matrix that only does the translation of this.
 */
template <class T>
Matrix4<T> Matrix4<T>::reduceToTranslation() const
{
    pd::Matrix4<T> M = Matrix4<T>::IDENTITY();
    for (int row=0; row<3; row++)
        M(row,3) = (*this)(row,3);
    return M;
}

/**
 * @brief Return euler angles (XYZ) of this matrix.
 */
template <class T>
Vec3<T> Matrix4<T>::extractEulerAngles() const
{
    auto M = reduceToRotation();
    T ex, ey, ez;
    ey = asinf(-M(2,0));
    if (ey < __PD__PIHALF)
    {
        if (ey > -__PD__PIHALF)
        {
            ex = atan2f(M(1,0),M(0,0));
            ez = atan2f(M(2,1),M(2,2));
        }
        else
        {
            // WARNING.  Not a unique solution.
            T fRmY = atan2f(-M(0,1),M(0,2));
            ez = 0.f;  // any angle works
            ex = ez - fRmY;
        }
    }
    else
    {
        // WARNING.  Not a unique solution.
        T fRpY = atan2f(-M(0,1),M(0,2));
        ez = 0.f;  // any angle works
        ex = fRpY - ez;
    }
    return Vec3<T>( ex, ey, ez );
}

/**
 * @brief Return matrix that only does the rotation of this.
 */
template <class T>
Matrix4<T> Matrix4<T>::reduceToRotation() const
{
    pd::Matrix4<T> R = Matrix4<T>::IDENTITY();
    for (int col=0; col<3; col++)
        for (int row=0; row<3; row++)
            R(row,col) = (*this)(row,col) / getColumn(col).calcLength();
    return R;
}

/**
 * @brief Return scale of this matrix.
 */
template <class T>
Vec3<T> Matrix4<T>::extractScale() const
{
    auto sx = getColumn(0).calcLength();
    auto sy = getColumn(1).calcLength();
    auto sz = getColumn(2).calcLength();
    auto s = Vec3<T>( sx, sy, sz );
    return s;
}

/**
 * @brief Return matrix that only does the scaling of this.
 */
template <class T>
Matrix4<T> Matrix4<T>::reduceToScale() const
{
    pd::Matrix4<T> S = Matrix4<T>::IDENTITY();
    for (int i=0; i<3; i++)
        S(i,i) = getColumn(i).calcLength();
    return S;
}

// ============================================================================
// ==== Multiplication ========================================================
// ============================================================================

/**
 * @brief Multiplies this matrix with vector src and stores the result in dst.
 * @param src The vector to multiply this matrix with.
 * @param dst The pointer to the vector to save the result in.
 */
template <class T>
void Matrix4<T>::multFullMatrixPnt( Vec3<T> src, Vec3<T> &dst ) const
{
T w = src.x*(*this)(3,0) + src.y*(*this)(3,1) + src.z*(*this)(3,2) + (*this)(3,3);

if (w == 0)
{
dst.x = 0;
dst.y = 0;
dst.z = 0;
printf( "FAIL: Can't multiply non-homogeneous 4x4 matrix!\n" );
this->print();
//assert(0);
}

w = 1.f/w;

dst.x = (src.x*(*this)(0,0) + src.y*(*this)(0,1) + src.z*(*this)(0,2) + (*this)(0,3)) * w;
dst.y = (src.x*(*this)(1,0) + src.y*(*this)(1,1) + src.z*(*this)(1,2) + (*this)(1,3)) * w;
dst.z = (src.x*(*this)(2,0) + src.y*(*this)(2,1) + src.z*(*this)(2,2) + (*this)(2,3)) * w;
}

/**
 * @brief Multiplies the upper left 3x3 part of this matrix with vector src and
 *        store the result in dst.
 * @param src The vector to multiply this matrix with.
 * @param dst The pointer to the vector to save the result in.
 */
template <class T>
void Matrix4<T>::multRotMatrixPnt( Vec3<T> src, Vec3<T> &dst ) const
{
dst.x = src.x*(*this)(0,0) + src.y*(*this)(0,1) + src.z*(*this)(0,2);
dst.y = src.x*(*this)(1,0) + src.y*(*this)(1,1) + src.z*(*this)(1,2);
dst.z = src.x*(*this)(2,0) + src.y*(*this)(2,1) + src.z*(*this)(2,2);
}

template class Matrix4<float>;
template class Matrix4<double>;

}