//
// Created by Tanakit Sachati on 2019-07-11.
//

#include "../header/SphereNode.h"
#include <math.h>
#include <float.h>
#include <math.h>

namespace pd
{

// ============================================================================
// ==== De-/Constructor(s) ====================================================
// ============================================================================

/**
 * @brief Initialize all attributes with sensible default values.
 */
    template <class T>
    SphereNode<T>::SphereNode()
            : m_center( 0,0,0 )
            , m_radius( 0 )
            , m_isLeaf( false )
            , m_nodeId( 0 )
    {
    }

// ============================================================================
// ==== Queries ===============================================================
// ============================================================================

/**
 * @brief Calculate the euclidean distance between the two SphereNode's.
 *
 * The distance between the two SphereNode's is simply the distance of the
 * two centers minus both radii. This gives the distance between the two sphere
 * surfaces on the straight line spanned from one center to the other center.
 *
 * @param other The other SphereNode, to which we calculate the distance.
 * @return The euclidean distance between this SphereNode and other SphereNode.
 */
    template <class T>
    T SphereNode<T>::distance( const SphereNode<T> &other )
    {
        return (m_center-other.m_center).calcLength() - m_radius - other.m_radius;
    }

/**
 * @brief Return the distance but don't allow negative values.
 * @see SphereNode::distance(const SphereNode&) for details.
 */
    template <class T>
    T SphereNode<T>::distanceCapped( const SphereNode<T> &other )
    {
        T dist = distance( other );

        if (dist > 0)
            return dist;

        return 0;
    }

/**
 * @brief Return the volume of this sphere.
 *
 * This is simply done by using the classic formula ¾πr³.
 *
 * @return The volume of this SphereNode's sphere geometry.
 */
    template <class T>
    T SphereNode<T>::calcFullVolume()
    {
        return 4.f/3.f * __PD__PI * m_radius * m_radius * m_radius;
    }

/**
 * @brief Calculate the volume that is intersected by the sphere and plane.
 *
 * We firstly convert the point + normal plane to Hesse normal form by
 * calculating the distance of the plane from origin (stored in D). With that
 * we calculate the signed distance of the sphere's center to the plane.
 *
 * - If the distance of the center to the plane is less than the radius of the
 * sphere, it means the plane intersects the sphere. In that case we calculate
 * the volume of a spherical cap that is of height h=r-d (r being the radius,
 * and d being the distance), given by ⅓πh²(3r-h).
 * - In case the distance is positive and bigger than r, we know the sphere is
 *   completely in front of the plane, without intersection, so we the volume
 *   is 0.
 * - If the distance is negative, but bigger than r, we know the sphere is
 *   completely behind the plane, so the full volume of the sphere is returned.
 *
 * @param planeOrigin A point on the plane.
 * @param planeNormal The normal of the plane.
 * @return The volume of the sphere that is behind the plane.
 *
 * @warning planeNormal is assumed to be normalized.
 */
    template <class T>
    T SphereNode<T>::calcIntersectedVolume( const Vec3<T> &planeOrigin, const Vec3<T> &planeNormal )
    {
        // Convert to Hesse normal form
        T D = -planeNormal.dot( planeOrigin );
        T signed_distance = planeNormal.dot( m_center ) + D;

        // Sphere intersects plane? -> Calculate spherical cap volume
        if (fabsf(signed_distance) < m_radius)
        {
            T height = m_radius - signed_distance;
            return ((__PD__PI*height*height)/3.f) * (3.f*m_radius - height);
        }

        // Sphere in front of plane? -> No volume
        if (signed_distance >= 0)
            return 0;

        // Sphere behind plane? -> Use complete sphere volume
        return 4.f/3.f * __PD__PI * m_radius * m_radius * m_radius;
    }

/**
 * @brief Calculate the intersection circle between this sphere and given plane.
 *
 * Convert the plane to Hesse normal form and calculate the signed distance,
 * if the distance is between less than m_radius, the plane intersects the
 * sphere. In that case we calculate the intersection circle in Circle result
 * and return true. Otherwise result is untouched and false is returned.
 *
 * @param planeOrigin A point on the plane.
 * @param planeNormal The normal of the plane.
 * @param result The resulting intersection circle.
 * @return Whether a intersection was found.
 * @retval true   If a intersection was found.
 * @retval false  If there was no intersection.
 *
 * @warning planeNormal is assumed to be normalized.
 */
    template <class T>
    bool SphereNode<T>::intersectionCircle( Vec3<T> planeOrigin, Vec3<T> planeNormal, Circle<T> &result )
    {
        // Convert to Hesse normal form
        T D = -planeNormal.dot( planeOrigin );
        T signed_distance = planeNormal.dot( m_center ) + D;

        if (signed_distance < -m_radius || m_radius < signed_distance)
            return false;

        // Calculate intersection circle
        result.m_radius = sqrtf((m_radius*m_radius) - (signed_distance*signed_distance));
        result.m_center = m_center - (planeNormal*signed_distance);
        result.m_normal = planeNormal;
        return true;
    }

}

// ============================================================================
// ==== Misc. =================================================================
// ============================================================================

/**
 * @brief Print a SphereNode by printing each of the geometry attributes.
 */
template <class T>
std::ostream& operator<<(std::ostream &o, const pd::SphereNode<T> &sn)
{
    o << "SphereNode x=" << sn.m_center.x << " y=" << sn.m_center.y << " z=" << sn.m_center.z << " r=" << sn.m_radius;
    return o;
}

template class pd::SphereNode<float>;
template class pd::SphereNode<double>;
