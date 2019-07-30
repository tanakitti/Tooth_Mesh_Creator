//
// Created by Tanakit Sachati on 2019-07-11.
//

#ifndef UNTITLED4_SPHERENODE_H
#define UNTITLED4_SPHERENODE_H

#include "Resources.h"
#include "Vec3.h"

namespace pd
{

/** The maximum degree of any SphereNode. This value can't just be changed here. */
    const int ISTTreeDegree = 4;

/**
 * @struct pd::Circle
 * @brief Simple data-structure to return a circle as an intersection result.
 */
    template <class T>
    struct Circle
    {
        /** The center of the circle. */
        Vec3<T> m_center;

        /** The normal that defines the plane in which the circle lies. */
        Vec3<T> m_normal;

        /** The radius of the circle. */
        T m_radius;
    };

    template <class T>
    class SphereNode
    {
    public:
        // ************************************************************************
        // *************************** public Functions ***************************
        // ************************************************************************

        // ---- De-/Constructor(s) ------------------------------------------------

        SphereNode();

        // ---- Queries -----------------------------------------------------------

        T distance( const SphereNode<T> &other );
        T distanceCapped( const SphereNode<T> &other );
        T calcFullVolume();
        T calcIntersectedVolume( const Vec3<T> &rayOrigin, const Vec3<T> &rayDir );
        bool intersectionCircle( Vec3<T> planeOrigin, Vec3<T> planeNormal, Circle<T> &result );

        // ************************************************************************
        // *************************** public Attributes **************************
        // ************************************************************************

        // ---- Geometry ----------------------------------------------------------

        /** The center of the sphere in 3D space. */
        Vec3<T> m_center;

        /** The radius of the sphere, the distance from the center to the surface. */
        T m_radius;

        // ---- Hierarchy ---------------------------------------------------------

        /** Flag that indicates if this is a leaf or an inner node with children. */
        bool m_isLeaf;

        /** The id of this node. */
        unsigned int m_nodeId;


        /** The id of this node's parent. */
        unsigned int m_parentId;

        /** The position of this node in its parent's children-array. */
        unsigned int m_childId;

        /** The pointers to potential child nodes, nullptr if not set. */
        SphereNode<T> *m_children[ISTTreeDegree];

        // ---- Collision ---------------------------------------------------------

        /** The density of the material that the sphere represents */
        T m_density;

    };

}

template <class T>
std::ostream& operator<<( std::ostream &o, const pd::SphereNode<T> &sn );

#endif //UNTITLED4_SPHERENODE_H
