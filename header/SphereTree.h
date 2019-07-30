//
// Created by Tanakit Sachati on 2019-07-11.
//

#ifndef UNTITLED4_SPHERETREE_H
#define UNTITLED4_SPHERETREE_H

#include <vector>
#include <map>
#include <float.h>

#include "Resources.h"
#include "Matrix4.h"
#include "SphereNode.h"
#include "Vec3.h"

namespace pd
{

    template <class T>
    class SphereTree
    {
    public:
        // ************************************************************************
        // *************************** public Functions ***************************
        // ************************************************************************

        // ---- De-/Constructor(s) ------------------------------------------------

        SphereTree();
        SphereTree( const char *fileName, int version=1, const Matrix4<T> M=Matrix4<T>::IDENTITY(), T density=1 );
        virtual ~SphereTree();

        // ---- Import ------------------------------------------------------------

        void convertToNewVersion();
        void generateRegions( Vec3<T> center, Matrix4<T> transform );

        // ---- Export ------------------------------------------------------------

        void exportHierarchyLegacy( const char *fileName );
        void exportHierarchy( const char *fileName );

        void exportCcsObj( const char *objFileName, const char *fileName );
        bool rayTriangleTest(Vec3<T> O, Vec3<T> d, T dMax, Vec3<T> v0, Vec3<T> v1, Vec3<T> v2) const;
        void accumulateNeighboursDown( SphereNode<T> *leaf, SphereNode<T> *root, std::vector<std::vector<unsigned int> > *neighbours2d ) const;

        // ---- Queries -----------------------------------------------------------

        bool isChildOf( SphereNode<T> *child, int parentId ) const;

        // ---- ISTGraph ----------------------------------------------------------

        void depthFirstSearch( unsigned int c, std::vector<std::vector<unsigned int> > *neighbourIds, std::map<unsigned int, int> *id_to_cc, SphereNode<T> *node );
        std::vector<std::vector<unsigned int> > computeConnectedComponents( std::vector<std::vector<unsigned int> > neighbourIds );
        void insertBridge( std::vector<std::vector<unsigned int> > *cc, std::vector<std::vector<unsigned int> > *neighbours );

        // ************************************************************************
        // *************************** public Attributes **************************
        // ************************************************************************

        // ---- Geometry ----------------------------------------------------------

        /** The inverse of the transformation of the IST (to transform PCDs with). */
        Matrix4<T> m_invTransform;

        /** The volume of the IST, simply the sum of each sphere volume. */
        T m_volume = 0.0;

        /**
         * The center of mass of the IST. This is simply the sum of all SphereNode
         * centers, weighted by their volume and later normalized by the overall
         * volume.
         */
        Vec3<T> m_centerOfMass = Vec3<T>(0, 0, 0);

        /** TODO */
        Vec3<T> m_aabbLo = Vec3<T>(FLT_MAX, FLT_MAX, FLT_MAX);

        /** TODO */
        Vec3<T> m_aabbHi = Vec3<T>(-FLT_MAX, -FLT_MAX, -FLT_MAX);

        /** TODO */
        Vec3<T> m_aabbCorners[8];

        /**
         * The inertia tensor of the IST. This is used for the physics simulation
         * to apply the different torques to, to visualize them.
         *
         * @see SphereTree::convertToNewVersion() for details on creation.
         */
        Matrix4<T> m_inertiaTensor = Matrix4<T>::IDENTITY();

        /**
         * The inverse of the inertia tensor of the IST. This is used for the
         * physics simulation to apply the different torques to, to visualize them.
         *
         * @see SphereTree::convertToNewVersion() for details on creation.
         */
        Matrix4<T> m_invInertiaTensor = Matrix4<T>::IDENTITY();

        // ---- Hierarchy ---------------------------------------------------------

        bool m_hasHierarchy = false;

        /** Amount of SphereNode's in m_nodes & m_dev_nodes (inner nodes & leaves). */
        unsigned int m_countSpheres = 0;

        /** The amount of inner node SphereNode's. */
        unsigned int m_countNodes = 0;

        /** The amount of leaf SphereNode's in m_leaves & m_dev_leaves. */
        unsigned int m_countLeaves = 0;

        /** Pointer to the root SphereNode hierarchy. */
        SphereNode<T> *m_rootNode = nullptr;

        /** The memory that holds all the SphereNode's inner nodes. */
        SphereNode<T> *m_nodes = nullptr;

        /**
         * The memory that holds all the SphereNode's leaves.
         *
         * @warning m_nodeId is adjusted in these SphereNode's, so that they match
         *          their position in this array m_leaves, not m_nodes.
         * @see SphereTree::convertToNewVersion() for details on how the nodeIds
         *      are adjusted for this.
         */
        SphereNode<T> *m_leaves = nullptr;

        /** @see SphereTree::m_rootNode, since this is analog, but in device memory. */
        SphereNode<T> *m_dev_rootNode = nullptr;

        /** @see SphereTree::m_nodes, since this is analog, but in device memory. */
        SphereNode<T> *m_dev_nodes = nullptr;

        /** @see SphereTree::m_leaves, since this is analog, but in device memory. */
        SphereNode<T> *m_dev_leaves = nullptr;

    private:
        // ---- Import ------------------------------------------------------------
        void loadHierarchyLegacy( const char *fileName, Matrix4<T> transform=Matrix4<T>::IDENTITY, T density=1 );
        template <int N=0>
        void loadHierarchy( const char *fileName, Matrix4<T> transform=Matrix4<T>::IDENTITY, T density=1 );
    };

}

#endif //UNTITLED4_SPHERETREE_H
