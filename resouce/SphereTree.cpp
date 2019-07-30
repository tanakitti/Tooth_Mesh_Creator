//
// Created by Tanakit Sachati on 2019-07-11.
//

#include "../header/SphereTree.h"

#include "../header/Resources.h"

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <queue>
#include <set>
#include <iomanip>
#include <float.h>
#include <stdio.h>
#include "../header/Bankhelper.h"

#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#include <fstream>
#include <math.h>

#define GetCurrentDir getcwd
#endif

namespace pd
{

// ============================================================================
// ==== De-/Constructor(s) ====================================================
// ============================================================================

/**
* @brief Initialize all members with sensible default values.
*/
template <class T>
SphereTree<T>::SphereTree()
{
}

/**
* @brief Calls the default constructor, then parses the given .tree-file.
*
* @param fileName The path to the .tree-file that holds the IST & hierarchy.
* @param version An number that indicates whether it's a legacy .tree-file or
*                one that was created with SphereTree::exportHierarchy().
*                (0 means legacy, anything else means modern)
*/
template <class T>
SphereTree<T>::SphereTree( const char *fileName, int version, const Matrix4<T> M, T density )
	: SphereTree<T>()
{
	if (version == 0) {
		loadHierarchyLegacy(fileName, M, density);
        convertToNewVersion();
    }
	else if (version == 1)
		loadHierarchy<1>(fileName, M, density);
	else if (version == 2)
		loadHierarchy<2>(fileName, M, density);

    // Bounding box
    Vec3<T> diff = m_aabbHi - m_aabbLo;
    m_aabbCorners[0] = m_aabbLo;
    m_aabbCorners[1] = m_aabbLo + (diff*Vec3<T>(0,0,1).normalizedCopy());
    m_aabbCorners[2] = m_aabbLo + (diff*Vec3<T>(0,1,0).normalizedCopy());
    m_aabbCorners[3] = m_aabbLo + (diff*Vec3<T>(0,1,1).normalizedCopy());
    m_aabbCorners[4] = m_aabbLo + (diff*Vec3<T>(1,0,0).normalizedCopy());
    m_aabbCorners[5] = m_aabbLo + (diff*Vec3<T>(1,0,1).normalizedCopy());
    m_aabbCorners[6] = m_aabbLo + (diff*Vec3<T>(1,1,0).normalizedCopy());
    m_aabbCorners[7] = m_aabbHi;

    // Printing
    std::cout << "[SphereTree] Done reading " << fileName << ", sphere count: " << m_countSpheres << std::endl;
    std::cout << "[SphereTree] bounding box: " << m_aabbLo << " -> " << m_aabbHi << std::endl;
}

/**
* @brief Release previously acquired resources.
*
* @todo "Rule of five"
*/
template <class T>
SphereTree<T>::~SphereTree()
{
    // Free host memory
	delete[] m_nodes;
    delete[] m_leaves;
}

// ============================================================================
// ==== Import ================================================================
// ============================================================================

/**
* @brief Read & parse a .tree-file in the legacy format.
*
* Go through every line and parse the line depending on with what token it
* starts out with. NrOfSpheres is the only keyword in legacy format, so that
* is handled in a special case. Other lines start with 'H' for inner nodes
* ("Hierarchy") or 'L' for leaf nodes. Otherwise these two types of lines are
* processed very similarly.
*
* @param fileName The path to the .tree-file to be read & parsed.
*
* @warning __PD__CONVERT_IST_HIERARCHIES defined to 1 is needed to have the
*          complete hierarchy reconstructed.
*/
template <class T>
void SphereTree<T>::loadHierarchyLegacy(const char *fileName, const Matrix4<T> transform, T density)
{
    char lines[2048];
	std::ifstream istr_in;
	istr_in.open(fileName);

	if (!istr_in.is_open())
	{
		printf("%s can't be opened.\n", fileName);
        throw( "Error" );
	}

	int node_id = -1;
	int parent_id = -1;
	int child_id = -1;

	SphereNode<T> *tmp_nodes;
	int count_spheres;
	T x, y, z, r;

	int i = 0;
	while (istr_in && (!istr_in.eof()) && istr_in.getline(lines, 2048))
	{
//	    print(i++);
		if (strncmp(lines, "NrOfSpheres ", 12) == 0)
		{
			// Allocate memory for all spheres
			sscanf(lines, "NrOfSpheres %d", &count_spheres);
//			__PD__CUDA_CHECK(cudaMalloc(&m_dev_nodes, count_spheres * sizeof(SphereNode<T>)));
			m_nodes = new SphereNode<T>[count_spheres];
			tmp_nodes = new SphereNode<T>[count_spheres];
			m_countSpheres = count_spheres;

			printf("[pd] reading IST: %s size: %i (device)\n", fileName, count_spheres);
		}
		else if (lines[0] == 'L' || lines[0] == 'H')
		{
			// Reading new sphere entry
			if (std::is_same<T, float>::value)
				sscanf(lines, "%*[LH] %d %d %d %f %f %f %f %*f %*f %*f %*f %*f %*f", &node_id, &parent_id, &child_id, &x, &y, &z, &r );
			else if (std::is_same<T, double>::value)
				sscanf(lines, "%*[LH] %d %d %d %lf %lf %lf %lf %*f %*f %*f %*f %*f %*f", &node_id, &parent_id, &child_id, &x, &y, &z, &r );

			m_nodes[node_id].m_nodeId = node_id;

#if __PD__CONVERT_IST_HIERARCHIES
			m_nodes[node_id].m_parentId = parent_id;
			m_nodes[node_id].m_childId = child_id;
#endif
//			print(x,y,z);
			m_nodes[node_id].m_center.x = x;
			m_nodes[node_id].m_center.y = y;
			m_nodes[node_id].m_center.z = z;
			transform.multFullMatrixPnt( m_nodes[node_id].m_center, m_nodes[node_id].m_center );
			m_nodes[node_id].m_radius = r * transform.extractScale().x;
			x = m_nodes[node_id].m_center.x;
			y = m_nodes[node_id].m_center.y;
			z = m_nodes[node_id].m_center.z;
			r = m_nodes[node_id].m_radius;
            m_nodes[node_id].m_density = density;
			tmp_nodes[node_id].m_nodeId = node_id;
#if __PD__CONVERT_IST_HIERARCHIES
			tmp_nodes[node_id].m_parentId = parent_id;
			tmp_nodes[node_id].m_childId = child_id;
#endif
			tmp_nodes[node_id].m_center = m_nodes[node_id].m_center;
			tmp_nodes[node_id].m_radius = r;
            tmp_nodes[node_id].m_density = density;

			if (lines[0] == 'H')
			{
				m_countNodes++;
				m_nodes[node_id].m_isLeaf = false;
				tmp_nodes[node_id].m_isLeaf = false;
				m_hasHierarchy = true;
			}
			else if (lines[0] == 'L')
			{
				m_countLeaves++;
				m_nodes[node_id].m_isLeaf = true;
				tmp_nodes[node_id].m_isLeaf = true;
                m_aabbLo.x = fminf( m_aabbLo.x, fminf(x-abs(r),x+abs(r)) );
                m_aabbLo.y = fminf( m_aabbLo.y, fminf(y-abs(r),y+abs(r)) );
                m_aabbLo.z = fminf( m_aabbLo.z, fminf(z-abs(r),z+abs(r)) );
                m_aabbHi.x = fmaxf( m_aabbHi.x, fmaxf(x-abs(r),x+abs(r)) );
                m_aabbHi.y = fmaxf( m_aabbHi.y, fmaxf(y-abs(r),y+abs(r)) );
                m_aabbHi.z = fmaxf( m_aabbHi.z, fmaxf(z-abs(r),z+abs(r)) );
			}

			// Insert pointer to the newly created sphere
			if (node_id == 0)
			{
				m_rootNode = &m_nodes[node_id];
				m_dev_rootNode = &m_dev_nodes[node_id];
			}
			else
			{
				m_nodes[parent_id].m_children[child_id] = &m_nodes[node_id];
				tmp_nodes[parent_id].m_children[child_id] = &m_dev_nodes[node_id];
			}
		}
	}

//	__PD__CUDA_CHECK(cudaMemcpy(m_dev_nodes, tmp_nodes, count_spheres * sizeof(SphereNode<T>), __PD__H2D));
//	__PD__CUDA_CHECK(cudaDeviceSynchronize());

	delete(tmp_nodes);
}

/**
* @brief Parse the given .tree-file and populate this SphereTree instance.
*
* Go through each line of the given .tree-file and parse it.
* To parse a single line, we utilize boost::tokenizer, which will simply find
* all individual tokens, that are separated by a delimiter in the current
* line. Then we handle each line depending on the first token, if there is a
* keyword, we handle it accordingly, otherwise we assume it is a new sphere.
* Those lines should not start with a '#' and all necessary keywords should
* already have been read.
* If all those criteria are met we check whether it is a leaf or an inner node
* and process accordingly. Leaves have neighbors, which are appended at the
* end of the line and need to be parsed additionally.
* We hold an additional copy of the IST hierarchy here (tmp_nodes) so that we
* can fill that up with GPU pointers for m_children pointing to GPU memory.
* If the file is parsed we copy the CPU data to the GPU.
*
* @param fileName The path to the .tree-file to be read & parsed.
*/
template<class T>
template<int N>
void SphereTree<T>::loadHierarchy(const char *fileName, const Matrix4<T> transform, T density)
{
    char buff[FILENAME_MAX];
    GetCurrentDir( buff, FILENAME_MAX );
    std::string current_working_dir(buff);
    std::cout << "[SphereTree] reading " << fileName << " from " << current_working_dir << std::endl;
	std::ifstream in(fileName);
	std::string line;

	// This copy will be used to temporarily hold GPU memory pointers
	SphereNode<T> *tmp_nodes;

	while (std::getline(in, line))
	{


		// Read current line
		std::istringstream iss(line);
		std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss},
			std::istream_iterator<std::string>{} };
		// --- Processing this line depends on first token --------------------

		if (tokens[0] == "NodeCount")
		{
			m_countNodes = std::stoul(tokens[1]);
            std::cout << "[SphereTree] node count " << m_countNodes << std::endl;
			m_nodes = new SphereNode<T>[m_countNodes];
			tmp_nodes = new SphereNode<T>[m_countNodes];
//			__PD__CUDA_CHECK(cudaMalloc(&m_dev_nodes, m_countNodes * sizeof(SphereNode<T>)));
		}
		else if (tokens[0] == "LeafCount")
		{
			m_countLeaves = std::stoul(tokens[1]);
            std::cout << "[SphereTree] leaf count " << m_countLeaves << std::endl;
			m_leaves = new SphereNode<T>[m_countLeaves];
//			__PD__CUDA_CHECK(cudaMalloc(&m_dev_leaves, m_countLeaves * sizeof(SphereNode<T>)));
		}
        else if (tokens[0] == "NeighbourArraySize")
            continue;
		else if (tokens[0] == "TotalLeafVolume")
		{
			m_volume = std::stod(tokens[1]);
		}
		else if (tokens[0] == "CenterOfMass")
		{
			m_centerOfMass.x = std::stof(tokens[1]);
			m_centerOfMass.y = std::stof(tokens[2]);
			m_centerOfMass.z = std::stof(tokens[3]);
		}
		else if (tokens[0] == "InertiaTensor")
		{
			for (int j = 0; j<4; j++)
				for (int i = 0; i<4; i++)
					m_inertiaTensor(j, i) = std::stof(tokens[j * 4 + i + 1]);
			m_invInertiaTensor = m_inertiaTensor;
			m_invInertiaTensor.invert();
		}
		else if (m_nodes != nullptr && m_leaves != nullptr && (tokens[0][0] != '#'))
		{
			// Reading new sphere entry
//			bool leaf = std::stoi(tokens[0]) != 0;
//			unsigned int node_id = std::stoul(tokens[1]);
//			unsigned int parent_id = std::stoul(tokens[2]);
//			unsigned int child_id = std::stoul(tokens[3]);
//			Vec3<T> center(std::stof(tokens[4]), std::stof(tokens[5]), std::stof(tokens[6]));
//            transform.multFullMatrixPnt(center, center);
//			T radius = std::stof(tokens[7]);
//			radius *= transform.extractScale().x;

            bool leaf = std::stoi(tokens[0]) != 0;
            Vec3<T> center(-1.f*std::stof(tokens[2]), -1.f*std::stof(tokens[3]), -1.f*std::stof(tokens[4]));
            transform.multFullMatrixPnt(center, center);
            unsigned int node_id = std::stoul(tokens[1]);
            T radius = std::stof(tokens[5]);



			if (leaf)
			{
				// Leaf
				m_leaves[node_id].m_nodeId = node_id;

				m_leaves[node_id].m_center = center;
				m_leaves[node_id].m_radius = radius;
                m_leaves[node_id].m_density = density;
				m_leaves[node_id].m_isLeaf = leaf;

                m_aabbLo.x = fminf( m_aabbLo.x, center.x-abs(radius) );
                m_aabbLo.y = fminf( m_aabbLo.y, center.y-abs(radius) );
                m_aabbLo.z = fminf( m_aabbLo.z, center.z-abs(radius) );
                m_aabbHi.x = fmaxf( m_aabbHi.x, center.x+abs(radius) );
                m_aabbHi.y = fmaxf( m_aabbHi.y, center.y+abs(radius) );
                m_aabbHi.z = fmaxf( m_aabbHi.z, center.z+abs(radius) );
			}
			else
			{
				// Inner Node
				m_nodes[node_id].m_nodeId = node_id;
				m_nodes[node_id].m_center = center;
				m_nodes[node_id].m_radius = radius;
				m_nodes[node_id].m_isLeaf = leaf;
				tmp_nodes[node_id] = m_nodes[node_id];
				m_hasHierarchy = true;
			}

			// Insert pointer to the newly created sphere in his parent
			if (leaf)
			{
//				m_nodes[parent_id].m_children[child_id] = &m_leaves[node_id];
//				tmp_nodes[parent_id].m_children[child_id] = &m_dev_leaves[node_id];
			}
            else
            {
                if (node_id == 0)
                {
                    m_rootNode = &m_nodes[node_id];
                    m_dev_rootNode = &m_dev_nodes[node_id];
                }
                else
                {
//                    m_nodes[parent_id].m_children[child_id] = &m_nodes[node_id];
//                    tmp_nodes[parent_id].m_children[child_id] = &m_dev_nodes[node_id];
                }
            }
		}
	}

//	m_countSpheres = m_countNodes + m_countLeaves;

	delete(tmp_nodes);
}

/**
 * @brief TODO
 */
template <class T>
void SphereTree<T>::exportCcsObj( const char *objFileName, const char *fileName )
{
    // ---- Read triangle mesh ----
    std::vector<SphereNode<T>> inside_leaves;
    std::vector<Vec3<T>> vertices;
    std::vector<std::vector<unsigned int>> indices;
    std::ifstream obj;
    std::string line;
    obj.open( objFileName, std::ios::in );

    while (std::getline(obj, line))
    {
        // Read current line
		std::istringstream iss( line );
		std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss},
			                             std::istream_iterator<std::string>{} };
        if (tokens[0] == "v") {
            vertices.push_back( pd::Vec3<T>(std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])) );
        } else if (tokens[0] == "f") {
            std::vector<unsigned int> face;
            face.push_back( std::stoul(tokens[1])-1 );
            face.push_back( std::stoul(tokens[2])-1 );
            face.push_back( std::stoul(tokens[3])-1 );
            indices.push_back( face );
        }
    }

    // Check which spheres are inside
    for (unsigned int i=0; i<m_countLeaves; i++)
    {
        if (i%(m_countLeaves/100)==0)
            std::cout << "[pd::SphereTree] inside check working " << roundf((T)i/m_countLeaves*100.f) << "%" << std::endl;
        pd::SphereNode<T> *leaf = &m_leaves[i];
        pd::Vec3<T> line( 1.f, 0.2f, -0.2f );
        line.normalize();
        Vec3<T> O = leaf->m_center;
        int intersections = 0;
        for (auto f : indices) {
            if (rayTriangleTest(O,line,0,vertices[f[0]],vertices[f[1]],vertices[f[2]]))
                intersections++;
        }

        if (intersections > 0 && intersections % 2 == 1)
            inside_leaves.push_back( *leaf );
    }

    const std::vector<Vec3<T>> sphere_vertices {
        Vec3<T>(-0.309016973, -0.9510566, 8.3144009e-08),
        Vec3<T>(-0.809016824, -0.587785423, 5.13858325e-08),
        Vec3<T>(-0.525731206, -0.723606944, -0.447213203),
        Vec3<T>(-1, 0, 0),
        Vec3<T>(-0.809016824, 0.587785423, -5.13858325e-08),
        Vec3<T>(-0.850650847, 0.276393622, -0.447213203),
        Vec3<T>(-0.309016973, 0.9510566, -8.3144009e-08),
        Vec3<T>(0.309016973, 0.9510566, -8.3144009e-08),
        Vec3<T>(0, 0.894427359, -0.447213203),
        Vec3<T>(0.809016824, 0.587785423, -5.13858325e-08),
        Vec3<T>(1, 0, 0),
        Vec3<T>(0.850650847, 0.276393622, -0.447213203),
        Vec3<T>(0.809016824, -0.587785423, 5.13858325e-08),
        Vec3<T>(0.309016973, -0.9510566, 8.3144009e-08),
        Vec3<T>(0.525731206, -0.723606944, -0.447213203),
        Vec3<T>(0, -0.850651026, -0.52573061),
        Vec3<T>(-0.809017241, -0.262865454, -0.525730789),
        Vec3<T>(-0.49999997, 0.688191354, -0.525730729),
        Vec3<T>(0.49999997, 0.688191354, -0.525730729),
        Vec3<T>(0.809017241, -0.262865454, -0.525730789),
        Vec3<T>(0, -8.74227766e-08, -1),
        Vec3<T>(0.309017003, -0.425325483, -0.850650847),
        Vec3<T>(-0.309017003, -0.425325483, -0.850650847),
        Vec3<T>(0.49999997, 0.162460044, -0.850650787),
        Vec3<T>(0, 0.525731087, -0.850650787),
        Vec3<T>(-0.49999997, 0.162460044, -0.850650787)
    };

    const std::vector<Vec3<T>> sphere_vnormals {
        Vec3<T>(-0.294135749, -0.929526031, -0.22240822),
        Vec3<T>(-0.793138981, -0.566979051, -0.222408101),
        Vec3<T>(-0.525731146, -0.723606884, -0.447213322),
        Vec3<T>(-0.974924862, -0.00749956351, -0.222407967),
        Vec3<T>(-0.7843225, 0.57911396, -0.222408041),
        Vec3<T>(-0.850650847, 0.276393473, -0.447213322),
        Vec3<T>(-0.308400929, 0.924891055, -0.222408265),
        Vec3<T>(0.308400899, 0.924891114, -0.22240828),
        Vec3<T>(-1.57931215e-08, 0.894427299, -0.447213531),
        Vec3<T>(0.784322441, 0.57911396, -0.222408041),
        Vec3<T>(0.974924862, -0.00749956351, -0.222407967),
        Vec3<T>(0.850650907, 0.276393473, -0.447213322),
        Vec3<T>(0.793138981, -0.566979051, -0.222408086),
        Vec3<T>(0.294135749, -0.929526031, -0.22240822),
        Vec3<T>(0.525731146, -0.723606884, -0.447213382),
        Vec3<T>(-2.11272813e-08, -0.850650966, -0.525730848),
        Vec3<T>(-0.809017241, -0.262865484, -0.525730908),
        Vec3<T>(-0.49999997, 0.688191116, -0.525730968),
        Vec3<T>(0.49999997, 0.688191175, -0.525730968),
        Vec3<T>(0.809017241, -0.262865484, -0.525730908),
        Vec3<T>(0, 4.31507772e-08, -1),
        Vec3<T>(0.309017062, -0.425325513, -0.850650728),
        Vec3<T>(-0.309017062, -0.425325513, -0.850650728),
        Vec3<T>(0.50000006, 0.162460059, -0.850650728),
        Vec3<T>(2.1127283e-08, 0.525731266, -0.850650668),
        Vec3<T>(-0.500000119, 0.162460059, -0.850650728),
    };

    int w = sphere_vertices.size();
    const int sphere_indices[] = {
        -(w+1-3),  -(w+1-1),  -(w+1-2),
        -(w+1-6),  -(w+1-4),  -(w+1-5),
        -(w+1-9),  -(w+1-7),  -(w+1-8),
        -(w+1-12), -(w+1-10), -(w+1-11),
        -(w+1-15), -(w+1-13), -(w+1-14),
        -(w+1-16), -(w+1-1),  -(w+1-3),
        -(w+1-15), -(w+1-14), -(w+1-16),
        -(w+1-14), -(w+1-1),  -(w+1-16),
        -(w+1-17), -(w+1-4),  -(w+1-6),
        -(w+1-3),  -(w+1-2),  -(w+1-17),
        -(w+1-2),  -(w+1-4),  -(w+1-17),
        -(w+1-18), -(w+1-7),  -(w+1-9),
        -(w+1-6),  -(w+1-5),  -(w+1-18),
        -(w+1-5),  -(w+1-7),  -(w+1-18),
        -(w+1-19), -(w+1-10), -(w+1-12),
        -(w+1-9),  -(w+1-8),  -(w+1-19),
        -(w+1-8),  -(w+1-10), -(w+1-19),
        -(w+1-20), -(w+1-13), -(w+1-15),
        -(w+1-12), -(w+1-11), -(w+1-20),
        -(w+1-11), -(w+1-13), -(w+1-20),
        -(w+1-23), -(w+1-21), -(w+1-22),
        -(w+1-16), -(w+1-22), -(w+1-15),
        -(w+1-3),  -(w+1-23), -(w+1-16),
        -(w+1-23), -(w+1-22), -(w+1-16),
        -(w+1-22), -(w+1-21), -(w+1-24),
        -(w+1-20), -(w+1-24), -(w+1-12),
        -(w+1-15), -(w+1-22), -(w+1-20),
        -(w+1-22), -(w+1-24), -(w+1-20),
        -(w+1-24), -(w+1-21), -(w+1-25),
        -(w+1-19), -(w+1-25), -(w+1-9),
        -(w+1-12), -(w+1-24), -(w+1-19),
        -(w+1-24), -(w+1-25), -(w+1-19),
        -(w+1-25), -(w+1-21), -(w+1-26),
        -(w+1-18), -(w+1-26), -(w+1-6),
        -(w+1-9),  -(w+1-25), -(w+1-18),
        -(w+1-25), -(w+1-26), -(w+1-18),
        -(w+1-26), -(w+1-21), -(w+1-23),
        -(w+1-17), -(w+1-23), -(w+1-3),
        -(w+1-6),  -(w+1-26), -(w+1-17),
        -(w+1-26), -(w+1-23), -(w+1-17),
    };

    std::ofstream file;
    file.open( fileName, std::ios::out | std::ios::trunc );
    std::ofstream file_sphere;
    std::string inside_file( objFileName );
    inside_file += "." + std::to_string( inside_leaves.size() );
    inside_file += ".spheres";
    file_sphere.open( inside_file, std::ios::out | std::ios::trunc );
    if (file.is_open() && file_sphere.is_open())
    {
        printf( "[pd] Exporting new sphere packing (inside spheres) to %s.\n", inside_file.c_str() );
        printf( "[pd] Exporting CCs as obj to %s.\n", fileName );

        file << "g CC" << std::setfill('0') << std::setw(10) << 0 << std::endl;
        for (int i=0; i<inside_leaves.size(); i++) {
            auto S = inside_leaves[i].m_center;
            auto r = inside_leaves[i].m_radius;
            file_sphere << S.x << " " << S.y << " " << S.z << " " << r;
            for (int j=0; j<30; j++)
                file_sphere << " 0";
            file_sphere << std::endl;

            for (auto v : sphere_vertices) {
                v *= r;
                v += S;
                file << "v " << v.x << " " << v.y << " " << v.z << std::endl;
            }
            for (auto vn : sphere_vnormals)
                file << "vn " << vn.x << " " << vn.y << " " << vn.z << std::endl;
            for (int k=0; k<(120/3-1); k++) {
                auto *f = &sphere_indices[3*k];
                file << "f " << f[0]<<"/"<<f[0]<<"/"<<f[0] << " " << f[1]<<"/"<<f[1]<<"/"<<f[1] << " " << f[2]<<"/"<<f[2]<<"/"<<f[2] << std::endl;
            }
            file << "" << std::endl;
        }
    }
    else
        printf( "[pd] ERROR: Couldn't export to target path!\n" );

}


/**
 * @brief TODO
 */
template <class T>
bool SphereTree<T>::rayTriangleTest(Vec3<T> O, Vec3<T> d, T dMax, Vec3<T> v0, Vec3<T> v1, Vec3<T> v2) const
{
    // Moller Trumbore method
    Vec3<T> v0v1 = v1 - v0;
    Vec3<T> v0v2 = v2 - v0;
    Vec3<T> pvec = d.crossProduct(v0v2);
    T det = v0v1.dot(pvec);
    // ray and triangle are parallel if det is close to 0
    if (fabs(det) < 0.00001f) return false;
    T invDet = 1 / det;

    Vec3<T> tvec = O - v0;
    auto u = tvec.dot(pvec) * invDet;
    if (u < 0 || u > 1) return false;

    Vec3<T> qvec = tvec.crossProduct(v0v1);
    auto v = d.dot(qvec) * invDet;
    if (v < 0 || u + v > 1) return false;

    auto t = v0v2.dot(qvec) * invDet;

    if (t > 0 && (dMax == 0 || t < dMax))
        return true;

    return false;
}


/**
 * @brief TODO
 */
template <class T>
void SphereTree<T>::accumulateNeighboursDown( SphereNode<T> *leaf, SphereNode<T> *root, std::vector<std::vector<unsigned int>> *neighbours2d ) const
{
	static const T epsilon = 1.0;
    for (int i=0; i<4; i++) {
        auto child = root->m_children[i];
        if (child == nullptr)
            continue;
        if (leaf->distance(*child) < epsilon) {
            if (child->m_isLeaf) {
                if (leaf->m_nodeId != child->m_nodeId) {
                    (*neighbours2d)[leaf->m_nodeId].push_back( child->m_nodeId );
                }
            } else {
                accumulateNeighboursDown( leaf, child, neighbours2d );
            }
        } else {
        }
    }
}

/**
 * @brief Convert a legacy SphereTree to a modern one.
 *
 * This takes the legacy SphereTree and converts it to the most recent version
 * which is used in this application.
 * - The biggest difference between them is
 *   that leaves are stored separately from inner nodes, but they are still
 *   part of the hierarchy because there are cross-references between hierarchy
 *   and leaf-memory.
 * - Another difference is that legacy SphereTree's have no ISTGraph, and so
 *   SphereNode's don't have neighbour data, which is also generated here.
 * - Also this will calculate some physics characteristics of the IST for the
 *   torque simulation.
 *
 * @warning This SphereTree needs to be a legacy one for this to make sense.
 *
 * @todo Test if this really works without calling export & load afterwards,
 *       this mostly depends on if the GPU memory is correctly filled.
 * @todo This function could probably be factorized in auxiliary functions for
 *       better overview.
 * @todo Maybe introduce m_legacy bool member to check if this call is legal.
 */
template <class T>
void SphereTree<T>::convertToNewVersion()
{
#if __PD__CONVERT_IST_HIERARCHIES
	// ---- Separate nodes & leaves -------------------------------------------

    // Fix node-less ISTs
    m_countNodes = std::max( m_countNodes, static_cast<unsigned int>(1) );

	// - CPU storage
	std::map<int, std::pair<bool, int>> map;
	SphereNode<T> *inodes = new SphereNode<T>[m_countNodes];
	SphereNode<T> *tmp_inodes = new SphereNode<T>[m_countNodes];
	SphereNode<T> *leaves = new SphereNode<T>[m_countLeaves];

	// - GPU storage
	__PD__CUDA_CHECK(cudaFree(m_dev_nodes));
	__PD__CUDA_CHECK(cudaMalloc(&m_dev_nodes, m_countNodes * sizeof(SphereNode<T>)));
	__PD__CUDA_CHECK(cudaMalloc(&m_dev_leaves, m_countLeaves * sizeof(SphereNode<T>)));
	unsigned int leaf_counter = 0;
	unsigned int node_counter = 0;

	// Convert to temporary data structure that maps old-id ==> (leaf?, new-id).
	// This is so we have new nodeIds that are separate for leaves & nodes.
	for (unsigned int i = 0; i<m_countSpheres; i++)
	{
		pd::SphereNode<T> *old_node(&m_nodes[i]);

		if (old_node->m_isLeaf)
			map[old_node->m_nodeId] = std::pair<bool, int>(true, leaf_counter++);
		else
			map[old_node->m_nodeId] = std::pair<bool, int>(false, node_counter++);
	}

	// Convert to new m_nodes & m_leaves (temporarily in inodes & leaves)
	// Analog to convert SphereTree::loadHierarchy(const char*), we hold an
	// additional temporary array for GPU pointers which we will later copy.
	for (unsigned int i = 0; i<m_countSpheres; i++)
	{
		pd::SphereNode<T> *old_node(&m_nodes[i]);

		if (old_node->m_nodeId == 0)
			m_rootNode = &inodes[node_counter];

		int id = map[old_node->m_nodeId].second;
		int parent_id = map[old_node->m_parentId].second;

		if (old_node->m_isLeaf)
		{
			inodes[parent_id].m_children[old_node->m_childId] = &leaves[id];
			leaves[id].m_nodeId = id;
			leaves[id].m_center = old_node->m_center;
			leaves[id].m_radius = old_node->m_radius;
			leaves[id].m_isLeaf = true;
			leaves[id].m_childId = old_node->m_childId;
			leaves[id].m_parentId = parent_id;
			leaves[id].m_density = old_node->m_density;

            // Fix node-less ISTs
            if (node_counter == 0)
            {
                inodes[node_counter] = leaves[id];
                inodes[node_counter].m_isLeaf = false;
                inodes[node_counter].m_children[0] = &leaves[id];
                node_counter++;
            }
		}
		else
		{
			inodes[parent_id].m_children[old_node->m_childId] = &inodes[id];
			inodes[id].m_nodeId = id;
			inodes[id].m_center = old_node->m_center;
			inodes[id].m_radius = old_node->m_radius;
			inodes[id].m_isLeaf = false;
			inodes[id].m_childId = old_node->m_childId;
			inodes[id].m_parentId = parent_id;
			tmp_inodes[id] = inodes[id];
			tmp_inodes[parent_id].m_children[old_node->m_childId] = &m_dev_nodes[id];
		}
	}

	// Overwrite attributes
	delete m_nodes;
	m_nodes = inodes;
	m_leaves = leaves;

	// Copy to GPU
	__PD__CUDA_CHECK(cudaMemcpy(m_dev_nodes, tmp_inodes, m_countNodes * sizeof(SphereNode<T>), __PD__H2D));
	__PD__CUDA_CHECK(cudaMemcpy(m_dev_leaves, m_leaves, m_countLeaves * sizeof(SphereNode<T>), __PD__H2D));
	__PD__CUDA_CHECK(cudaDeviceSynchronize());
	delete tmp_inodes;

	// ---- Physics Calculations ----------------------------------------------

	// Calculate combined inertia tensor, volume & center of mass
	m_inertiaTensor = Matrix4<T>::ZERO();
	for (unsigned int i = 0; i<m_countLeaves; i++)
	{
		pd::SphereNode<T> *leaf(&m_leaves[i]);

		// Calculate total volume to scale force with later
		T volume = leaf->calcFullVolume();
		m_volume += volume;
		m_centerOfMass += leaf->m_center * volume;

		// Calculate local tensor for the current sphere
		T m(leaf->calcFullVolume() * Density);
		T r(leaf->m_radius);
		T inertia(2.f / 5.f * m * r*r);
		Matrix4<T> inertia_tensor_local(0.f);
		inertia_tensor_local(0, 0) = inertia;
		inertia_tensor_local(1, 1) = inertia;
		inertia_tensor_local(2, 2) = inertia;

		// Parallel axis theorem to have tensor relative to the center of mass of the whole IST
		// J_ij = I_ij + m c, with c = R^2 d_ij - R_i R_j and d_ij = 1 iff i=j else 0
		// https://en.wikipedia.org/wiki/Parallel_axis_theorem#Tensor_generalization [Date: 2016.05.23]
		Vec3<T> R(m_centerOfMass - leaf->m_center);
		Matrix4<T> c(R.y*R.y + R.z*R.z, -R.x*R.y, -R.x*R.z, 0.0,
			-R.x*R.y, R.x*R.x + R.z*R.z, -R.y*R.z, 0.0,
			-R.x*R.z, -R.y*R.z, R.x*R.x + R.y*R.y, 0.0,
			0.0, 0.0, 0.0, 0.0);
		Matrix4<T> inertia_tensor_com(inertia_tensor_local + (m*c));

		// Add to total tensor
		m_inertiaTensor = m_inertiaTensor + inertia_tensor_com;
	}

	std::cout << "[pd::SphereTree] Packing volume: " << m_volume << std::endl;
	std::cout << "[pd::SphereTree] Center of mass: " << m_centerOfMass << std::endl;

	// Normalize center of mass by total volume
	m_centerOfMass /= m_volume;

	// Inverted inertia tensor
	m_invInertiaTensor = m_inertiaTensor;
	m_invInertiaTensor.invert();

    // Statistics
#endif
}

/**
 * TODO
 */
template <class T>
void SphereTree<T>::generateRegions( Vec3<T> center, Matrix4<T> transform )
{
#if __IHC_MAHIDOL
    Vec3<T> normal_north_west( 0, 0, 1 );
    Vec3<T> normal_north_east( 0, 0, 1 );
    Matrix4<T> M_west; M_west.identity();
    Matrix4<T> M_east; M_east.identity();
    T angle_west = (45.0/180.0) * __PD__PI;
    T angle_east = (315.0/180.0) * __PD__PI;
    M_west(0,0) =  cos(angle_west); M_east(0,0) =  cos(angle_east);
    M_west(0,2) =  sin(angle_west); M_east(0,2) =  sin(angle_east);
    M_west(2,0) = -sin(angle_west); M_east(2,0) = -sin(angle_east);
    M_west(2,2) =  cos(angle_west); M_east(2,2) =  cos(angle_east);
    M_west = transform * M_west;
    M_east = transform * M_east;
    M_west.multRotMatrixPnt( normal_north_west, normal_north_west );
    M_east.multRotMatrixPnt( normal_north_east, normal_north_east );
    T wNx = normal_north_west.x;
    T wNy = normal_north_west.y;
    T wNz = normal_north_west.z;
    T eNx = normal_north_east.x;
    T eNy = normal_north_east.y;
    T eNz = normal_north_east.z;
    T wNd = center.dot( normal_north_west );
    T eNd = center.dot( normal_north_east );

    for (int i=0; i<m_countLeaves; i++)
    {
        SphereNode<T> *leaf = &m_leaves[i];
        T x = leaf->m_center.x;
        T y = leaf->m_center.y;
        T z = leaf->m_center.z;
        bool north_west = (wNx*x + wNy*y + wNz*z + wNd) < 0;
        bool north_east = (eNx*x + eNy*y + eNz*z + eNd) < 0;
        int region = 0;
        if (north_west && north_east)
            region = 0;
        else if (!north_west && north_east)
            region = 1;
        else if (!north_west && !north_east)
            region = 2;
        else if (north_west && !north_east)
            region = 3;
        leaf->m_region = region;
        leaf->m_score = 0;
        leaf->m_weight = 0;
    }
#endif
}

// ============================================================================
// ==== Export ================================================================
// ============================================================================


/**
* @brief Write the current SphereTree to a new .tree-file.
*
* @param fileName Path to the file that is to be written to.
*
* @warning This exports old version type hierarchies.
*/
template <class T>
void SphereTree<T>::exportHierarchyLegacy(const char *fileName)
{
	std::ofstream file;
	file.open(fileName, std::ios::out | std::ios::trunc);

	if (file.is_open())
	{
		printf("[pd] Exporting hierarchy (legacy format) to %s.\n", fileName);

		// Comment that describes the layout of the file
		file << "#LeafOrHierarchy [L or H] NodeID ParentID ChildNr nodeCenter[x, y, z] nodeRadius VolumeCardinal[x, y, z, -x, -y, -z] nodeVolumeInner nodeVolumeRatio normalDirection[x, y, z] normalAngle, IfLeaf: closestTriangleP1[x, y, z] closestTriangleP2[x, y, z] closestTriangleP3[x, y, z]\n";

		// Print keywords
		file << "NrOfSpheres " << m_countSpheres << "\n";

        // Traverse spheres in breads-first-style
        int node_counter = 0;
        std::queue<SphereNode<T>*> queue;
        queue.push( m_rootNode );
        SphereNode<T> *current;
        while (!queue.empty())
        {
            // Get first element of queue
            current = queue.front();
            queue.pop();

            // Push children
            for (int i=0; i<4 && current->m_children[i]; i++)
                queue.push( current->m_children[i] );

		    // Print sphere
			SphereNode<T> *node = current;
			file << (node->m_isLeaf ? "L" : "H") << " ";
			file << node_counter++ << " ";
			file << node->m_parentId << " " << node->m_childId << " ";
			file << node->m_center.x << " " << node->m_center.y << " " << node->m_center.z << " ";

			file << node->m_radius << " ";


            for (int i=0; i<12; i++)
			    file << 0.f << " ";
            if (node->m_isLeaf)
                for (int i=0; i<9; i++)
                    file << 0.f << " ";
            file << "\n";
        }

		file.close();
	}
	else
		printf("[pd] ERROR: Couldn't export to target path!\n");
}

/**
* @brief Write the current SphereTree to a new .tree-file.
*
* First we write important keywords at the top of the file, then we go through
* every inner node and write an appropriate line, then all leaves and do the
* same. For leaves we need to additionally write the neighbors of the leaf at
* the end of the line, since the length is variable without a set maximum.
*
* @param fileName Path to the file that is to be written to.
*
* This exports new version type hierarchies.
* @warning Only works with new version SphereTree's, legacy ones are untested.
*          Use SphereTree::loadHierarchyLegacy(const char*) =>
*          SphereTree::convertToNewVersion() =>
*          SphereTree::exportHierarchy(const char*) to convert your legacy
*          type hierarchy files to new version & export to a file.
*          To simply load a new file exported with this just load via
*          SphereTree::loadHierarchy(const char*).
*/
template <class T>
void SphereTree<T>::exportHierarchy(const char *fileName)
{

	std::ofstream file;
	file.open(fileName, std::ios::out | std::ios::trunc);

	if (file.is_open())
	{
		printf("[pd] Exporting hierarchy (new format) to %s.\n", fileName);

		// Comment that describes the layout of the file
		file << "# isLeaf? {Node,Leaf}ID ParentID ChildID Center[x y z] Radius";
#if __IHC_MAHIDOL
        file << " Leaf[Region Score Weight]";
#endif
		file << "\n";

		// Print keywords
		file << "NodeCount " << m_countNodes << "\n";
		file << "LeafCount " << m_countLeaves << "\n";
		file << "TotalLeafVolume " << m_volume << "\n";
		file << "CenterOfMass " << m_centerOfMass.x << " " << m_centerOfMass.y << " " << m_centerOfMass.z << "\n";
		file << "InertiaTensor";
		for (int i = 0; i<4; i++)
			for (int j = 0; j<4; j++)
				file << " " << m_inertiaTensor(i, j);
		file << "\n";

		// Print all nodes
		for (unsigned int i = 0; i<m_countNodes; i++)
		{
			const SphereNode<T> &node = m_nodes[i];
			file << node.m_isLeaf << " ";
			file << node.m_nodeId << " ";
			file << node.m_parentId << " " << node.m_childId << " ";
			file << node.m_center.x << " " << node.m_center.y << " " << node.m_center.z << " ";
			file << node.m_radius << " ";
			file << "\n";
		}

		// Lastly, print all leaves
		for (unsigned int i = 0; i<m_countLeaves; i++)
		{
		    print(i);
			const SphereNode<T> &leaf = m_leaves[i];
			file << leaf.m_isLeaf << " ";
			file << leaf.m_nodeId << " ";
			file << leaf.m_parentId << " " << leaf.m_childId << " ";
			file << leaf.m_center.x << " " << leaf.m_center.y << " " << leaf.m_center.z << " ";
			file << leaf.m_radius << " ";
			file << "\n";
		}

		file.close();
	}
	else
		printf("[pd] ERROR: Couldn't export to target path!\n");
}

// ============================================================================
// ==== Copy ==================================================================
// ============================================================================

///**
//* @brief Create a copy of this SphereTree and safe it in ptr.
//*
//* In case ptr does not point to a memory address, we allocate new memory for
//* it. Next we delete old data in ptr, if there is any and allocate new memory
//* for the heap storage & copy data over.
//* The same is done for GPU attributes, but it is unclear if that is needed,
//* right now this is just to debug the GPU data of this by copying them to
//* ptr's CPU attributes to be able to read & visualize them.
//*
//* @param ptr A reference of a pointer, this will point to the new copy afterwards.
//* @param stream A CUDA stream on which to perform GPU data transfers.
//*
//* @warning Copies created with this function don't have their m_children
//*          members in m_nodes updated, so they point to the old instance.
//* @warning m_dev_nodes & m_dev_leaves have no meaningful data afterwards.
//* @warning The created copy is only useful for debugging data on CPU.
//* @warning In the current implementation, this is only useful when previously
//*          SphereTree::copyGPUtoCPU() was called, which will copy all
//*          SphereNode's from GPU memory to CPU memory, which are used here.
//* @todo Create complete & valid copies here? Depends on the use-case.
//*/
//template <class T>
//void SphereTree<T>::createCopy(SphereTree<T> *&ptr, cudaStream_t stream)
//{
//	// Allocate memory for SphereTree in case it is null
//	if (!ptr)
//	{
//		ptr = new SphereTree();
//		ptr->m_nodes = nullptr;
//		ptr->m_rootNode = nullptr;
//		ptr->m_dev_nodes = nullptr;
//		ptr->m_dev_rootNode = nullptr;
//		ptr->m_countNodes = m_countNodes;
//		ptr->m_countLeaves = m_countLeaves;
//		ptr->m_volume = m_volume;
//		ptr->m_centerOfMass = m_centerOfMass;
//	}
//
//	// ---- CPU Pointers ------------------------------------------------------
//
//	// Copy hierarchy nodes of IST
//	delete ptr->m_nodes;
//	ptr->m_nodes = new SphereNode<T>[ptr->m_countNodes];
//	memcpy(ptr->m_nodes, m_nodes, m_countNodes * sizeof(SphereNode<T>));
//	ptr->m_rootNode = &ptr->m_nodes[m_rootNode->m_nodeId];
//	// TODO FIXME: All SphereNode's in ptr->m_nodes have m_children pointers
//	//             pointing into m_nodes's m_children still. If this copy
//	//             is going to be used to traverse the hierarchy, fix this.
//
//	// Copy leaves of IST
//	delete ptr->m_leaves;
//	ptr->m_leaves = new SphereNode<T>[ptr->m_countLeaves];
//	memcpy(ptr->m_leaves, m_leaves, m_countLeaves * sizeof(SphereNode<T>));
//
//	// ---- GPU Pointers ------------------------------------------------------
//
//	// Copy hierarchy nodes of IST
//	__PD__CUDA_CHECK(cudaFree(ptr->m_dev_nodes));
//	__PD__CUDA_CHECK(cudaMalloc(&ptr->m_dev_nodes, ptr->m_countNodes * sizeof(SphereNode<T>)));
//	__PD__CUDA_CHECK(cudaMemcpy(ptr->m_dev_nodes, m_dev_nodes, ptr->m_countNodes * sizeof(SphereNode<T>), __PD__D2D));
//	ptr->m_dev_rootNode = &ptr->m_dev_nodes[m_rootNode->m_nodeId];
//	// TODO FIXME: The problem stated above is the case here too, m_dev_nodes
//	//             entires don't have useful m_children pointers.
//
//	// Copy leaves of IST
//	__PD__CUDA_CHECK(cudaFree(ptr->m_dev_leaves));
//	__PD__CUDA_CHECK(cudaMalloc(&ptr->m_dev_leaves, ptr->m_countLeaves * sizeof(SphereNode<T>)));
//	__PD__CUDA_CHECK(cudaMemcpy(ptr->m_dev_leaves, m_dev_leaves, ptr->m_countLeaves * sizeof(SphereNode<T>), __PD__D2D));
//}

///**
//* @brief Overwrite the CPU SphereNode's with those of the GPU.
//*
//* This is done, since when we run the GPU algorithms, additional data is
//* collected and accumulated per SphereNode, which we can then debug on the CPU.
//*
//* @param stream The CUDA stream on which to perform the GPU data transfers.
//*
//* @warning This will overwrite the CPU copy, meaning m_children (which used to
//*          be CPU pointers) will point to GPU memory, so they will be invalid.
//*/
//template <class T>
//void SphereTree<T>::copyGPUtoCPU(cudaStream_t stream)
//{
//	// TODO FIXME: Save m_nodes to fix children pointers after GPU to CPU copy
//	__PD__CUDA_CHECK(cudaMemcpyAsync(m_nodes, m_dev_nodes, m_countNodes * sizeof(SphereNode<T>), __PD__D2H, stream));
//	__PD__CUDA_CHECK(cudaMemcpyAsync(m_leaves, m_dev_leaves, m_countLeaves * sizeof(SphereNode<T>), __PD__D2H, stream));
//	__PD__CUDA_CHECK(cudaStreamSynchronize(stream));
//	// TODO: Maybe also copy neighbourRefs (if it is necessary in the future..)
//}

// ============================================================================
// ==== Queries ==============================================================
// ============================================================================

/**
* @brief Recursively find out if node child is a child of node with parentId.
*
* The check is performed by checking whether the current node child is the
* parent indicated by parentId, in which case we return true, otherwise we
* recursively call this function for the direct parent of child. In that
* fashion we traverse the hierarchy in reverse up until the root.
* We have some trivial cases to check to make the recursion stop. If parentId
* is the node-ID of the root, we just return true, since every node is a child
* of the root, unless the child is itself the root, in which case we return
* false. In that case we also have to check whether child is not a leaf,
* because we have separate node-IDs for nodes and leaves, which means there
* might be a leaf with the same node-ID as the root node, which should be OK
* to continue to traverse the hierarchy.
*
* @param child The node that is to be checked for being the child.
* @param parentId The node-ID of the node that is supposed to be the parent.
* @return Whether child is a child of node with parentId as id, in terms of
*         the given hierarchy.
* @retval true   Iff node with id parentId is the parent of the node child.
* @retval false  Otherwise.
*
* @warning Only works correctly when compiled with the following define set to
*          1: __PD__CONVERT_IST_HIERARCHIES.
*/
template <class T>
bool SphereTree<T>::isChildOf(SphereNode<T> *child, int parentId) const
{
	if (child->m_nodeId == m_rootNode->m_nodeId && !child->m_isLeaf)
		return false;
	if (parentId == m_rootNode->m_nodeId)
		return true;
	if (child->m_nodeId == parentId)
		return true;
#if __PD__CONVERT_IST_HIERARCHIES
	return isChildOf(&m_nodes[child->m_parentId], parentId);
#else
	return false;
#endif
}

// ============================================================================
// ==== ISTGraph ==============================================================
// ============================================================================

/**
* @brief Recursively find all leaves that are currently connected to leaf.
*
* We mark the current leaf as part of the connected component with id c and
* recursively call this function for all it's neighbors (based on the current
* iteration of the ISTGraph (described in neighbourIds)).
*
* @param c The current connected component id, to mark found leaves with.
* @param neighbourIds The preliminary ISTGraph iteration.
* @param id_to_cc Mapping leaf-IDs to connected components.
* @param leaf The current leaf from which to proceed recursion.
*/
template <class T>
void SphereTree<T>::depthFirstSearch(unsigned int c,
	std::vector<std::vector<unsigned int>> *neighbourIds,
	std::map<unsigned int, int> *id_to_cc,
	SphereNode<T> *leaf)
{
	(*id_to_cc)[leaf->m_nodeId] = c;
	for (unsigned int i = 0; i<(*neighbourIds)[leaf->m_nodeId].size(); i++)
	{
		pd::SphereNode<T> *neighbour = &m_leaves[(*neighbourIds)[leaf->m_nodeId][i]];

		if ((*id_to_cc)[neighbour->m_nodeId] == 0)
			depthFirstSearch(c, neighbourIds, id_to_cc, neighbour);
	}
}

/**
* @brief Find and return all connected components based on current ISTGraph.
*
* Initialize a map from leaf-ID to connected component ID and fill it by
* starting graph traversal from each leaf and combining all found leaves to
* the same connected component, that way all existing connected components are
* found. Afterwards we convert the results into a more useful data-structure
* than a map.
*
* @param neighbourIds The current (preliminary) iteration of the ISTGraph.
* @return A 2D-Array that puts all leaves that share a connected component in
*         to the same list.
*/
template <class T>
std::vector<std::vector<unsigned int>> SphereTree<T>::computeConnectedComponents(std::vector<std::vector<unsigned int>> neighbourIds)
{
	unsigned int c = 0;
	std::map<unsigned int, int> id_to_cc;

	// Init maps
	for (unsigned int i = 0; i<m_countLeaves; i++)
	{
		pd::SphereNode<T> *leaf = &m_leaves[i];
		id_to_cc[leaf->m_nodeId] = 0;
	}

	// Find connected components
	for (unsigned int i = 0; i<m_countLeaves; i++)
	{
		pd::SphereNode<T> *leaf = &m_leaves[i];
		if (id_to_cc[leaf->m_nodeId] == 0)
		{
			c++;
			depthFirstSearch(c, &neighbourIds, &id_to_cc, leaf);
		}
	}

	printf("[pd::SphereTree] found connected components: %u\n", c);

	// Create data-structure to represent the found connected components
	std::vector<std::vector<unsigned int>> result;
	for (unsigned int i = 0; i<c; i++)
	{
		std::vector<unsigned int> single_cc;
		for (unsigned int j = 0; j<m_countLeaves; j++)
		{
			pd::SphereNode<T> *leaf = &m_leaves[j];

			if (id_to_cc[leaf->m_nodeId] == i + 1)
				single_cc.push_back(leaf->m_nodeId);
		}

		if (!single_cc.empty())
			result.push_back(single_cc);
	}

	return result;
}
//
///**
//* @brief Insert bridge into the preliminary ISTGraph.
//*
//* This will simply go through every connected component and every leaf in
//* those and compare each to every other connected component and each leaf in
//* those. The pair of leaves that have the smallest distance and are from
//* different connected components will be connected by an edge, this is the
//* definition of bridge.
//*
//* @param ccs 2D-array of Leaf-IDs, where all IDs in the same list are in the
//*            same connected component (a group of nodes where all are
//*            visitable based on the edges between them).
//* @param neighbors A 2D-array of Leaf-IDs, where the entry of ID i at position
//*                  j indicates a edge between leaf j => i.
//*/
//template <class T>
//void SphereTree<T>::insertBridge(std::vector<std::vector<unsigned int>> *ccs,
//	std::vector<std::vector<unsigned int>> *neighbors)
//{
//	if ((*ccs).size() <= 1)
//	{
//		std::cout << "[pd] FAIL: tried to insert a bridge in a connected graph!" << std::endl;
//		return;
//	}
//
//	unsigned int global_nearest_id1, global_nearest_id2;
//	unsigned int global_nearest_cc1, global_nearest_cc2;
//	T global_nearest_distance = FLT_MAX;
//
//	for (unsigned int c1 = 0; c1<(*ccs).size(); c1++)
//	{
//		for (unsigned int i = 0; i<(*ccs)[c1].size(); i++)
//		{
//			// Check from each node
//			SphereNode<T> *leaf = &m_leaves[(*ccs)[c1][i]];
//
//			for (unsigned int c2 = 0; c2<(*ccs).size(); c2++)
//			{
//				if (c1 == c2)
//					continue;
//
//				for (unsigned int j = 0; j<(*ccs)[c2].size(); j++)
//				{
//					// To every node from another connected component
//					SphereNode<T> *leaf_ = &m_leaves[(*ccs)[c2][j]];
//					T distance(leaf->distance(*leaf_));
//					if (distance < global_nearest_distance)
//					{
//						global_nearest_distance = distance;
//						global_nearest_id1 = (*ccs)[c1][i];
//						global_nearest_id2 = (*ccs)[c2][j];
//						global_nearest_cc1 = c1;
//						global_nearest_cc2 = c2;
//					}
//				}
//			}
//		}
//	}
//
//	printf("[pd::SphereTree] inserting bridge: %u <---> %u\n", global_nearest_id1, global_nearest_id2);
//
//	// Insert undirected edge between nearest leaves
//	(*neighbors)[global_nearest_id1].push_back(global_nearest_id2);
//	(*neighbors)[global_nearest_id2].push_back(global_nearest_id1);
//
//	// Update connected components
//	(*ccs)[global_nearest_cc1].insert((*ccs)[global_nearest_cc1].end(), (*ccs)[global_nearest_cc2].begin(), (*ccs)[global_nearest_cc2].end());
//	ccs->erase(ccs->begin() + global_nearest_cc2);
//}

template class SphereTree<float>;
template class SphereTree<double>;

}
