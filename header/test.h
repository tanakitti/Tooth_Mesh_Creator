//
// Created by Tanakit Sachati on 2019-07-15.
//

#ifndef UNTITLED4_TEST_H
#define UNTITLED4_TEST_H
#include <fstream>
#include <math.h>
#include "Vec3.h"
#include "Misc.h"
#include "Matrix4.h"
#include "SphereTree.h"
#include "Bankhelper.h"
#include "Config.h"
#include <unordered_map>

using namespace pd;

namespace test
{

    std::vector<SphereTree<float>*> g_idealEnvIsts;
    unsigned int g_idealEnvLeafCount = 0;
    float *g_idealEnvLeafX = nullptr;
    float *g_idealEnvLeafY = nullptr;
    float *g_idealEnvLeafZ = nullptr;
    float *g_idealEnvLeafR = nullptr;

    int numTriangle = 0;

    struct FCell
    {
        pd::Vec3<float> p[8];
        float val[8];
        int iFlagIndex = 0;
        int iEdgeFlags = 0;
    };

    std::unordered_map<int, std::vector<Vec3<float>>> m_gridToTriangles;
    Vec3<float>CubeSize;

    void registerIdealEnvIST( const char *fileName, int version, pd::Matrix4<float> M, float density, const char * fileName2){

        print("[registerIdealEnvIST] Register All Leaves");

        SphereTree<float> *env_ist = new SphereTree<float> ( fileName, version, M, density );

        if(fileName2 != ""){
            std::string filePath = fileName2;
            std::ifstream in(filePath);
            std::string line;
            std::ofstream out;
            int count = 0;

            while (std::getline(in, line))
            {
                if(count<3){
                    count++;
                    continue;
                }

                std::istringstream iss(line);
                std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss},std::istream_iterator<std::string>{} };

                unsigned int node_id = std::stoul(tokens[1]);
//                Vec3<float> center(std::stof(tokens[2]), std::stof(tokens[3]), std::stof(tokens[4]));
//                float radius = std::stof(tokens[5]);
                env_ist->m_leaves[node_id].m_radius = 0;

            }
        }



        g_idealEnvIsts.push_back( env_ist );

        // Write to CPU globals, to accumulate
        g_idealEnvLeafCount += env_ist->m_countLeaves;
        delete[] g_idealEnvLeafX;
        delete[] g_idealEnvLeafY;
        delete[] g_idealEnvLeafZ;
        delete[] g_idealEnvLeafR;

        g_idealEnvLeafX = new float[g_idealEnvLeafCount];
        g_idealEnvLeafY = new float[g_idealEnvLeafCount];
        g_idealEnvLeafZ = new float[g_idealEnvLeafCount];
        g_idealEnvLeafR = new float[g_idealEnvLeafCount];



        // Copy elemental components as SoA
        int cnt = 0;
        for (int i=0; i<g_idealEnvIsts.size(); i++) {
            for (int j=0; j<g_idealEnvIsts[i]->m_countLeaves; j++) {
                g_idealEnvLeafX[cnt] = g_idealEnvIsts[i]->m_leaves[j].m_center.x;
                g_idealEnvLeafY[cnt] = g_idealEnvIsts[i]->m_leaves[j].m_center.y;
                g_idealEnvLeafZ[cnt] = g_idealEnvIsts[i]->m_leaves[j].m_center.z;
                g_idealEnvLeafR[cnt] = g_idealEnvIsts[i]->m_leaves[j].m_radius;
                cnt++;
            }
        }
        print("[registerIdealEnvIST] Finished Register All Leaves\"");
    }

    void cpukern_scalarField( Vec3<float> bbStart, Vec3<float> bbEnd, Vec3<float> bbDiff, Vec3<float> cellDim, Vec3<float> min, Vec3<float> max, Vec3<float> gridRes, float *scalarField, bool full, unsigned long leafCount, float *leafX, float *leafY, float *leafZ, float *leafR ,bool logSpherePacking)
    {

        print("[cpukern_scalarField] Create Scalar Field");

        if (!full)
        {
            for (int z=min.z; z<=max.z-1; z++)
            {
                for (int y=min.y; y<=max.y-1; y++)
                {
                    for (int x=min.x; x<=max.x-1; x++)
                    {

                        int idx = x + y*gridRes.x + z*gridRes.x*gridRes.y;
                        scalarField[idx] = -0.1f;
//                        print(idx);
                    }
                }
            }
        } else {
            std::fill_n(scalarField, gridRes.x*gridRes.y*gridRes.z, -0.1f);
        }


        std::ofstream packing_user_txt;
        if (logSpherePacking) {
            packing_user_txt.open(std::string( "packing_user__" + GetRuntimeDate() + ".log" ), std::ofstream::trunc);
            WriteVersionInfo( packing_user_txt );
        }

        // Initalize back transformation for spheres (once)
        static Matrix4<float> M;
        static float scale;
        if (logSpherePacking) {
            static bool once = [](){
                M = Matrix4<float >::IDENTITY();
                if (false)
                    M(0,0) = -M(0,0);
                else
                    M(1,1) = -M(1,1);
                M.invert();
                scale = (M.extractScale().x+M.extractScale().y+M.extractScale().z) / 3.f;
                return true;
            } ();
        }

        print("[cpukern_scalarField] Precessing Scalar Field");

        // Thread input stride
        for (int s=0; s<leafCount; s++)
        {
            // Thread-specific sphere
            Vec3<float> ec( leafX[s], leafY[s], leafZ[s] );
            float er = leafR[s];

            // Find all grid cells of sphere intersection
            auto bbStartSphere = ec - Vec3<float >(1,1,1)*(er*2.f);
            auto bbEndSphere = ec + Vec3<float >(1,1,1)*(er*2.f);
            auto offsetStart = bbStartSphere - bbStart;
            auto offsetEnd = bbEndSphere - bbStart;
            Vec3<float> gridCellStart( offsetStart.x/cellDim.x,
                                       offsetStart.y/cellDim.y,
                                       offsetStart.z/cellDim.z );
            Vec3<float> gridCellEnd( offsetEnd.x/cellDim.x,
                                     offsetEnd.y/cellDim.y,
                                     offsetEnd.z/cellDim.z );

            if (!full) {
                gridCellStart.x = fmaxf( gridCellStart.x, static_cast<float >(min.x) );
                gridCellStart.y = fmaxf( gridCellStart.y, static_cast<float >(min.y) );
                gridCellStart.z = fmaxf( gridCellStart.z, static_cast<float >(min.z) );
                gridCellEnd.x = fminf( gridCellEnd.x, static_cast<float >(max.x) );
                gridCellEnd.y = fminf( gridCellEnd.y, static_cast<float >(max.y) );
                gridCellEnd.z = fminf( gridCellEnd.z, static_cast<float >(max.z) );
            }

            // Cap cells
            if (gridCellStart.x < 0) gridCellStart.x = 0;
            if (gridCellStart.y < 0) gridCellStart.y = 0;
            if (gridCellStart.z < 0) gridCellStart.z = 0;
            if (gridCellStart.x >= gridRes.x) gridCellStart.x = static_cast<float>(gridRes.x-1);
            if (gridCellStart.y >= gridRes.y) gridCellStart.y = static_cast<float >(gridRes.y-1);
            if (gridCellStart.z >= gridRes.z) gridCellStart.z = static_cast<float >(gridRes.z-1);
            if (gridCellEnd.x < 0) gridCellEnd.x = 0;
            if (gridCellEnd.y < 0) gridCellEnd.y = 0;
            if (gridCellEnd.z < 0) gridCellEnd.z = 0;
            if (gridCellEnd.x >= gridRes.x) gridCellEnd.x = static_cast<float >(gridRes.x-1);
            if (gridCellEnd.y >= gridRes.y) gridCellEnd.y = static_cast<float >(gridRes.y-1);
            if (gridCellEnd.z >= gridRes.z) gridCellEnd.z = static_cast<float >(gridRes.z-1);

            // Iterate over all grid cells in BB and write contribution
            for (int z=static_cast<int>(gridCellStart.z); z<=gridCellEnd.z; z++)
            {
                for (int y=static_cast<int>(gridCellStart.y); y<=gridCellEnd.y; y++)
                {
                    for (int x=static_cast<int>(gridCellStart.x); x<=gridCellEnd.x; x++)
                    {

                        if (er <= 0)
                            continue;

                        Vec3<float> gridDiff( static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) );
                        Vec3<float> worldVoxelSample = bbStart + gridDiff*cellDim;

                        float d = (worldVoxelSample-ec).calcLength()-er;
                        float p = d - (cellDim.calcLength()/1.73205080757f); //sqrt(3) =~ 1.73205080757


                        int idx = x + y*gridRes.x + z*gridRes.x*gridRes.y;
                        scalarField[idx] = fmaxf( -p, scalarField[idx]);


                    }
                }
            }
        }
        if (logSpherePacking)
            packing_user_txt.close();

        print("[cpukern_scalarField] Finished Precessing Scalar Field");
    }

    void CalculateFlagIndex(FCell& OutCell, const Vec3<float>& InCurrentCube, const Vec3<float>& InNumberOfCubes, const float* InScalarField, const int InLOD)
    {
        Vec3<float> BoundsMin( -90, -90, -90);

        OutCell.iFlagIndex = 0;
        for (int i = 0; i < 8; ++i)
        {
            // Calculate new index
            const int& Index = (InCurrentCube.z + cubeOffsets[i][2]) * (InNumberOfCubes.y) * (InNumberOfCubes.x)
                               + (InCurrentCube.y + cubeOffsets[i][1]) * (InNumberOfCubes.x)
                               + (InCurrentCube.x + cubeOffsets[i][0]);



            // The absolute position of all cubes corners
            OutCell.p[i].x = BoundsMin.x + CubeSize.x * (float)(InCurrentCube.x + cubeOffsets[i][0]);
            OutCell.p[i].y = BoundsMin.y + CubeSize.y * (float)(InCurrentCube.y + cubeOffsets[i][1]);
            OutCell.p[i].z = BoundsMin.z + CubeSize.z * (float)(InCurrentCube.z + cubeOffsets[i][2]);
            OutCell.val[i] = InScalarField[Index];



            // Generate a look up index for this specific edge intersection combination
            if (OutCell.val[i] <= 0.f )
                OutCell.iFlagIndex |= 1 << i;
        }

    }

    void Triangulate(FCell& OutCell, Vec3<float> OutTriangles[15], int& OutNumberOfTriangles, const int InLOD) {
        // Set edge flags for the current flag index
        OutCell.iEdgeFlags = aiCubeEdgeFlags[OutCell.iFlagIndex];

        // I don't exactly know in which case this happens, but let's just assume that we need these lines
        if (OutCell.iEdgeFlags == 0)
            return;

        // Compute actual edge vertices
        Vec3<float> EdgeVertices[12];

        for (int iEdge = 0; iEdge < 12; ++iEdge) {
            if (OutCell.iEdgeFlags & (1 << iEdge)) {

                // Get edge vertices indices
                const int iEdgeStart = a2iEdgeConnection[iEdge][0];
                const int iEdgeEnd = a2iEdgeConnection[iEdge][1];

                // Calculate difference between the two vertices
                const float DeltaVertices = OutCell.val[iEdgeEnd] - OutCell.val[iEdgeStart];

                // Calculate offset for vertices
                float VertexOffset = DeltaVertices == 0.0f ? 0.5f : ((0.f - OutCell.val[iEdgeStart]) / DeltaVertices);

//                float VertexOffset = 0.5f;
                // Set edge vertex
                EdgeVertices[iEdge] = OutCell.p[iEdgeStart] + VertexOffset * (OutCell.p[iEdgeEnd] - OutCell.p[iEdgeStart]);
            }
        }

        // Calculate triangles from edge data
        OutNumberOfTriangles = 0;
        for (int iTriangle = 0; iTriangle < 5; ++iTriangle)
        {
            // Not a valid triangle, break
            if (a2iTriangleConnectionTable[OutCell.iFlagIndex][3 * iTriangle] < 0)
                break;

            // Increase triangle count
            ++OutNumberOfTriangles;

            // Get vertex for each triangle
            for (int iCorner = 0; iCorner < 3; ++iCorner)
            {
                int iVertex = a2iTriangleConnectionTable[OutCell.iFlagIndex][3 * iTriangle + iCorner];
                OutTriangles[iTriangle * 3 + iCorner] = EdgeVertices[iVertex];
            }
        }
    }

    void MarchCubes(const float *scalarField, Vec3<float> GridRes, Vec3<float> BoundsSize, Vec3<float> GridDiff)
    {
        print("[MarchCubes] Creating Mesh ");
        CubeSize = Vec3<float>(GridDiff.x/(float)GridRes.x, GridDiff.y/(float)GridRes.y, GridDiff.z/(float)GridRes.z);
        Vec3<float> NumberOfCubes(
                BoundsSize.x / CubeSize.x,
                BoundsSize.y / CubeSize.y,
                BoundsSize.z / CubeSize.z
        );

        // Iterate over all
        Vec3<float> min = Vec3<float>( 0, 0, 0 );
        Vec3<float> max = NumberOfCubes - Vec3<float>(1,1,1);

        // March cubes
        Vec3<float> iCube( min );
        for (iCube.z = min.z; iCube.z < max.z; ++iCube.z)
        {
            for (iCube.y = min.y; iCube.y < max.y; ++iCube.y)
            {
                for (iCube.x = min.x; iCube.x < max.x; ++iCube.x)
                {
                    const int& Index = iCube.z * (NumberOfCubes.y) * (NumberOfCubes.z)
                                         + iCube.y * (NumberOfCubes.x)
                                         + iCube.x;

                    // Prepare for triangulation
                    FCell Cell;
                    int NumberOfTriangles = 0;
                    Vec3<float> CurrentTriangles[5 * 3];

                    // Triangulate cell
                    CalculateFlagIndex(Cell, iCube, NumberOfCubes, scalarField, 0);
                    Triangulate(Cell, CurrentTriangles, NumberOfTriangles, 0);

                    // Iterate over all triangles to fill vertex buffer and the index buffer
                    for (int iTriangle = 0; iTriangle < NumberOfTriangles; ++iTriangle)
                    {
                        // Iterate over the vertices of the triangle
                        for (int iVertex = 0; iVertex < 3; ++iVertex)
                        {
                            // Find index and normal
                            const Vec3<float>& VertexKey = CurrentTriangles[iTriangle * 3 + iVertex];
                            m_gridToTriangles[Index].push_back( VertexKey );

                            // count triangle
                            numTriangle++;

                        }
                    }
                }
            }
        }

        print("[MarchCubes] Finished Marchin Cueb");
        printf("[MarchCubes] Number of Generated Cube: ");
        print(test::numTriangle);

    }

    void writeFile(std::unordered_map<int, std::vector<Vec3<float>>> m_gridToTriangles){

        print("[writeFile] Saving Mesh");
        // open output file
        std::ofstream out("../output/tooth.obj");
        out<< "# ";
        out<< test::numTriangle;
        out<< "\n";
        int row = 1;

        // save vertrex into the output file
        for (auto index : m_gridToTriangles) {
            for (int i=0; i<index.second.size(); i++) {
                out<<"v "<<index.second[i].x << " " << index.second[i].y << " " << index.second[i].z << "\n";
                row++;
            }
        }

        // save face into the output file
        for(int f = 0; f< row-1;f+=3){
            out<<"f " << f+3 << " " << f+2 << " " << f+1 << "\n";
        }

        // save output file
        out.close();

        print("[writeFile] Finished Saving Mesh");
    }

    void convertObjtoDea(){


        time_t rawtime;
        struct tm * timeinfo;
        char buffer[80];

        time (&rawtime);
        timeinfo = localtime(&rawtime);

        strftime(buffer,sizeof(buffer),"%Y-%m-%dT%H:%M:%SZ",timeinfo);
        std::string time(buffer);
        std::string header = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                             "<COLLADA xmlns=\"http://www.collada.org/2005/11/COLLADASchema\" version=\"1.4.1\">\n"
                             " <asset>\n"
                             "  <contributor>\n"
                             "   <authoring_tool>SceneKit Collada Exporter v1.0</authoring_tool>\n"
                             "  </contributor>\n"
                             "  <created>"+time+"</created>\n"
                                                "  <modified>"+time+"</modified>\n"
                                                                    "  <unit meter=\"1.000000\"/>\n"
                                                                    "  <up_axis>Y_UP</up_axis>\n"
                                                                    " </asset>\n"
                                                                    " <library_materials>\n"
                                                                    "  <material id=\"\" name=\"\">\n"
                                                                    "   <instance_effect url=\"#effect_\"/>\n"
                                                                    "  </material>\n"
                                                                    " </library_materials>\n"
                                                                    " <library_effects>\n"
                                                                    "  <effect id=\"effect_\">\n"
                                                                    "   <profile_COMMON>\n"
                                                                    "    <technique sid=\"common\">\n"
                                                                    "     <phong>\n"
                                                                    "      <diffuse>\n"
                                                                    "       <color>1 1 1 1</color>\n"
                                                                    "      </diffuse>\n"
                                                                    "      <reflective>\n"
                                                                    "       <color>0 0 0 1</color>\n"
                                                                    "      </reflective>\n"
                                                                    "      <transparent opaque=\"A_ONE\">\n"
                                                                    "       <color>1 1 1 1</color>\n"
                                                                    "      </transparent>\n"
                                                                    "      <transparency>\n"
                                                                    "       <float>1</float>\n"
                                                                    "      </transparency>\n"
                                                                    "      <index_of_refraction>\n"
                                                                    "       <float>1</float>\n"
                                                                    "      </index_of_refraction>\n"
                                                                    "     </phong>\n"
                                                                    "    </technique>\n"
                                                                    "   </profile_COMMON>\n"
                                                                    "   <extra>\n"
                                                                    "    <technique profile=\"SceneKit\">\n"
                                                                    "     <litPerPixel>1</litPerPixel>\n"
                                                                    "     <ambient_diffuse_lock>1</ambient_diffuse_lock>\n"
                                                                    "    </technique>\n"
                                                                    "   </extra>\n"
                                                                    "  </effect>\n"
                                                                    " </library_effects>\n"
                                                                    " <library_geometries>\n"
                                                                    "  <geometry id=\"geometry1\">\n"
                                                                    "   <mesh>\n"
                                                                    "    <source id=\"geometrySource2\">\n"
                                                                    "     <float_array id=\"ID3-array\" count=\"";

        std::string header2 = "\">";


        std::string body = "</float_array>\n"
                           "<technique_common>\n"
                           "      <accessor source=\"#ID3-array\" count=\"";

        std::string body2 = "\" stride=\"3\">\n"
                            "       <param name=\"X\" type=\"float\"/>\n"
                            "       <param name=\"Y\" type=\"float\"/>\n"
                            "       <param name=\"Z\" type=\"float\"/>\n"
                            "      </accessor>\n"
                            "     </technique_common>\n"
                            "    </source>\n"
                            "    <vertices id=\"geometrySource2-vertices\">\n"
                            "     <input semantic=\"POSITION\" source=\"#geometrySource2\"/>\n"
                            "    </vertices>\n"
                            "    <triangles count=\"";
        std::string body3 = "\" material=\"geometryElement4\">\n"
                            "     <input semantic=\"VERTEX\" offset=\"0\" source=\"#geometrySource2-vertices\"/>\n"
                            "     <p>";

        std::string bottom = "</p>\n"
                             "    </triangles>\n"
                             "   </mesh>\n"
                             "  </geometry>\n"
                             " </library_geometries>\n"
                             " <library_visual_scenes>\n"
                             "  <visual_scene id=\"scene5\">\n"
                             "   <node id=\"MDL_OBJ\" name=\"MDL_OBJ\">\n"
                             "    <instance_geometry url=\"#geometry1\">\n"
                             "     <bind_material>\n"
                             "      <technique_common>\n"
                             "       <instance_material symbol=\"geometryElement4\" target=\"#\"/>\n"
                             "      </technique_common>\n"
                             "     </bind_material>\n"
                             "    </instance_geometry>\n"
                             "   </node>\n"
                             "  </visual_scene>\n"
                             " </library_visual_scenes>\n"
                             " <scene>\n"
                             "  <instance_visual_scene url=\"#scene5\"/>\n"
                             " </scene>\n"
                             "</COLLADA>";

        print("[convertObjtoDea] Convert Object to Dea file format");
        std::string filePath = "../output/tooth.obj";

        std::ifstream in(filePath);
        std::string line;
        std::ofstream out("../output/tooth.dae");
        int count = 0;



        int i = 0;

        while (std::getline(in, line))
        {
            std::istringstream iss(line);
            std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss},std::istream_iterator<std::string>{} };

            if(tokens[0]=="f" ){
                if(count==0){
                    out << body;
                    count++;
                }
                out<< i++;
                out<< " ";
                out<< i++;
                out<< " ";
                out<< i++;
                out<< " ";
            }else if(tokens[0]=="v"){
                out<< (tokens[3])+" ";
                out<< (tokens[2])+" ";
                out<< (tokens[1])+" ";
            }else{

                int numVertrex = std::stoi(tokens[1]);
                int numXYZ = numVertrex*3;
                int numTriangle = numVertrex/3;

                header = header+toString(numXYZ)+header2;
                body = body+toString(numXYZ)+body2+toString(numTriangle)+body3;
                out << header;
            }
        }

        out<< bottom;
        out.close();

        print("[convertObjtoDea] Finished Converting Object to Dea File Format");
    }
}
#endif //UNTITLED4_TEST_H
