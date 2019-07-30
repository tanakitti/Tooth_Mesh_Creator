#include <iostream>
#include "header/Vec3.h"
#include "header/SphereTree.h"
#include "header/Matrix4.h"
#include "header/SphereNode.h"
#include "header/Bankhelper.h"
#include "header/test.h"

int main(int argc, char *argv[]) {

    char * fileName = "../input/framediff46.log";
    if(argv[1]!= nullptr) fileName = argv[1];

    // register tooth
    auto M = pd::Matrix4<float>::IDENTITY();
    test::registerIdealEnvIST("../input/framefull.log",1,M,1,fileName);

    // Dimentional Setting
    Vec3<float> m_min( -90, -90, -90);
    Vec3<float> m_max( 90, 90, 90);

    Vec3<float> m_wd( m_max - m_min );
    Vec3<float> m_d(  100,100,100 );
    Vec3<float> m_df( m_d.x, m_d.y, m_d.z );
    Vec3<float> m_cellDim( m_wd/m_df );
    float *m_scalarFieldRaw = new float[m_d.x*m_d.y*m_d.z];
    Vec3<float> max(m_df);
    Vec3<float> min( 0, 0, 0);

    // Compute ScalarField
    test::cpukern_scalarField(m_min,m_max,m_wd,
                              m_cellDim, min, max, m_df,
                              m_scalarFieldRaw, true, test::g_idealEnvLeafCount, test::g_idealEnvLeafX,
                              test::g_idealEnvLeafY, test::g_idealEnvLeafZ, test::g_idealEnvLeafR,true);

    //
    int n = sizeof(m_scalarFieldRaw) / sizeof(m_scalarFieldRaw);
    std::reverse(m_scalarFieldRaw,m_scalarFieldRaw+n);

    test::MarchCubes(m_scalarFieldRaw,m_df,m_wd,m_wd);

    test::writeFile(test::m_gridToTriangles);

    test::convertObjtoDea();


    return 0;
}





