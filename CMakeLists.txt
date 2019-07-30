cmake_minimum_required(VERSION 3.14)
project(untitled8)

set(CMAKE_CXX_STANDARD 14)

add_executable(untitled8 main.cpp header/Vec3.h resouce/Vec3.cpp header/Resources.h header/SphereNode.h resouce/SphereNode.cpp header/SphereTree.h header/Matrix4.h resouce/Matrix4.cpp resouce/SphereTree.cpp header/Misc.h header/buildparams.h header/Bankhelper.h header/Kernels.h)