/*
 * MyTraits.h
 *
 *  Created on: Mar 11, 2016
 *      Author: darryl
 */

#ifndef MYTRAITS_H_
#define MYTRAITS_H_


#include <iostream>
#include "priorityqueue.h"
#include <functional>
#include <math.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

using namespace std;
// Define my personal traits

struct MyTraits : OpenMesh::DefaultTraits
{
    // Let Point and Normal be a vector of doubles
    typedef OpenMesh::Vec3d Point;
    typedef OpenMesh::Vec3d Normal;
    typedef OpenMesh::Vec4f Color;


    FaceAttributes(OpenMesh::Attributes::Color);
    EdgeAttributes(OpenMesh::Attributes::Color);
    // Already defined in OpenMesh::DefaultTraits
    // HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );

    // Uncomment next line to disable attribute PrevHalfedge
    // HalfedgeAttributes( OpenMesh::Attributes::None );
    //
    // or
    //
    // HalfedgeAttributes( 0 );
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> Mesh;
typedef OpenMesh::EPropHandleT< vector<Mesh::Point>*  > Prop;
//because bool defaults to true, we use int, which defaults to 0, and we can define
//however we like
typedef OpenMesh::EPropHandleT<int> Delaunay_indicator;
typedef OpenMesh::EPropHandleT<int> Flippable;



#endif /* MYTRAITS_H_ */
