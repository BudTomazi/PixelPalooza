#ifndef MARCHINGCUBES_H_
#define MARCHINGCUBES_H_

#include <vector>
#include "CGL/CGL.h"
#include "MCTable.h"		//tables used by Marching Cubes (edgeTable and triTable)
#include "CGL/vector3D.h"
#include "particle.h"

using namespace CGL;

typedef struct {
	Vector3D p[3];
    Vector3D norm[3];
    Vector3D colors[3];
    double shaderTypes[3];
} MeshTriangle;

struct ScalarLoc {
    Vector3D pos;
    float value;
    Vector3D color;
    int shaderType;
    
    ScalarLoc() : pos(Vector3D(0)), value(0), color(Vector3D(0)) { }
    
    ScalarLoc(double x, double y, double z, double value, Vector3D color) : pos(Vector3D(x, y, z)), value(value), color(color) { }
    
    ScalarLoc(Vector3D pos, double value, Vector3D color) : pos(pos), value(value), color(color) { }
    
    ScalarLoc(ScalarLoc& loc) : pos(loc.pos), value(loc.value), color(loc.color) { }
};

//does Linear Interpolation between points p1 and p2 (they already contain their computed values)
Vector3D LinearInterp(Vector4D p1, Vector4D p2, float value);
Vector3D LinearInterp(ScalarLoc p1, ScalarLoc p2, float value);

////////////////////////////////////////////////////////////////////////////////////////
//POINTERS TO FUNCTIONS																////
//pointer to function which computes the value at a point							////
typedef float (*FORMULA)(Vector3D);													////
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
///// the MARCHING CUBES algorithm /////											////
////////////////////////////////////////////////////////////////////////////////////////

//	Version 1
//Runs Marching Cubes on 3D array of points which already defined coordinates and values.
// ncellsX, ncellsY, ncellsZ):  number of cells to subdivide on each axis
// minValue: minimum to which value at every vertex is compared to
// points: array of length (ncellsX+1)(ncellsY+1)(ncellsZ+1) of mp4Vector points containing coordinates and values
//returns pointer to triangle array and the number of triangles in numTriangles
//note: array of points is first taken on z axis, then y and then x. So for example, if u iterate through it in a
//       for loop, have indexes i, j, k for x, y, z respectively, then to get the point you will have to make the
//		 following index: i*(ncellsY+1)*(ncellsZ+1) + j*(ncellsZ+1) + k .
//		Also, the array starts at the minimum on all axes. (point at indices (i,j,k) is the minimum)
MeshTriangle* MarchingCubesCross(int ncellsX, int ncellsY, int ncellsZ, float minValue, ScalarLoc* points, int& numTriangles);

//	Version 2
//First computes the 3D array of coordintes and values and then runs the above function on the data
//takes dimensions (minx,maxx,miny,...) and the number of cells (ncellsX,...) to subdivide on each axis
// minValue: minimum to which value at every vertex is compared to
// function of type float (mpVector p) formula, which computes value of p at its coordinates
// function of type mpVector (mp4Vector p1, mp4Vector p2) intersection, which determines the 
//  point of intersection of the surface and the edge between points p1 and p2
// saves number of triangles in numTriangles and the pointer to them is returned
// (note: mins' and maxs' are included in the algorithm)
MeshTriangle* MarchingCubesCross(float mcMinX, float mcMaxX, float mcMinY, float mcMaxY, float mcMinZ, float mcMaxZ,
	int ncellsX, int ncellsY, int ncellsZ, float minValue,
	FORMULA formula, int& numTriangles);


// 'STRAIGHT' Marching Cubes Algorithm //////////////////////////////////////////////////
//takes number of cells (ncellsX, ncellsY, ncellsZ) to subdivide on each axis
// minValue used to pass into LinearInterp
// gradFactor for each axis (multiplies each component of gradient vector by 1/(2*gradFactor) ).
//        Should be set to the length of a side (or close to it)
// array of length (ncellsX+1)(ncellsY+1)(ncellsZ+1) of mp4Vector points containing coordinates and values
//returns pointer to triangle array and the number of triangles in numTriangles
//note: array of points is first taken on z axis, then y and then x. So for example, if you iterate through it in a
//       for loop, have indexes i, j, k for x, y, z respectively, then to get the point you will have to make the
//         following index: i*(ncellsY+1)*(ncellsZ+1) + j*(ncellsZ+1) + k .
//        Also, the array starts at the minimum on all axes.
MeshTriangle* MarchingCubes(int ncellsX, int ncellsY, int ncellsZ, float gradFactorX, float gradFactorY, float gradFactorZ, float minValue, ScalarLoc* points, int &numTriangles);
/////////////////////////////////////////////////////////////////////////////////////////


#endif
