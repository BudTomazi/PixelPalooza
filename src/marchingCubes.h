#ifndef MARCHINGCUBES_H_
#define MARCHINGCUBES_H_

#include <vector>
#include "CGL/CGL.h"
#include "MCTable.h"		//tables used by Marching Cubes (edgeTable and triTable)
#include "CGL/vector3D.h"

using namespace CGL;

//used to save triangles - 3 vertices and a normal vector
typedef struct {
	Vector3D p[3];
    Vector3D norm[3];
} MeshTriangle;

struct ScalarLoc {
    Vector3D pos;
    float value;
    
    ScalarLoc() : pos(Vector3D(0)), value(0) { }
    
    ScalarLoc(double x, double y, double z, double value) : pos(Vector3D(x, y, z)), value(value) { }
    
    ScalarLoc(Vector3D pos, double value) : pos(pos), value(value) { }
    
    ScalarLoc(ScalarLoc& loc) : pos(loc.pos), value(loc.value) { }
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


// RECURSIVE Marching Cubes    /////////////////////////////////////////////////////////////
//DrawBacks: must know how many objects are drawn

//This function starts the Recursive Marching Cubes
// numCubes: number of intersecting cubes passed as starting points
// ii, jj, kk: arrays of size numCubes. Contain respective indexes of cubes that are intersected
MeshTriangle* MarchingCubesRec(int ncellsX, int ncellsY, int ncellsZ,
                            float gradFactorX, float gradFactorY, float gradFactorZ,
                            int numCubes, int *ii, int *jj, int *kk,
                            float minValue, ScalarLoc* points, int &numTriangles);

//Next 6 functions are called by the corresponding face of each cube (e.g. face 0 calls MCFace0 etc...)
// Each function accepts the following information as arguments:
//    First 3 arguments: number of cells to subdivide on each axis
//    gradFactor: factor to scale the gradient vectors (multiplies by 1/(2*gradFactor) )
//    ind: index of this cube in the points array (this is so it doesnt have to be computed again)
//    i,j,k: indexes of the this cube on x,y,z axis respectively
//    minValue: minimum used in LinearInterp
//    points: array of type mp4Vector and size (ncellsX+1)*(ncellsY+1)*(ncellsZ+1) that contains
//        coordinates and values at each point
//    triangle: array of triangles which is being built and returned at the end of recursion
//    numTriangles: number of triangles formed
//    prevVerts: adjacent 4 vertices of the previous cube, passed in the special order which the correspoding
//        MCFace function recognizes. For specificatino on which indexes passed from the last cube go to
//        which vertexes of the current cube visit my website at www.angelfire.com/linux/myp
//    prevIntVerts: array of 4 linearly interpolated vertices on 4 adjacent edges
//    edgeIndex: integer, bits of which correspond to the edges of the current cube which have been computed
//        from the previous cube
//    prevGradVerts: array of 4 gradient vectors at 4 adjacent vertices
//    prevGrads: linearly interpolated gradient vectors on 4 adjacent edges
//    gradIndex: integer bits of which correspond to already computed vertices of the current cube
//    marchedCubes: bool array of size ncellsX*ncellsY*ncellsZ which stores which cubes already have been marched
//        through. Initialized to FALSE. This is not the most effective way in terms of memory management, but
//        faster than using STLs vector<bool> which stores booleans as bits.
//* NOTE *: Each of these 6 functions run marching cubes on the 'next'cube adjacent to the surface number
//        that is specified in their names (for numbering see my webpage). Each assigns the previous computed
//        values to the face 'opposite' of that specified. Then each one runs Marching Cubes for the current cube.
//        It returns doing nothing if the cube is not intersected, however, if it is then other recursive
//        functions are called using macros MC_FACE defined for each face. See MarchingCubes.cpp
//        For example MCFace0 will initialize surface 2 using the previous values that are passed to it. Then if
//        the current cube is intersected it will run MC_FACE0, MC_FACE1, MC_FACE3, MC_FACE4, and MC_FACE5. Each
//        returnes the pointer to the array of formed triangles.

//FACE 0 Marching Cubes
MeshTriangle* MCFace0(int ncellsX, int ncellsY, int ncellsZ,
                        float gradFactorX, float gradFactorY, float gradFactorZ,
                        int ind, int i, int j, int k,
                        float minValue, ScalarLoc* points, MeshTriangle* triangles, int &numTriangles,
                      ScalarLoc prevVerts[4], Vector3D prevIntVerts[4], int edgeIndex,
                        Vector4D prevGradVerts[4], Vector3D prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 1 Marching Cubes
MeshTriangle* MCFace1(int ncellsX, int ncellsY, int ncellsZ,
                        float gradFactorX, float gradFactorY, float gradFactorZ,
                        int ind, int i, int j, int k,
                        float minValue, ScalarLoc* points, MeshTriangle* triangles, int &numTriangles,
                      ScalarLoc prevVerts[4], Vector3D prevIntVerts[4], int edgeIndex,
                      Vector4D prevGradVerts[4], Vector3D prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 2 Marching Cubes
MeshTriangle* MCFace2(int ncellsX, int ncellsY, int ncellsZ,
                        float gradFactorX, float gradFactorY, float gradFactorZ,
                        int ind, int i, int j, int k,
                        float minValue, ScalarLoc* points, MeshTriangle* triangles, int &numTriangles,
                      ScalarLoc prevVerts[4], Vector3D prevIntVerts[4], int edgeIndex,
                      Vector4D prevGradVerts[4], Vector3D prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 3 Marching Cubes
MeshTriangle* MCFace3(int ncellsX, int ncellsY, int ncellsZ,
                        float gradFactorX, float gradFactorY, float gradFactorZ,
                        int ind, int i, int j, int k,
                        float minValue, ScalarLoc* points, MeshTriangle* triangles, int &numTriangles,
                      ScalarLoc prevVerts[4], Vector3D prevIntVerts[4], int edgeIndex,
                      Vector4D prevGradVerts[4], Vector3D prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 4 Marching Cubes
MeshTriangle* MCFace4(int ncellsX, int ncellsY, int ncellsZ,
                        float gradFactorX, float gradFactorY, float gradFactorZ,
                        int ind, int i, int j, int k,
                        float minValue, ScalarLoc* points, MeshTriangle* triangles, int &numTriangles,
                      ScalarLoc prevVerts[4], Vector3D prevIntVerts[4], int edgeIndex,
                      Vector4D prevGradVerts[4], Vector3D prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 5 Marching Cubes
MeshTriangle* MCFace5(int ncellsX, int ncellsY, int ncellsZ,
                        float gradFactorX, float gradFactorY, float gradFactorZ,
                        int ind, int i, int j, int k,
                        float minValue, ScalarLoc* points, MeshTriangle* triangles, int &numTriangles,
                      ScalarLoc prevVerts[4], Vector3D prevIntVerts[4], int edgeIndex,
                      Vector4D prevGradVerts[4], Vector3D prevGrads[4], int gradIndex, bool* marchedCubes);



//Does Marching Cubes on cube (i, j, k) in points
// Requirements: needs index ind which specifies which cube it is (its 0th point in points array)
//    If cube is intersected the triangles are added to triangles array passed to it, and numTriangles is
//    incremented appropriately. verts, array of vertices of this cube should be initialized.
//    intVerts, interpolated vertices on edges, gradVerts, gradients on vertices, and grads, interpolated
//    gradient vectors on edges dont have to be initialized, but if they are corresponding bits in edgeIndex
//    and indGrad index should show it. (for numbering see my web page)
//    Global variables YtimeZ should be initialized to (ncellsZ+1)*(ncellsY+1) and pointsZ should be ncellsZ+1
MeshTriangle* MarchOneCube(int ncellsX, int ncellsY, int ncellsZ,
                        float gradFactorX, float gradFactorY, float gradFactorZ,
                        int ind, int i, int j, int k,
                        float minValue, ScalarLoc* points, MeshTriangle* triangles, int &numTriangles,
                        ScalarLoc verts[8], Vector3D intVerts[12], int &edgeIndex,
                           Vector4D gradVerts[8], Vector3D grads[12], int &indGrad);


//Find the first intersecting cube (if it exists).
//Starts at (i,j,k)=(0,0,0) at the minimum x,y,z and then goes linearly searching for an intersecting cube.
//Returns array of indices i,j,k of the first found cube or -1,-1,-1 if cube was not found
float* MCFind(int ncellsX, int ncellsY, int ncellsZ, float minValue, ScalarLoc* points);

//Calls the MCFind and if found calls MarchingCubesRec with the returned indices
//If not found NULL is returned and numTriangles is set to ZERO.
MeshTriangle* MCRecFind(int ncellsX, int ncellsY, int ncellsZ,
                        float gradFactorX, float gradFactorY, float gradFactorZ,
                        float minValue, ScalarLoc* points, int &numTriangles);

#endif
