#include <iostream>
#include <vector>

#include <CGL/vector3D.h>
#include <CGL/vector4D.h>
#include "marchingCubes.h"

Vector3D LinearInterp(ScalarLoc p1, ScalarLoc p2, float value) {
	if (p1.value != p2.value)
		return p1.pos + (p2.pos - p1.pos) / (p2.value - p1.value) * (value - p1.value);
	else
		return p1.pos;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//	MARCHING CUBES	//
//////////////////////////////////////////////////////////////////////////////////////////////////////

//  VERSION  1  //
MeshTriangle* MarchingCubesCross(int ncellsX, int ncellsY, int ncellsZ, float minValue, ScalarLoc* points, int& numTriangles)
{
	//this should be enough space, if not change 3 to 4
    MeshTriangle* triangles = new MeshTriangle[3 * ncellsX * ncellsY * ncellsZ];
	numTriangles = int(0);

	int YtimeZ = (ncellsY + 1) * (ncellsZ + 1);	//for little extra speed
	int ni, nj;

	//go through all the points
	for (int i = 0; i < ncellsX; i++) {		//x axis
		ni = i * YtimeZ;
		for (int j = 0; j < ncellsY; j++) {	//y axis
			nj = j * (ncellsZ + 1);
			for (int k = 0; k < ncellsZ; k++)	//z axis
			{
				//initialize vertices
				ScalarLoc verts[8];
				int ind = ni + nj + k;
				/*(step 3)*/ verts[0] = points[ind];
				verts[1] = points[ind + YtimeZ];
				verts[2] = points[ind + YtimeZ + 1];
				verts[3] = points[ind + 1];
				verts[4] = points[ind + (ncellsZ + 1)];
				verts[5] = points[ind + YtimeZ + (ncellsZ + 1)];
				verts[6] = points[ind + YtimeZ + (ncellsZ + 1) + 1];
				verts[7] = points[ind + (ncellsZ + 1) + 1];

				//get the index
				int cubeIndex = int(0);
				for (int n = 0; n < 8; n++)
					/*(step 4)*/		if (verts[n].value <= minValue) cubeIndex |= (1 << n);

				//check if its completely inside or outside
				/*(step 5)*/ if (!edgeTable[cubeIndex]) continue;

				//get linearly interpolated vertices on edges and save into the array
				Vector3D intVerts[12];
				/*(step 6)*/ if (edgeTable[cubeIndex] & 1) intVerts[0] = LinearInterp(verts[0], verts[1], minValue);
				if (edgeTable[cubeIndex] & 2) intVerts[1] = LinearInterp(verts[1], verts[2], minValue);
				if (edgeTable[cubeIndex] & 4) intVerts[2] = LinearInterp(verts[2], verts[3], minValue);
				if (edgeTable[cubeIndex] & 8) intVerts[3] = LinearInterp(verts[3], verts[0], minValue);
				if (edgeTable[cubeIndex] & 16) intVerts[4] = LinearInterp(verts[4], verts[5], minValue);
				if (edgeTable[cubeIndex] & 32) intVerts[5] = LinearInterp(verts[5], verts[6], minValue);
				if (edgeTable[cubeIndex] & 64) intVerts[6] = LinearInterp(verts[6], verts[7], minValue);
				if (edgeTable[cubeIndex] & 128) intVerts[7] = LinearInterp(verts[7], verts[4], minValue);
				if (edgeTable[cubeIndex] & 256) intVerts[8] = LinearInterp(verts[0], verts[4], minValue);
				if (edgeTable[cubeIndex] & 512) intVerts[9] = LinearInterp(verts[1], verts[5], minValue);
				if (edgeTable[cubeIndex] & 1024) intVerts[10] = LinearInterp(verts[2], verts[6], minValue);
				if (edgeTable[cubeIndex] & 2048) intVerts[11] = LinearInterp(verts[3], verts[7], minValue);

				//now build the triangles using triTable
				for (int n = 0; triTable[cubeIndex][n] != -1; n += 3) {
					/*(step 7)*/ 	triangles[numTriangles].p[0] = intVerts[triTable[cubeIndex][n + 2]];
					triangles[numTriangles].p[1] = intVerts[triTable[cubeIndex][n + 1]];
					triangles[numTriangles].p[2] = intVerts[triTable[cubeIndex][n]];
					//Computing normal as cross product of triangle's edges
					/*(step 8)*/ 	triangles[numTriangles].norm =
						cross(triangles[numTriangles].p[1] - triangles[numTriangles].p[0], triangles[numTriangles].p[2] -
							triangles[numTriangles].p[0]).unit();
					numTriangles++;
				}

			}	//END OF Z FOR LOOP
		}	//END OF Y FOR LOOP
	}	//END OF X FOR LOOP

	//free all the wasted space
    MeshTriangle* retTriangles = new MeshTriangle[numTriangles];
	for (int i = 0; i < numTriangles; i++)
		retTriangles[i] = triangles[i];
	delete[] triangles;

	return retTriangles;
}

