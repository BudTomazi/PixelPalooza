/////////////////////////////////////////////////////////////////////////////////////////////
//    FileName:    MarchingCubes.cpp
//    Author    :    Michael Y. Polyakov
//    email    :    myp@andrew.cmu.edu  or  mikepolyakov@hotmail.com
//    website    :    www.angelfire.com/linux/myp
//    date    :    July 2002
//
//    Description:    'Straight' and Recursive Marching Cubes Algorithms
//                    Recursive method is faster than the 'straight' one, especially when intersection does not
//                        have to be searched for every time.
//                Normal vectors are defined for each vertex as a gradients
//                For function definitions see MarchingCubes.h
//                For a tutorial on Marching Cubes please visit www.angelfire.com/myp/linux
//
//    Please email me with any suggestion/bugs.
/////////////////////////////////////////////////////////////////////////////////////////////

#include "marchingCubes.h"
#include "MCTable.h"
#include <math.h>

#include <vector>
#include "CGL/vector3D.h"
#include "CGL/vector4D.h"

#include <iostream>

//Linear Interpolation between two points
Vector3D LinearInterp(Vector4D p1, Vector4D p2, float value)
{
    if (abs(p1.w - p2.w) > 0.00001)
        return p1.to3D() + (p2.to3D() - p1.to3D()) / (p2.w - p1.w) * (value - p1.w);
    else
        return p1.to3D();
}

Vector3D LinearInterp(ScalarLoc p1, ScalarLoc p2, float value)
{
    if (abs(p1.value - p2.value) > 0.00001)
        return p1.pos + (p2.pos - p1.pos) / (p2.value - p1.value) * (value - p1.value);
    else
        return p1.pos;
}

Vector3D ColorInterp(ScalarLoc p1, ScalarLoc p2, float value) {
    if (p1.value > 0 && p2.value > 0) {
        return p1.color + (p2.color - p1.color) / (p2.value - p1.value) * (value - p1.value);
    } else if (p1.value > 0) {
        return p1.color;
    } else { // doens't matter
        return p2.color;
    }
}

int GetShaderType(ScalarLoc p1, ScalarLoc p2, float value) {
    if (p1.shaderType != p2.shaderType) return p1.shaderType;
    
    if (p1.value > 0 && p2.value > 0) {
        if (abs(value - p1.value) < abs(value - p2.value)) {
            return p1.shaderType; // closer to p1
        } else {
            return p2.shaderType;
        }
    } else if (p1.value > 0) {
        return p1.shaderType;
    } else { // doens't matter
        return p2.shaderType;
    }
}

//Macros used to compute gradient vector on each vertex of a cube
//argument should be the name of array of vertices
//can be verts or *verts if done by reference
#define CALC_GRAD_VERT_0(verts) Vector4D(points[ind-YtimeZ].value-(verts[1]).value,points[ind-pointsZ].value-(verts[4]).value,points[ind-1].value-(verts[3]).value, (verts[0]).value);
#define CALC_GRAD_VERT_1(verts) Vector4D((verts[0]).value-points[ind+2*YtimeZ].value,points[ind+YtimeZ-pointsZ].value-(verts[5]).value,points[ind+YtimeZ-1].value-(verts[2]).value, (verts[1]).value);
#define CALC_GRAD_VERT_2(verts) Vector4D((verts[3]).value-points[ind+2*YtimeZ+1].value,points[ind+YtimeZ-ncellsZ].value-(verts[6]).value,(verts[1]).value-points[ind+YtimeZ+2].value, (verts[2]).value);
#define CALC_GRAD_VERT_3(verts) Vector4D(points[ind-YtimeZ+1].value-(verts[2]).value,points[ind-ncellsZ].value-(verts[7]).value,(verts[0]).value-points[ind+2].value, (verts[3]).value);
#define CALC_GRAD_VERT_4(verts) Vector4D(points[ind-YtimeZ+ncellsZ+1].value-(verts[5]).value,(verts[0]).value-points[ind+2*pointsZ].value,points[ind+ncellsZ].value-(verts[7]).value, (verts[4]).value);
#define CALC_GRAD_VERT_5(verts) Vector4D((verts[4]).value-points[ind+2*YtimeZ+ncellsZ+1].value,(verts[1]).value-points[ind+YtimeZ+2*pointsZ].value,points[ind+YtimeZ+ncellsZ].value-(verts[6]).value, (verts[5]).value);
#define CALC_GRAD_VERT_6(verts) Vector4D((verts[7]).value-points[ind+2*YtimeZ+ncellsZ+2].value,(verts[2]).value-points[ind+YtimeZ+2*ncellsZ+3].value,(verts[5]).value-points[ind+YtimeZ+ncellsZ+3].value, (verts[6]).value);
#define CALC_GRAD_VERT_7(verts) Vector4D(points[ind-YtimeZ+ncellsZ+2].value-(verts[6]).value,(verts[3]).value-points[ind+2*ncellsZ+3].value,(verts[4]).value-points[ind+ncellsZ+3].value, (verts[7]).value);

///////////////////////////////////////////////////////////////////////////////////////////////////////
// GLOBAL //
//Global variables - so they dont have to be passed into functions
int pointsZ;    //number of points on Z zxis (equal to ncellsZ+1)
int YtimeZ;        //'plane' of cubes on YZ (equal to (ncellsY+1)*pointsZ )
///////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
//    'STRAIGHT' MARCHING CUBES    ALGORITHM  ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//for gradients at the edges values of 1.0, 1.0, 1.0, 1.0  are given
MeshTriangle* MarchingCubes(int ncellsX, int ncellsY, int ncellsZ,
                            float gradFactorX, float gradFactorY, float gradFactorZ,
                        float minValue, ScalarLoc* points, int &numTriangles)
{
    //this should be enough space, if not change 3 to 4
    MeshTriangle* triangles = new MeshTriangle[3 * ncellsX * ncellsY * ncellsZ];
    numTriangles = int(0);

    pointsZ = ncellsZ + 1;            //initialize global variable (for extra speed)
    YtimeZ = (ncellsY + 1) * pointsZ;
    int lastX = ncellsX;            //left from older version
    int lastY = ncellsY;
    int lastZ = ncellsZ;

    ScalarLoc* verts[8];            //vertices of a cube (array of pointers for extra speed)
    Vector3D intVerts[12];            //linearly interpolated vertices on each edge
    Vector3D intColors[12];
    int shaderTypes[12];
    int cubeIndex;                    //shows which vertices are outside/inside
    int edgeIndex;                    //index returned by edgeTable[cubeIndex]
    Vector4D gradVerts[8];            //gradients at each vertex of a cube
    Vector3D grads[12];                //linearly interpolated gradients on each edge
    int indGrad;                    //shows which gradients already have been computed
    int ind, ni, nj;                //ind: index of vertex 0
    //factor by which corresponding coordinates of gradient vectors are scaled
    Vector3D factor = Vector3D(1.0 / (2.0 * gradFactorX),
                               1.0 / (2.0 * gradFactorY),
                               1.0 / (2.0 * gradFactorZ));

    //DebugWrite("1");

    //MAIN LOOP: goes through all the points
    for(int i=0; i < lastX; i++) {            //x axis
        ni = i*YtimeZ;
        for(int j=0; j < lastY; j++) {        //y axis
            nj = j*pointsZ;
            for(int k=0; k < lastZ; k++, ind++)    //z axis
            {
                //initialize vertices
                ind = ni + nj + k;
                verts[0] = &points[ind];
                verts[1] = &points[ind + YtimeZ];
                verts[4] = &points[ind + pointsZ];
                verts[5] = &points[ind + YtimeZ + pointsZ];
                verts[2] = &points[ind + YtimeZ + 1];
                verts[3] = &points[ind + 1];
                verts[6] = &points[ind + YtimeZ + pointsZ + 1];
                verts[7] = &points[ind + pointsZ + 1];

                //get the index
                cubeIndex = int(0);
                for(int n=0; n < 8; n++)
                    if(verts[n]->value <= minValue) cubeIndex |= (1 << n);

                //check if its completely inside or outside
                if(!edgeTable[cubeIndex]) continue;

                indGrad = int(0);
                edgeIndex = edgeTable[cubeIndex];

                if(edgeIndex & 1) {
                    intVerts[0] = LinearInterp(*verts[0], *verts[1], minValue);
                    intColors[0] = ColorInterp(*verts[0], *verts[1], minValue);
                    shaderTypes[0] = GetShaderType(*verts[0], *verts[1], minValue);
                    
                    if(i != 0 && j != 0 && k != 0) gradVerts[0] = CALC_GRAD_VERT_0(*verts)
                    else gradVerts[0] = Vector4D(1.0, 1.0, 1.0, 1.0);
                    if(i != lastX-1 && j != 0 && k != 0) gradVerts[1] = CALC_GRAD_VERT_1(*verts)
                    else gradVerts[1] = Vector4D(1.0, 1.0, 1.0, 1.0);
                    indGrad |= 3;
                    grads[0] = LinearInterp(gradVerts[0], gradVerts[1], minValue);
                    grads[0].x *= factor.x; grads[0].y *= factor.y; grads[0].z *= factor.z;
                }
                if(edgeIndex & 2) {
                    intVerts[1] = LinearInterp(*verts[1], *verts[2], minValue);
                    intColors[1] = ColorInterp(*verts[1], *verts[2], minValue);
                    shaderTypes[1] = GetShaderType(*verts[1], *verts[2], minValue);
                    
                    if(! (indGrad & 2)) {
                        if(i != lastX-1 && j != 0 && k != 0) gradVerts[1] = CALC_GRAD_VERT_1(*verts)
                        else gradVerts[1] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 2;
                    }
                    if(i != lastX-1 && j != 0 && k != 0) gradVerts[2] = CALC_GRAD_VERT_2(*verts)
                    else gradVerts[2] = Vector4D(1.0, 1.0, 1.0, 1.0);
                    indGrad |= 4;
                    grads[1] = LinearInterp(gradVerts[1], gradVerts[2], minValue);
                    grads[1].x *= factor.x; grads[1].y *= factor.y; grads[1].z *= factor.z;
                }
                if(edgeIndex & 4) {
                    intVerts[2] = LinearInterp(*verts[2], *verts[3], minValue);
                    intColors[2] = ColorInterp(*verts[2], *verts[3], minValue);
                    shaderTypes[2] = GetShaderType(*verts[2], *verts[3], minValue);

                    if(! (indGrad & 4)) {
                        if(i != lastX-1 && j != 0 && k != 0) gradVerts[2] = CALC_GRAD_VERT_2(*verts)
                        else gradVerts[2] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 4;
                    }
                    if(i != 0 && j != 0 && k != lastZ-1) gradVerts[3] = CALC_GRAD_VERT_3(*verts)
                    else gradVerts[3] = Vector4D(1.0, 1.0, 1.0, 1.0);
                    indGrad |= 8;
                    grads[2] = LinearInterp(gradVerts[2], gradVerts[3], minValue);
                    grads[2].x *= factor.x; grads[2].y *= factor.y; grads[2].z *= factor.z;
                }
                if(edgeIndex & 8) {
                    intVerts[3] = LinearInterp(*verts[3], *verts[0], minValue);
                    intColors[3] = ColorInterp(*verts[3], *verts[0], minValue);
                    shaderTypes[3] = GetShaderType(*verts[3], *verts[0], minValue);

                    if(! (indGrad & 8)) {
                        if(i != 0 && j != 0 && k != lastZ-1) gradVerts[3] = CALC_GRAD_VERT_3(*verts)
                        else gradVerts[3] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 8;
                    }
                    if(! (indGrad & 1)) {
                        if(i != 0 && j != 0 && k != 0) gradVerts[0] = CALC_GRAD_VERT_0(*verts)
                        else gradVerts[0] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 1;
                    }
                    grads[3] = LinearInterp(gradVerts[3], gradVerts[0], minValue);
                    grads[3].x *= factor.x; grads[3].y *= factor.y; grads[3].z *= factor.z;
                }
                if(edgeIndex & 16) {
                    intVerts[4] = LinearInterp(*verts[4], *verts[5], minValue);
                    intColors[4] = ColorInterp(*verts[4], *verts[5], minValue);
                    shaderTypes[4] = GetShaderType(*verts[4], *verts[5], minValue);

                    if(i != 0 && j != lastY-1 && k != 0) gradVerts[4] = CALC_GRAD_VERT_4(*verts)
                    else gradVerts[4] = Vector4D(1.0, 1.0, 1.0, 1.0);

                    if(i != lastX-1 && j != lastY-1 && k != 0) gradVerts[5] = CALC_GRAD_VERT_5(*verts)
                    else gradVerts[5] = Vector4D(1.0, 1.0, 1.0, 1.0);

                    indGrad |= 48;
                    grads[4] = LinearInterp(gradVerts[4], gradVerts[5], minValue);
                    grads[4].x *= factor.x; grads[4].y *= factor.y; grads[4].z *= factor.z;
                }
                if(edgeIndex & 32) {
                    intVerts[5] = LinearInterp(*verts[5], *verts[6], minValue);
                    intColors[5] = ColorInterp(*verts[5], *verts[6], minValue);
                    shaderTypes[5] = GetShaderType(*verts[5], *verts[6], minValue);

                    if(! (indGrad & 32)) {
                        if(i != lastX-1 && j != lastY-1 && k != 0) gradVerts[5] = CALC_GRAD_VERT_5(*verts)
                        else gradVerts[5] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 32;
                    }

                    if(i != lastX-1 && j != lastY-1 && k != lastZ-1) gradVerts[6] = CALC_GRAD_VERT_6(*verts)
                    else gradVerts[6] = Vector4D(1.0, 1.0, 1.0, 1.0);
                    indGrad |= 64;
                    grads[5] = LinearInterp(gradVerts[5], gradVerts[6], minValue);
                    grads[5].x *= factor.x; grads[5].y *= factor.y; grads[5].z *= factor.z;
                }
                if(edgeIndex & 64) {
                    intVerts[6] = LinearInterp(*verts[6], *verts[7], minValue);
                    intColors[6] = ColorInterp(*verts[6], *verts[7], minValue);
                    shaderTypes[6] = GetShaderType(*verts[6], *verts[7], minValue);

                    if(! (indGrad & 64)) {
                        if(i != lastX-1 && j != lastY-1 && k != lastZ-1) gradVerts[6] = CALC_GRAD_VERT_6(*verts)
                        else gradVerts[6] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 64;
                    }

                    if(i != 0 && j != lastY-1 && k != lastZ-1) gradVerts[7] = CALC_GRAD_VERT_7(*verts)
                    else gradVerts[7] = Vector4D(1.0, 1.0, 1.0, 1.0);
                    indGrad |= 128;
                    grads[6] = LinearInterp(gradVerts[6], gradVerts[7], minValue);
                    grads[6].x *= factor.x; grads[6].y *= factor.y; grads[6].z *= factor.z;
                }
                if(edgeIndex & 128) {
                    intVerts[7] = LinearInterp(*verts[7], *verts[4], minValue);
                    intColors[7] = ColorInterp(*verts[7], *verts[4], minValue);
                    shaderTypes[7] = GetShaderType(*verts[7], *verts[4], minValue);

                    if(! (indGrad & 128)) {
                        if(i != 0 && j != lastY-1 && k != lastZ-1) gradVerts[7] = CALC_GRAD_VERT_7(*verts)
                        else gradVerts[7] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 128;
                    }
                    if(! (indGrad & 16)) {
                        if(i != 0 && j != lastY-1 && k != 0) gradVerts[4] = CALC_GRAD_VERT_4(*verts)
                        else gradVerts[4] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 16;
                    }
                    grads[7] = LinearInterp(gradVerts[7], gradVerts[4], minValue);
                    grads[7].x *= factor.x; grads[7].y *= factor.y; grads[7].z *= factor.z;
                }
                if(edgeIndex & 256) {
                    intVerts[8] = LinearInterp(*verts[0], *verts[4], minValue);
                    intColors[8] = ColorInterp(*verts[0], *verts[4], minValue);
                    shaderTypes[8] = GetShaderType(*verts[0], *verts[4], minValue);

                    if(! (indGrad & 1)) {
                        if(i != 0 && j != 0 && k != 0) gradVerts[0] = CALC_GRAD_VERT_0(*verts)
                        else gradVerts[0] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 1;
                    }
                    if(! (indGrad & 16)) {
                        if(i != 0 && j != lastY-1 && k != 0) gradVerts[4] = CALC_GRAD_VERT_4(*verts)
                        else gradVerts[4] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 16;
                    }
                    grads[8] = LinearInterp(gradVerts[0], gradVerts[4], minValue);
                    grads[8].x *= factor.x; grads[8].y *= factor.y; grads[8].z *= factor.z;
                }
                if(edgeIndex & 512) {
                    intVerts[9] = LinearInterp(*verts[1], *verts[5], minValue);
                    intColors[9] = ColorInterp(*verts[1], *verts[5], minValue);
                    shaderTypes[9] = GetShaderType(*verts[1], *verts[5], minValue);

                    if(! (indGrad & 2)) {
                        if(i != lastX-1 && j != 0 && k != 0) gradVerts[1] = CALC_GRAD_VERT_1(*verts)
                        else gradVerts[1] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 2;
                    }
                    if(! (indGrad & 32)) {
                        if(i != lastX-1 && j != lastY-1 && k != 0) gradVerts[5] = CALC_GRAD_VERT_5(*verts)
                        else gradVerts[5] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 32;
                    }
                    grads[9] = LinearInterp(gradVerts[1], gradVerts[5], minValue);
                    grads[9].x *= factor.x; grads[9].y *= factor.y; grads[9].z *= factor.z;
                }
                if(edgeIndex & 1024) {
                    intVerts[10] = LinearInterp(*verts[2], *verts[6], minValue);
                    intColors[10] = ColorInterp(*verts[2], *verts[6], minValue);
                    shaderTypes[10] = GetShaderType(*verts[2], *verts[6], minValue);

                    if(! (indGrad & 4)) {
                        if(i != lastX-1 && j != 0 && k != 0) gradVerts[2] = CALC_GRAD_VERT_2(*verts)
                        else gradVerts[5] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 4;
                    }
                    if(! (indGrad & 64)) {
                        if(i != lastX-1 && j != lastY-1 && k != lastZ-1) gradVerts[6] = CALC_GRAD_VERT_6(*verts)
                        else gradVerts[6] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 64;
                    }
                    grads[10] = LinearInterp(gradVerts[2], gradVerts[6], minValue);
                    grads[10].x *= factor.x; grads[10].y *= factor.y; grads[10].z *= factor.z;
                }
                if(edgeIndex & 2048) {
                    intVerts[11] = LinearInterp(*verts[3], *verts[7], minValue);
                    intColors[11] = ColorInterp(*verts[3], *verts[7], minValue);
                    shaderTypes[11] = GetShaderType(*verts[3], *verts[7], minValue);

                    if(! (indGrad & 8)) {
                        if(i != 0 && j != 0 && k != lastZ-1) gradVerts[3] = CALC_GRAD_VERT_3(*verts)
                        else gradVerts[3] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 8;
                    }
                    if(! (indGrad & 128)) {
                        if(i != 0 && j != lastY-1 && k != lastZ-1) gradVerts[7] = CALC_GRAD_VERT_7(*verts)
                        else gradVerts[7] = Vector4D(1.0, 1.0, 1.0, 1.0);
                        indGrad |= 128;
                    }
                    grads[11] = LinearInterp(gradVerts[3], gradVerts[7], minValue);
                    grads[11].x *= factor.x; grads[11].y *= factor.y; grads[11].z *= factor.z;
                }

                //DebugWrite("2");

                //now build the triangles using triTable
                for (int n=0; triTable[cubeIndex][n] != -1; n+=3) {
                    int index[3] = {triTable[cubeIndex][n+2], triTable[cubeIndex][n+1], triTable[cubeIndex][n]};
                    int cur;
                    for(int h=0; h < 3; h++) {    //copy vertices and normals into triangles array
                        cur = index[h];
                        triangles[numTriangles].p[h] = intVerts[cur];
                        triangles[numTriangles].norm[h] = grads[cur];
                        triangles[numTriangles].colors[h] = intColors[cur];
                        triangles[numTriangles].shaderTypes[h] = shaderTypes[cur];
                    }
                    numTriangles++;    //one more triangle has been added
                }
            }    //END OF FOR LOOP ON Z AXIS
        }    //END OF FOR LOOP ON Y AXIS
    }    //END OF FOR LOOP ON X AXIS

    //DebugWrite("3");

    //free all wasted space - no need since we delete this right after
//    MeshTriangle* retTriangles = new MeshTriangle[numTriangles];
//    for(int i=0; i < numTriangles; i++)
//        retTriangles[i] = triangles[i];
//    delete [] triangles;

    return triangles;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// THE END ///////////////////////////////////////////////////////////////

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
					/*(step 8)*/ 	Vector3D norm =
						cross(triangles[numTriangles].p[1] - triangles[numTriangles].p[0], triangles[numTriangles].p[2] -
							triangles[numTriangles].p[0]).unit();
                    triangles[numTriangles].norm[0] = norm;
                    triangles[numTriangles].norm[1] = norm;
                    triangles[numTriangles].norm[2] = norm;

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

