#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

#include "marchingCubes.h"

//#include <nanoflann.hpp>

//some force laws stored here for now
Vector3D r2_law(PointMass* p1, PointMass* p2) {
    Vector3D pos1 = p1->position;
    Vector3D pos2 = p2->position;
    Vector3D dir = pos1 - pos2;
    float distSqr = dir.norm2();
    dir = dir.unit();
    return dir / distSqr;
};

Vector3D r4_law(PointMass* p1, PointMass* p2) {
    Vector3D pos1 = p1->position;
    Vector3D pos2 = p2->position;
    Vector3D dir = pos1 - pos2;
    float distSqr = dir.norm2();
    dir = dir.unit();
    return dir / (distSqr * distSqr);
};

Vector3D cross_law(PointMass* p1, PointMass* p2) {
    Vector3D pos1 = p1->position;
    Vector3D pos2 = p2->position;
    Vector3D vel1 = pos1 - p1->last_position;
    Vector3D vel2 = pos2 - p2->last_position;
    Vector3D dir = pos1 - pos2;
    float distSqr = dir.norm2();
    dir = dir.unit();
    return cross(vel1, vel2);
};

using namespace std;

Cloth::~Cloth() {
    point_masses.clear();

    particleProperties.clear();
    particleColors.clear();
    
    delete[] cells;
}

void Cloth::spawnParticles(int count, Vector3D spawnPos, double spawnRadius, ParticleProperties properties) {
    int particleType = particleProperties.size();
    Vector3D curSpawnPos;
    
    for (int i = 0; i < count; i++) {
        curSpawnPos = spawnPos + Vector3D(
                    -spawnRadius + (rand() % 100 / 100.0) * (2 * spawnRadius),
                    -spawnRadius + (rand() % 100 / 100.0) * (2 * spawnRadius),
                    -spawnRadius + (rand() % 100 / 100.0) * (2 * spawnRadius));
        
        point_masses.push_back(PointMass(curSpawnPos, particleType));
    }
    
    particleProperties.push_back(properties);
    particleColors.push_back(properties.color);
}

void Cloth::simulate(double frames_per_sec, double simulation_steps,
    vector<Vector3D> external_accelerations,
    vector<CollisionObject*>* collision_objects) {
    
    //Particle force calculations-----------------------------
    double delta_t = 1.0f / frames_per_sec / simulation_steps;
    double distSqr;
    Vector3D dir;
    ParticleProperties* curProperties;
    PointMass* curMass;
    PointMass* otherMass;
    Vector3D forces;
    double collDist;

    //loop over all pairs of particles and compute forces
    for (int i = 0; i < point_masses.size(); i++) {
        curMass = &point_masses[i];
        curProperties = &particleProperties[curMass->particle_type];
        int cType = curMass->particle_type;
        for (int j = 0; j < point_masses.size(); j++) {
            if (i == j) continue;
            forces = Vector3D(0);
            otherMass = &point_masses[j];
            int oType = otherMass->particle_type;

            dir = otherMass->position - curMass->position;
            distSqr = dir.norm2();
            dir = dir.unit();

            if (distSqr <= 0.000005) continue; //cancel computation if distance too small to avoid explosions

            for (int k = 0; k < curProperties->force_laws.size(); k++) {
                Vector3D force = -curProperties->force_laws[k](curMass, otherMass);
                float factor = curProperties->strengths[k][oType];
                forces += factor * force;
            }
        }

        if (!curMass->pinned) {
            curMass->forces = forces;
            //curMass->forces += 0.3* cross(curMass->position, Vector3D(1.0, 0.0, 0.0));
            curMass->forces += Vector3D(0);
        }
    }

    for (auto p = point_masses.begin(); p != point_masses.end(); p++) {
        for (auto prim = collision_objects->begin(); prim != collision_objects->end(); prim++) {
            (*prim)->collide(*p);
        }
    }

    Vector3D curPos;
    Vector3D newPos;
    double delta_t_sqr = delta_t * delta_t;
    double dampingFactor = 1 - (0.01 / 100.0); // damping is arbitrary
    double delta = 0.1;
    Vector3D normal;
    PointMass* curParticle;

    //loop through and update particle positions via verlet integration
    for (int i = 0; i < point_masses.size(); i++) {
        curParticle = &point_masses[i];
        curProperties = &particleProperties[curParticle->particle_type];

        curPos = curParticle->position;
        newPos = curPos + dampingFactor * (curPos - curParticle->last_position) + ((curParticle->forces / curProperties->mass) * delta_t_sqr);

        // we cannot move farther than "delta" per timestep
        if ((curPos - newPos).norm2() > delta * delta) {
            newPos = curPos + (newPos - curPos).unit() * delta;
        }

        curParticle->position = newPos;
        curParticle->last_position = curPos;
    }

    /*for (auto p = point_masses.begin(); p != point_masses.end(); p++) {
        if (p->pinned) continue;
        self_collide(*p, simulation_steps);
    }*/

    
}

void Cloth::build_spatial_map() {
    for (const auto& entry : map) {
        delete(entry.second);
    }
    map.clear();

    // TODO (Part 4): Build a spatial map out of all of the point masses.

}

void Cloth::self_collide(PointMass& pm, double simulation_steps) {
    // TODO (Part 4): Handle self-collision for a given point mass.

}

float Cloth::hash_position(Vector3D pos) {
    // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.

    return 0.f;
}

// Marching cubes

void Cloth::initMarchingCubes(int numCells, double cellSize) {
    this->sideCellCount = numCells;
    this->cellSize = cellSize;
    
    double width = numCells * cellSize; // how big the box is
        
    int numPoints = numCells + 1;
    int yz = numPoints * numPoints;        //for little extra speed
    
    double min = -width / 2.0; // centers the entire grid
    
    totalVertexCount = numPoints * numPoints * numPoints;
    cells = new ScalarLoc[totalVertexCount];
    Vector3D bottomLeft = Vector3D(min);
    
    for (int i = 0; i < totalVertexCount; i++) {
        cells[i].pos = cellSize * Vector3D(i / yz, (i / numPoints) % numPoints, i % numPoints) + bottomLeft;
    }
}

MeshTriangle* Cloth::getMarchingCubeMesh(int& numTriangles) {
    double width = sideCellCount * cellSize; // how big the box is
    double min = -width / 2.0; // centers the entire grid
    int n = sideCellCount + 1;
    int yz = n * n;
    
    for (int i = 0; i < totalVertexCount; i++) {
        cells[i].value = 0;
        cells[i].color = Vector3D(0);
    }
    
    // "accurate" positions
    /*for (int i = 0; i < point_masses.size(); i += 3) {
        Vector3D pos = point_masses[i].position;

        int index = (int)(floor((pos.x - min) / c) * yz + floor((pos.y - min) / c) * n + floor((pos.z - min) / c));
        if (index < 0 || index >= numCells) continue;
        cells[index].value += 0.3;
    }*/

    // nick version
    /*Vector3D color;
    Vector3D color2;
    Vector3D color3;
    for (int i = 0; i < point_masses.size(); i++) {
        color = particleColors[point_masses[i].particle_type];
        for (int k = 0; k < point_masses.size(); k++) {
            color2 = particleColors[point_masses[k].particle_type];

            for (int l = 0; l < point_masses.size(); l++) {
                color3 = particleColors[point_masses[l].particle_type];

                for (int j = 0; j < point_masses.size(); j++) {
                    if (i == j || i == k || i == l || j == k || j == l || k == l) continue;
                    Vector3D pos = 0.25 * (point_masses[i].position + point_masses[j].position + point_masses[k].position + point_masses[l].position);
                    int index = (int)(floor((pos.x - min) / c) * yz + floor((pos.y - min) / c) * n + floor((pos.z - min) / c));
                    if (index < 0 || index >= numCells) continue;
                    
                    Vector3D newColor = 0.25 * (color + color2 + color3 + particleColors[point_masses[j].particle_type]);

                    cells[index].color = cells[index].color * (cells[index].value / (0.01 + cells[index].value)) + newColor * (0.01 / (0.01 + cells[index].value));
                    cells[index].value += 0.01;
                }
            }
        }
    }*/
    
    double particleRadius;
    int gridRadius;
    
    Vector3D offset;
    Vector3D particlePos;
    int startIndex;
    int curIndex;
    double dist;
    
    double newVal;
    double particleStrength;
    Vector3D particleColor;
    
    int cellPos[3];
    
    for (int i = 0; i < point_masses.size(); i++) {
        particleRadius = particleProperties[point_masses[i].particle_type].radius;
        gridRadius = ceil(particleRadius / cellSize);
        particlePos = point_masses[i].position;
        
        cellPos[0] = floor((particlePos.x - min) / cellSize);
        if (cellPos[0] < 0 || cellPos[0] >= n) continue;

        cellPos[1] = floor((particlePos.y - min) / cellSize);
        if (cellPos[1] < 0 || cellPos[1] >= n) continue;

        cellPos[2] = floor((particlePos.z - min) / cellSize);
        if (cellPos[2] < 0 || cellPos[2] >= n) continue;
                    
        startIndex = (int)(cellPos[0] * yz + cellPos[1] * n + cellPos[2]);

        offset = particlePos - cells[startIndex].pos;
        particleColor = particleColors[point_masses[i].particle_type];
        
        //TODO: can probably be optimized!!!!
        for (int dx = -gridRadius; dx <= gridRadius; dx++) {
            if (dx < -cellPos[0] || dx >= n - cellPos[0]) continue;
            
            for (int dy = -gridRadius; dy <= gridRadius; dy++) {
                if (dy < -cellPos[1] || dy >= n - cellPos[1]) continue;

                curIndex = startIndex + dx * yz + dy * n - gridRadius;
                for (int dz = -gridRadius; dz <= gridRadius; dz++) {
                    if (dz >= -cellPos[2] && dz < n - cellPos[2]) {
                        dist = (cellSize * Vector3D(dx, dy, dz) - offset).norm();
                        
                        if (dist <= particleRadius) {
                            particleStrength = cellSize + 0.11;
                        } else {
                            particleStrength = (cellSize / 100) / (dist * dist * dist);
                        }
                        
                        newVal = cells[curIndex].value + particleStrength;
                        cells[curIndex].color = (cells[curIndex].value / newVal) * cells[curIndex].color + (particleStrength / newVal) * particleColor;
                        cells[curIndex].value = newVal;
                    }
                    curIndex++;
                }
            }
        }
        
        
//        for (int j = 0; j < point_masses.size(); j += 3) {
//            if (i == j) continue;
//            Vector3D pos = 0.5 * (point_masses[i].position + point_masses[j].position);
//            int index = (int)(floor((pos.x - min) / c) * yz + floor((pos.y - min) / c) * n + floor((pos.z - min) / c));
//            if (index < 0 || index >= numCells) continue;
//
//            cells[index].value += 0.01;
//        }
    }
    
    double minValue = cellSize / 2;
    MeshTriangle* triangles = MarchingCubes(sideCellCount, sideCellCount, sideCellCount, cellSize, cellSize, cellSize, minValue, cells, numTriangles);
    
//    MeshTriangle* triangles = MarchingCubesCross(size, size, size, c, cells, numTriangles);
    
//    MeshTriangle* triangles = MCRecFind(size, size, size, c, c, c, 0.1f, cells, numTriangles);
    return triangles;
}


///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
    PointMass* pm = &point_masses[0];
    for (int i = 0; i < point_masses.size(); i++) {
        pm->position = pm->start_position;
        pm->last_position = pm->start_position;
        pm++;
    }
}
