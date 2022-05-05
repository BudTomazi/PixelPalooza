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
Vector3D r2_law(Particle* p1, Particle* p2) {
    Vector3D dir = p1->position - p2->position;
    double distSqr = dir.norm2();
    dir = dir / sqrt(distSqr);
    return dir / distSqr;
};

Vector3D r4_law(Particle* p1, Particle* p2) {
    Vector3D dir = p1->position - p2->position;
    float distSqr = dir.norm2();
    dir = dir / sqrt(distSqr);
    return dir / (distSqr * distSqr);
};

Vector3D cross_law(Particle* p1, Particle* p2) {
    Vector3D pos1 = p1->position;
    Vector3D pos2 = p2->position;
    Vector3D vel1 = pos1 - p1->last_position;
    Vector3D vel2 = pos2 - p2->last_position;
    Vector3D dir = pos1 - pos2;
    float distSqr = dir.norm2();
    return cross(vel1, vel2);
};

void PlaneCollision(Vector3D planeLoc, Vector3D planeNorm, Particle& p) {
    double last_side, curr_side;
    last_side = dot(p.last_position - planeLoc, planeNorm);
    curr_side = dot(p.position - planeLoc, planeNorm);
    //abs(curr_side - last_side) > abs(curr_side)

    bool a = last_side > 0;
    bool b = curr_side > 0;
    if (a != b) {
        Vector3D diff = (p.position - p.last_position);
        p.position -= planeNorm * dot(diff, planeNorm);
    }
//    if (abs(curr_side) < 0.1) {
//        if (dot(p.forces, planeNorm) < 0.0) {
//            //p.forces -= dot(p.forces, planeNorm);
//        }
//    }
}

const double scalingFactor = 3000.0 / (64 * 3.14159);
double smoothingKernel(Vector3D x1, Vector3D x2, double h2, double h9) {
    double hr = (h2 - (x2 - x1).norm2());
    if (hr < 0) return 0;
    return scalingFactor * hr * hr * hr / h9;
}

using namespace std;

Cloth::~Cloth() {
    particles.clear();

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
        
        particles.push_back(Particle(curSpawnPos, particleType));
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
    Particle* curParticle;
    Particle* otherParticle;
    Vector3D forces;
    double collDist;

    //loop over all pairs of particles and compute forces
    build_spatial_map();
    for (auto p = particles.begin(); p != particles.end(); p++) {
        if (p->pinned) continue;
        self_collide(*p, simulation_steps);
    }
    build_spatial_map();
    for (int i = 0; i < particles.size(); i++) {
        curParticle = &particles[i];
        curProperties = &particleProperties[curParticle->particle_type];
        int cType = curParticle->particle_type;
        forces = Vector3D(0);
        
        vector<Particle*>* localCandidates = map[hash_position(curParticle->position)];
        
        for (int k = 0; k < curProperties->force_laws.size(); k++) {
            if (!curProperties->localized[k]) {
                for (int j = 0; j < particles.size(); j++) {
                    if (i == j) continue;
                    otherParticle = &particles[j];
                    int oType = otherParticle->particle_type;

                    dir = otherParticle->position - curParticle->position;
                    distSqr = dir.norm2();

                    if (distSqr <= 0.000005) continue; //cancel computation if distance too small to avoid explosions

                    Vector3D force = -curProperties->force_laws[k](curParticle, otherParticle);
                    float factor = curProperties->strengths[k][oType];
                    forces += factor * force;
                }
            }
            else {
                for (Particle* p : *localCandidates) {
                    if (p == curParticle) continue;
                    
                    int oType = p->particle_type;
                    dir = p->position - curParticle->position;
                    distSqr = dir.norm2();

                    if (distSqr <= 0.000005) continue; //cancel computation if distance too small to avoid explosions

                    Vector3D force = -curProperties->force_laws[k](curParticle, p);
                    float factor = curProperties->strengths[k][oType];
//                    forces += factor * force;
                }
            }
        }

        if (curProperties->external_forces) {
            forces += Vector3D(0.0, -100.0, 0.0);
        }
        if (cType == 2) {
            forces += Vector3D(0.0, 5.0, 0.0);
        }

        if (!curProperties->pinned) {
            curParticle->forces = forces;
            //curMass->forces += 0.3* cross(curMass->position, Vector3D(1.0, 0.0, 0.0));
        }

    }

    for (auto p = particles.begin(); p != particles.end(); p++) {
        for (auto prim = collision_objects->begin(); prim != collision_objects->end(); prim++) {
            (*prim)->collide(*p);
        }
    }
    
    for (int i = 0; i < planeNorms.size(); i++) {
        for (int j = 0; j < particles.size(); j++) {
            PlaneCollision(planeLocs[i], planeNorms[i], particles[j]);
        }
    }
    
    Vector3D curPos;
    Vector3D newPos;
    double delta_t_sqr = delta_t * delta_t;
    double dampingFactor = 1 - (0.01 / 100.0); // damping is arbitrary
    double delta = 0.1;
    Vector3D normal;

    //cerr << "help1\n";
    
    //cerr << "help2\n";

    //loop through and update particle positions via verlet integration
    for (int i = 0; i < particles.size(); i++) {
        curParticle = &particles[i];
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

    
    //build_spatial_map();
    /*for (auto p = particles.begin(); p != particles.end(); p++) {
        if (p->pinned) continue;
        self_collide(*p, simulation_steps);
    }*/

    
}

void Cloth::build_spatial_map() {
    //cerr << "help3\n";
    for (const auto& entry : map) {
        delete(entry.second);
    }
    map.clear();

    //cerr << "help4\n";

    // TODO (Part 4): Build a spatial map out of all of the point masses.
    for (Particle& curMass : particles) {
        float hashed_position = hash_position(curMass.position);
        if (map.find(hashed_position) == map.end()) {
            map[hashed_position] = new vector<Particle*>;
        }
        map[hashed_position]->push_back(&curMass);
    }
    //cerr << "help5\n";
    // TODO (Part 4): Build a spatial map out of all of the point masses.

}

void Cloth::self_collide(Particle& pm, double simulation_steps) {
    vector<Particle*>* candidates = map[hash_position(pm.position)];

    //cerr << "start self collide\n";
    Vector3D correction = Vector3D(0.0, 0.0, 0.0);
    double counter = 0.0;

    ParticleProperties* curProperties = &particleProperties[pm.particle_type];
    float r1 = curProperties->collRadius;

    for (Particle* p : *candidates) {
        ParticleProperties* otherProperties = &particleProperties[p->particle_type];
        float r2 = otherProperties->collRadius;
        if (p != &pm && (pm.position - p->position).norm() < (r1 + r2)) {
            counter += 1;
            correction += (pm.position - p->position).unit() * (r1 + r2) + p->position - pm.position;
            pm.particle_type = curProperties->collision_transformations[p->particle_type];
        }
    }


    //cerr << "end self collide\n";
    if (counter == 0) {
        correction = Vector3D(0.0, 0.0, 0.0);
    }
    else {
        correction = correction / counter;
    }

    pm.position = (correction / simulation_steps) + pm.position;
    //pm.last_position = pm.position;

}

float Cloth::hash_position(Vector3D pos) {
    // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.

    double num_width_points = 4;
    double num_height_points = 4;
    double w = 2 * physicsBorder / num_width_points;
    double h = 2 * physicsBorder / num_width_points;
    double z = 2 * physicsBorder / num_height_points;
    

    int x_coord;
    int y_coord;
    int z_coord;

    x_coord = (int)((pos.x + physicsBorder) / w);
    y_coord = (int)((pos.y + physicsBorder) / h);
    z_coord = (int)((pos.z + physicsBorder) / z);
 

    float output = x_coord + num_width_points * y_coord + num_width_points * num_height_points * z_coord;
   return output; 
}

// Marching cubes

void Cloth::initMarchingCubes(int numCells, double cellSize) {
    this->sideCellCount = numCells;
    this->cellSize = cellSize;
    
    double width = numCells * cellSize; // how big the box is
        
    int numPoints = numCells + 1;
    int yz = numPoints * numPoints;        //for little extra speed
    
    borderDist = width / 2.0; // centers the entire grid
    physicsBorder = borderDist / 3.0;// - 20 * cellSize;
    planeLocs.emplace_back(Vector3D(physicsBorder, 0.0, 0.0));
    planeLocs.emplace_back(Vector3D(-physicsBorder, 0.0, 0.0));
    //planeLocs.emplace_back(Vector3D(0.0, physicsBorder, 0.0));
    planeLocs.emplace_back(Vector3D(0.0, -physicsBorder, 0.0));
    planeLocs.emplace_back(Vector3D(0.0, 0.0, physicsBorder));
    planeLocs.emplace_back(Vector3D(0.0, 0.0, -physicsBorder));

    planeNorms.emplace_back(Vector3D(-1.0, 0.0, 0.0));
    planeNorms.emplace_back(Vector3D(1.0, 0.0, 0.0));
    //planeNorms.emplace_back(Vector3D(0.0, -1.0, 0.0));
    planeNorms.emplace_back(Vector3D(0.0, 1.0, 0.0));
    planeNorms.emplace_back(Vector3D(0.0, 0.0, -1.0));
    planeNorms.emplace_back(Vector3D(0.0, 0.0, 1.0));
    
    totalVertexCount = numPoints * numPoints * numPoints;
    cells = new ScalarLoc[totalVertexCount];
    Vector3D bottomLeft = Vector3D(-borderDist);
    
    for (int i = 0; i < totalVertexCount; i++) {
        cells[i].pos = cellSize * Vector3D(i / yz, (i / numPoints) % numPoints, i % numPoints) + bottomLeft;
    }
}

MeshTriangle* Cloth::getMarchingCubeMesh(int& numTriangles) {
    double width = sideCellCount * cellSize; // how big the box is
    int n = sideCellCount + 1;
    int yz = n * n;
    
    int numParticleTypes = particleProperties.size();
    for (int i = 0; i < totalVertexCount; i++) {
        cells[i].value = 0;
        cells[i].color = Vector3D(0);
        cells[i].shaderType = -1;
    }
    
    build_spatial_map();
    ParticleProperties* otherProperties;
    for (Particle& curParticle : particles) {
        curParticle.density = 0;
        vector<Particle*>* localCandidates = map[hash_position(curParticle.position)];

        for (Particle* p : *localCandidates) {
//            if (p == curParticle) continue;
            otherProperties = &particleProperties[p->particle_type];
            double h2 = otherProperties->radius * otherProperties->radius;
            double h9 = h2 * h2 * h2 * h2 * otherProperties->radius;
            curParticle.density += otherProperties->mass * smoothingKernel(curParticle.position, p->position, h2, h9);
        }
    }

    // nick version
    /*Vector3D color;
    Vector3D color2;
    Vector3D color3;
    for (int i = 0; i < particles.size(); i++) {
        color = particleColors[particles[i].particle_type];
        for (int k = 0; k < particles.size(); k++) {
            color2 = particleColors[particles[k].particle_type];

            for (int l = 0; l < particles.size(); l++) {
                color3 = particleColors[particles[l].particle_type];

                for (int j = 0; j < particles.size(); j++) {
                    if (i == j || i == k || i == l || j == k || j == l || k == l) continue;
                    Vector3D pos = 0.25 * (particles[i].position + particles[j].position + particles[k].position + particles[l].position);
                    int index = (int)(floor((pos.x - min) / c) * yz + floor((pos.y - min) / c) * n + floor((pos.z - min) / c));
                    if (index < 0 || index >= numCells) continue;
                    
                    Vector3D newColor = 0.25 * (color + color2 + color3 + particleColors[particles[j].particle_type]);

                    cells[index].color = cells[index].color * (cells[index].value / (0.01 + cells[index].value)) + newColor * (0.01 / (0.01 + cells[index].value));
                    cells[index].value += 0.01;
                }
            }
        }
    }*/
    
    int gridRadius;
    
    Vector3D offset;
    Vector3D particlePos;
    int startIndex;
    int curIndex;
    double dist;
    
    double newVal;
    double particleStrength;
    ParticleProperties* curProperties;
    
    int cellPos[3];
    
    for (int i = 0; i < particles.size(); i++) {
        curProperties = &particleProperties[particles[i].particle_type];
        
        gridRadius = ceil(curProperties->radius / cellSize);
        particlePos = particles[i].position;
        
        cellPos[0] = floor((particlePos.x + borderDist) / cellSize);
        if (cellPos[0] < 0 || cellPos[0] >= n) continue;

        cellPos[1] = floor((particlePos.y + borderDist) / cellSize);
        if (cellPos[1] < 0 || cellPos[1] >= n) continue;

        cellPos[2] = floor((particlePos.z + borderDist) / cellSize);
        if (cellPos[2] < 0 || cellPos[2] >= n) continue;
                    
        startIndex = (int)(cellPos[0] * yz + cellPos[1] * n + cellPos[2]);

        offset = particlePos - cells[startIndex].pos;
        
        //TODO: can probably be optimized!!!!
        double h2 = curProperties->radius;
        double h9 = h2 * h2 * h2 * h2 * curProperties->radius;
        
        for (int dx = -gridRadius; dx <= gridRadius; dx++) {
            if (dx < -cellPos[0] || dx >= n - cellPos[0]) continue;
            
            for (int dy = -gridRadius; dy <= gridRadius; dy++) {
                if (dy < -cellPos[1] || dy >= n - cellPos[1]) continue;

                curIndex = startIndex + dx * yz + dy * n - gridRadius;
                for (int dz = -gridRadius; dz <= gridRadius; dz++) {
                    if (dz >= -cellPos[2] && dz < n - cellPos[2]) {
                        particleStrength = curProperties->mass * smoothingKernel(cellSize * Vector3D(dx, dy, dz), offset, h2, h9) / particles[i].density;
//                        dist = (cellSize * Vector3D(dx, dy, dz) - offset).norm();
//
//                        if (dist <= curProperties->radius) {
//                            particleStrength = cellSize + 0.11;
//                        } else {
//                            particleStrength = (cellSize / 100) / (dist * dist * dist);
//                        }
                        if (particleStrength > 0) {
                            newVal = cells[curIndex].value + particleStrength;
                            cells[curIndex].color = (cells[curIndex].value / newVal) * cells[curIndex].color + (particleStrength / newVal) * curProperties->color;
                            cells[curIndex].value = newVal;
                            if (curProperties->shaderType > cells[curIndex].shaderType) {
                                cells[curIndex].shaderType = curProperties->shaderType;
                            }
                        }
                    }
                    curIndex++;
                }
            }
        }
    }
    
    double minValue = 0.02;
    MeshTriangle* triangles = MarchingCubes(sideCellCount, sideCellCount, sideCellCount, cellSize, cellSize, cellSize, minValue, cells, numTriangles);
    
//    MeshTriangle* triangles = MarchingCubesCross(size, size, size, c, cells, numTriangles);
    
//    MeshTriangle* triangles = MCRecFind(size, size, size, c, c, c, 0.1f, cells, numTriangles);
    return triangles;
}


///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
    Particle* particle = &particles[0];
    for (int i = 0; i < particles.size(); i++) {
        particle->position = particle->start_position;
        particle->last_position = particle->start_position;
        particle++;
    }
}
