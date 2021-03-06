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

Vector3D fire_force(Particle* p1, Particle* p2) {
    Vector3D pos1 = p1->position;
    Vector3D pos2 = p2->position;
    Vector3D dir = pos1 - pos2;
    float dist = dir.norm();
    dir = Vector3D(0.0, -1.0, 0.0);
    dir[0] = 0.2 * p1->position.x;
    dir[2] = 0.2 * p1->position.z;
    //float factor = (3.0 + p1->position.y)/(10.0);

    //dir[0] *= factor;
    //dir[2] *= factor;
    //dir[1] *= (1.2 - factor);
    //new commen
    //if(rand()%4 == 0) return Vector3D(0.0);
    return dir/dist;
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

double smoothingKernel(Vector3D x1, Vector3D x2, double h2, double h9) {
    double hr = (h2 - (x2 - x1).norm2());
    if (hr < 0) return 0;
    return 1.5 * hr * hr * hr / h9;
}

using namespace std;

Cloth::~Cloth() {
    particles.clear();

    particleProperties.clear();
    particleColors.clear();

    delete[] cells;
}

void Cloth::spawnParticles(int count, Vector3D spawnPos, Vector3D spawnExtents, ParticleProperties properties) {
    int particleType = particleProperties.size();
    Vector3D curSpawnPos;
    
    Vector3D velOffset = -properties.velocity * (1.0f / frames_per_sec / simulation_steps);

    int sideCount = round(pow(count, (1.0 / 3.0)));
    Vector3D spawnOffset;
    double r = 0.02;
    for (int x = 0; x < sideCount; x++) {
        spawnOffset.x = -spawnExtents.x + (2 * spawnExtents.x) * (x / (double)sideCount);
        for (int y = 0; y < sideCount; y++) {
            spawnOffset.y = -spawnExtents.y + (2 * spawnExtents.y) * (y / (double)sideCount);
            for (int z = 0; z < sideCount; z++) {
                spawnOffset.z = -spawnExtents.z + (2 * spawnExtents.z) * (z / (double)sideCount);
                curSpawnPos = spawnPos + spawnOffset +
                    Vector3D(-r * (rand() % 100 / 100.0) * (2 * r),
                             -r + (rand() % 100 / 100.0) * (2 * r),
                             -r + (rand() % 100 / 100.0) * (2 * r));

                Particle newParticle = Particle(curSpawnPos, particleType);
                newParticle.last_position += velOffset;
                particles.push_back(newParticle);
            }
        }
    }

    particleProperties.push_back(properties);
    particleColors.push_back(properties.color);
}

void Cloth::simulate(vector<CollisionObject*>* collision_objects) {
    //Particle force calculations-----------------------------
    double delta_t = 1.0f / frames_per_sec / simulation_steps;
    double distSqr;
    ParticleProperties* curProperties;
    Particle* curParticle;
    Particle* otherParticle;
    Vector3D forces;
    double collDist;
    
    Vector3D dir;

    //loop over all pairs of particles and compute forces
    build_spatial_map();
    for (auto p = particles.begin(); p != particles.end(); p++) {
        if (particleProperties[p->particle_type].pinned) continue;
        self_collide(*p);
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

                    if (distSqr <= 0.005) continue; //cancel computation if distance too small to avoid explosions

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

                    if (distSqr <= 0.005) continue; //cancel computation if distance too small to avoid explosions

                    Vector3D force = -curProperties->force_laws[k](curParticle, p);
                    float factor = curProperties->strengths[k][oType];
                    forces += factor * force;
                }
            }
        }

        forces += curProperties->external_forces;

        if (!curProperties->pinned) {
            curParticle->forces = forces;
            //curMass->forces += 0.3* cross(curMass->position, Vector3D(1.0, 0.0, 0.0));
        }

    }

    if (!isInSphere) {
        for (int i = 0; i < planeNorms.size(); i++) {
            for (int j = 0; j < particles.size(); j++) {
                PlaneCollision(planeLocs[i], planeNorms[i], particles[j]);
            }
        }
    } else {
        Vector3D posDir;

        for (Particle& curParticle : particles) {
            posDir = curParticle.position.unit();

            if (curParticle.position.norm() >= sphereRad) {
                
                curParticle.forces += dot(curParticle.forces, posDir) * posDir;
                curParticle.last_position -= dot(curParticle.last_position - curParticle.position, posDir) * posDir;
            }
        }
    }

    for (auto p = particles.begin(); p != particles.end(); p++) {
        for (auto prim = collision_objects->begin(); prim != collision_objects->end(); prim++) {
            (*prim)->collide(*p);
        }
    }

    Vector3D curPos;
    Vector3D newPos;
    double delta_t_sqr = delta_t * delta_t;
    double dampingFactor = 1 - (damping / 100.0);
    double delta = 0.1;
    Vector3D normal;

    //cerr << "help1\n";

    //cerr << "help2\n";

    //loop through and update particle positions via verlet integration
    for (int i = 0; i < particles.size(); i++) {
        curParticle = &particles[i];

        curProperties = &particleProperties[curParticle->particle_type];
        if (curProperties->pinned) continue;
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
        self_collide(*p);
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

void Cloth::self_collide(Particle& pm) {
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
            int last_type = pm.particle_type;
            pm.particle_type = curProperties->collision_transformations[p->particle_type];
            p->particle_type = otherProperties->collision_transformations[last_type];
        }
    }


    //cerr << "end self collide\n";
    if (counter == 0 || !curProperties->particle_collisions) {
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
    double w = 2 * physicsBorder.x / num_width_points;
    double h = 2 * physicsBorder.y / num_width_points;
    double z = 2 * physicsBorder.z / num_height_points;


    int x_coord;
    int y_coord;
    int z_coord;

    x_coord = (int)((pos.x + physicsBorder.x) / w);
    y_coord = (int)((pos.y + physicsBorder.y) / h);
    z_coord = (int)((pos.z + physicsBorder.z) / z);


    float output = x_coord + num_width_points * y_coord + num_width_points * num_height_points * z_coord;
   return output;
}

// Marching cubes

void Cloth::initMarchingCubes(int numCellsX, int numCellsY, int numCellsZ, double cellSize, double physicsBuffer, bool noTop) {
    sideCellCount[0] = numCellsX;
    sideCellCount[1] = numCellsY;
    sideCellCount[2] = numCellsZ;
    this->cellSize = cellSize;

    borderDist = cellSize * Vector3D(numCellsX, numCellsY, numCellsZ) / 2;

    int numPoints[3];
    numPoints[0] = sideCellCount[0]+1;
    numPoints[1] = sideCellCount[1]+1;
    numPoints[2] = sideCellCount[2]+1;

    int yz = numPoints[1] * numPoints[2];        //for little extra speed

    physicsBorder = borderDist - Vector3D(physicsBuffer);
    if (!isInSphere) {
        planeLocs.emplace_back(Vector3D(physicsBorder.x, 0.0, 0.0));
        planeLocs.emplace_back(Vector3D(-physicsBorder.x, 0.0, 0.0));
        if (!noTop) planeLocs.emplace_back(Vector3D(0.0, physicsBorder.y, 0.0));
        planeLocs.emplace_back(Vector3D(0.0, -physicsBorder.y, 0.0));
        planeLocs.emplace_back(Vector3D(0.0, 0.0, physicsBorder.z));
        planeLocs.emplace_back(Vector3D(0.0, 0.0, -physicsBorder.z));

        planeNorms.emplace_back(Vector3D(-1.0, 0.0, 0.0));
        planeNorms.emplace_back(Vector3D(1.0, 0.0, 0.0));
        if (!noTop) planeNorms.emplace_back(Vector3D(0.0, -1.0, 0.0));
        planeNorms.emplace_back(Vector3D(0.0, 1.0, 0.0));
        planeNorms.emplace_back(Vector3D(0.0, 0.0, -1.0));
        planeNorms.emplace_back(Vector3D(0.0, 0.0, 1.0));
    } else {
        sphereRad = min(physicsBorder.x, min(physicsBorder.y, physicsBorder.z));
    }

    totalVertexCount = numPoints[0] * numPoints[1] * numPoints[2];
    cells = new ScalarLoc[totalVertexCount];
    Vector3D bottomLeft = Vector3D(-borderDist);

    for (int i = 0; i < totalVertexCount; i++) {
        cells[i].pos = cellSize * Vector3D(i / yz, (i / numPoints[2]) % numPoints[1], i % numPoints[2]) + bottomLeft;
    }
}

MeshTriangle* Cloth::getTrianglesMesh(int& numTriangles) {
    numTriangles = particles.size();
    MeshTriangle* tempTriangles = new MeshTriangle[numTriangles];
    Vector3D centroid = Vector3D(0);
    
    for (Particle& curParticle : particles) {
        centroid += curParticle.position;
    }
    
    centroid /= numTriangles;
    
    for (int i = 0; i < numTriangles; i++) {
        ParticleProperties* curParticleProperties = &particleProperties[particles[i].particle_type];
        Vector3D normal = (centroid - particles[i].position).unit();
        Vector3D orthogonal = Vector3D(1, 1, -(normal.x + normal.y) / normal.z).unit();
        
        tempTriangles[i].p[0] = particles[i].position;
        tempTriangles[i].p[1] = particles[i].position + orthogonal * curParticleProperties->radius * 2;
        tempTriangles[i].p[2] = particles[i].position + cross(orthogonal, normal) * curParticleProperties->radius * 2;

        tempTriangles[i].norm[0] = normal;
        tempTriangles[i].norm[1] = normal;
        tempTriangles[i].norm[2] = normal;

        tempTriangles[i].colors[0] = curParticleProperties->color;
        tempTriangles[i].colors[1] = curParticleProperties->color;
        tempTriangles[i].colors[2] = curParticleProperties->color;

        tempTriangles[i].shaderTypes[0] = curParticleProperties->shaderType;
        tempTriangles[i].shaderTypes[1] = curParticleProperties->shaderType;
        tempTriangles[i].shaderTypes[2] = curParticleProperties->shaderType;
    }
    
    return tempTriangles;
}

MeshTriangle* Cloth::getMarchingCubeMesh(int& numTriangles) {
    int numX = sideCellCount[0] + 1;
    int numY = sideCellCount[1] + 1;
    int numZ = sideCellCount[2] + 1;
    int yz = numY * numZ;

    int numParticleTypes = particleProperties.size();
    for (int i = 0; i < totalVertexCount; i++) {
        cells[i].value = 0;
        cells[i].color = Vector3D(0);
        cells[i].shaderType = -1;
    }

    if (useDensity) {
        build_spatial_map();
        ParticleProperties* curProperties;
        for (Particle& curParticle : particles) {
            curParticle.density = 0;
            vector<Particle*>* localCandidates = map[hash_position(curParticle.position)];

            curProperties = &particleProperties[curParticle.particle_type];
            for (Particle* p : *localCandidates) {
                if (p->particle_type != curParticle.particle_type) continue;
                
                double h = curProperties->radius * 2;
                double h2 = h * h;
                double h9 = h2 * h2 * h2 * h2 * h;
                curParticle.density += curProperties->mass * smoothingKernel(curParticle.position, p->position, h2, h9);
            }
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
    Vector3D curParticleColor;

    double delta_t = 1.0 / frames_per_sec / simulation_steps;
    
    int cellPos[3];

    for (int i = 0; i < particles.size(); i++) {
        curProperties = &particleProperties[particles[i].particle_type];

        gridRadius = ceil(curProperties->radius / cellSize) * 2;
        particlePos = particles[i].position;

        cellPos[0] = floor((particlePos.x + borderDist.x) / cellSize);
        if (cellPos[0] < 0 || cellPos[0] >= numX) continue;

        cellPos[1] = floor((particlePos.y + borderDist.y) / cellSize);
        if (cellPos[1] < 0 || cellPos[1] >= numY) continue;

        cellPos[2] = floor((particlePos.z + borderDist.z) / cellSize);
        if (cellPos[2] < 0 || cellPos[2] >= numZ) continue;

        startIndex = (int)(cellPos[0] * yz + cellPos[1] * numZ + cellPos[2]);

        offset = particlePos - cells[startIndex].pos;
        curParticleColor = curProperties->color;
        if (curProperties->velocityColor.norm2() > 0) {
            curParticleColor += curProperties->velocityColor * particles[i].velocity(delta_t).norm2();
        }

        //TODO: can probably be optimized!!!!
        double h = curProperties->radius * 2;
        double h2 = h * h;
        double h9 = h2 * h2 * h;

        for (int dx = -gridRadius; dx <= gridRadius; dx++) {
            if (dx < -cellPos[0] || dx >= numX - cellPos[0]) continue;

            for (int dy = -gridRadius; dy <= gridRadius; dy++) {
                if (dy < -cellPos[1] || dy >= numY - cellPos[1]) continue;

                curIndex = startIndex + dx * yz + dy * numZ - gridRadius;
                for (int dz = -gridRadius; dz <= gridRadius; dz++) {
                    if (dz >= -cellPos[2] && dz < numZ - cellPos[2]) {
                        if (useDensity) {
                            particleStrength = smoothingKernel(cells[curIndex].pos, particlePos, h2, h9) * curProperties->mass / particles[i].density;
                        } else {
                            particleStrength = smoothingKernel(cells[curIndex].pos, particlePos, h2, h9);
                        }

//                        dist = (cellSize * Vector3D(dx, dy, dz) - offset).norm();
//
//                        if (dist <= curProperties->radius) {
//                            particleStrength = cellSize + 0.11;
//                        } else {
//                            particleStrength = (cellSize / 100) / (dist * dist * dist);
//                        }
                        
                        if (particleStrength > 0.000001) {
                            newVal = cells[curIndex].value + particleStrength;
                            cells[curIndex].color = (cells[curIndex].value / newVal) * cells[curIndex].color + (particleStrength / newVal) * curParticleColor;
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


        if (curProperties->particleAveragingFactor > 0) {
            int curParticleType = particles[i].particle_type;
            curParticleColor *= curProperties->particleAveragingBrightness;
            
            for (Particle& other : particles) {
                if (other.particle_type == curParticleType) {
                    if (curProperties->particleAveragingDist > 0 && (particles[i].position - other.position).norm2() > curProperties->particleAveragingDist) {
                        continue;
                    }
                    
                    Vector3D pos = 0.5 * (particles[i].position + other.position);

                    cellPos[0] = floor((pos.x + borderDist.x) / cellSize);
                    if (cellPos[0] < 0 || cellPos[0] >= numX) continue;

                    cellPos[1] = floor((pos.y + borderDist.y) / cellSize);
                    if (cellPos[1] < 0 || cellPos[1] >= numY) continue;

                    cellPos[2] = floor((pos.z + borderDist.z) / cellSize);
                    if (cellPos[2] < 0 || cellPos[2] >= numZ) continue;

                    curIndex = (int)(cellPos[0] * yz + cellPos[1] * numZ + cellPos[2]);

                    particleStrength = curProperties->particleAveragingFactor;

                    newVal = cells[curIndex].value + particleStrength;
                    cells[curIndex].color = (cells[curIndex].value / newVal) * cells[curIndex].color + (particleStrength / newVal) * curParticleColor;
                    cells[curIndex].value = newVal;
                    if (curProperties->shaderType > cells[curIndex].shaderType) {
                        cells[curIndex].shaderType = curProperties->shaderType;
                    }
                }
            }
        }
    }

    double minValue = 0.02;
    MeshTriangle* triangles = MarchingCubes(sideCellCount[0], sideCellCount[1], sideCellCount[2], cellSize, cellSize, cellSize, minValue, cells, numTriangles);

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
        Vector3D velOffset = -particleProperties[particles[i].particle_type].velocity * (1.0f / frames_per_sec / simulation_steps);
        
        particle->position = particle->start_position;
        particle->last_position = particle->start_position + velOffset;
        particle->particle_type = particle->start_type;
        particle++;
    }
}
