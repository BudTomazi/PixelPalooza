#ifndef CLOTH_H
#define CLOTH_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "clothMesh.h"
#include "collision/collisionObject.h"
#include "marchingCubes.h"
#include "particle.h"

using namespace CGL;
using namespace std;

enum e_orientation { HORIZONTAL = 0, VERTICAL = 1 };

Vector3D r2_law(Particle* p1, Particle* p2);
Vector3D r4_law(Particle* p1, Particle* p2);
Vector3D cross_law(Particle* p1, Particle* p2);
Vector3D fire_force(Particle* p1, Particle* p2);

struct ClothParameters {
    ClothParameters() {}
    ClothParameters(bool enable_structural_constraints,
        bool enable_shearing_constraints,
        bool enable_bending_constraints, double damping,
        double density, double ks)
        : enable_structural_constraints(enable_structural_constraints),
        enable_shearing_constraints(enable_shearing_constraints),
        enable_bending_constraints(enable_bending_constraints),
        damping(damping), density(density), ks(ks) {}
    ~ClothParameters() {}

    // Global simulation parameters

    bool enable_structural_constraints;
    bool enable_shearing_constraints;
    bool enable_bending_constraints;

    double damping;

    // Mass-spring parameters
    double density;
    double ks;
};

struct Cloth {
    Cloth() {}
    ~Cloth();

    void spawnParticles(int count, Vector3D spawnPos, Vector3D spawnExtents, ParticleProperties properties);

    void simulate(vector<CollisionObject*>* collision_objects);

    void reset();

    void build_spatial_map();
    void self_collide(Particle& pm);
    float hash_position(Vector3D pos);
    
    // Cloth components
    vector<Particle> particles;
    double damping;
    
    bool isInSphere;
    double sphereRad;

    // Particle properties
    vector<ParticleProperties> particleProperties;
    vector<Vector3D> particleColors;
    
    double simulation_steps;
    double frames_per_sec;
    
    // marching cubes
    ScalarLoc* cells;
    int totalVertexCount;
    int sideCellCount[3];
    double cellSize;
    Vector3D borderDist;
    Vector3D physicsBorder;
    vector<Vector3D> planeLocs;
    vector<Vector3D> planeNorms;
    
    void initMarchingCubes(int numCellsX, int numCellsY, int numCellsZ, double cellSize, double physicsBuffer, bool noTop);
    MeshTriangle* getMarchingCubeMesh(int& numTriangles);
    MeshTriangle* getTrianglesMesh(int& numTriangles);

    // Spatial hashing
    unordered_map<float, vector<Particle*>*> map;
};

#endif /* CLOTH_H */
