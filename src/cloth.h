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

    void spawnParticles(int count, Vector3D spawnPos, double spawnRadius, ParticleProperties properties);

    void simulate(double frames_per_sec, double simulation_steps,
        vector<Vector3D> external_accelerations,
        vector<CollisionObject*>* collision_objects);

    void reset();

    void build_spatial_map();
    void self_collide(Particle& pm, double simulation_steps);
    float hash_position(Vector3D pos);
    
    // Cloth components
    vector<Particle> particles;

    // Particle properties
    vector<ParticleProperties> particleProperties;
    vector<Vector3D> particleColors;
    
    // marching cubes
    ScalarLoc* cells;
    int totalVertexCount;
    int sideCellCount;
    double cellSize;
    double borderDist;
    double physicsBorder;
    vector<Vector3D> planeLocs;
    vector<Vector3D> planeNorms;
    
    void initMarchingCubes(int numCells, double cellSize);
    MeshTriangle* getMarchingCubeMesh(int& numTriangles);

    // Spatial hashing
    unordered_map<float, vector<Particle*>*> map;
};

#endif /* CLOTH_H */
