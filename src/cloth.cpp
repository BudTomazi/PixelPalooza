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

Cloth::Cloth(double width, double height, int num_width_points,
    int num_height_points, float thickness) {
    this->width = width;
    this->height = height;
    this->num_width_points = num_width_points;
    this->num_height_points = num_height_points;
    this->thickness = thickness;

    buildGrid();
    buildClothMesh();
}

Cloth::~Cloth() {
    point_masses.clear();
    springs.clear();

    particleProperties.clear();
    particleColors.clear();

    if (clothMesh) {
        delete clothMesh;
    }
}

void Cloth::buildGrid() {
    initializeProperties();
    Vector3D right = Vector3D(1, 0, 0);
    Vector3D up = Vector3D(0, 0, 1);

    for (int y = 0; y < num_height_points; y++) {
        for (int x = 0; x < num_width_points; x++) {
            Vector3D pos = Vector3D(0, 0, 0);
            int particleType = rand() % 2;
            pos += Vector3D(0.5 * (rand() % 100 / 100.0), 0.5 * (rand() % 100 / 100.0), 1.0 * (rand() % 100 / 100.0));

            double radius = particleProperties[particleType].radius;

            PointMass p = PointMass(pos, particleType);

            if (rand() % 10 == 0) {
                p.pinned = true;
            }
            else {
                p.pinned = false;
            }
            
            point_masses.push_back(p);
            point_masses.emplace_back(PointMass(pos + right * radius, particleType));
            point_masses.emplace_back(PointMass(pos + up * radius, particleType));
        }
    }
}

void Cloth::initializeProperties() {

    // example particles
    particleProperties.emplace_back(ParticleProperties(0.1, 0.5, Vector3D(1, 0, 0)));
    particleProperties.emplace_back(ParticleProperties(0.1, 0.5, Vector3D(0, 0.5, 1)));

    int numParticleTypes = particleProperties.size();

    //setting up forces and force strengths
    //particle type 1 laws
    particleProperties[0].strengths.resize(3);
    particleProperties[0].force_laws.push_back(&r2_law);
    particleProperties[0].strengths[0].insert(particleProperties[0].strengths[0].end(), { 0.2,0.0 });

    particleProperties[0].force_laws.push_back(&r4_law);
    particleProperties[0].strengths[1].push_back(-0.02);
    particleProperties[0].strengths[1].push_back(-0.2);

    particleProperties[0].force_laws.push_back(&cross_law);
    particleProperties[0].strengths[2].push_back(0.0);
    particleProperties[0].strengths[2].push_back(0.2);

    //particle type 2 laws
    particleProperties[1].strengths.resize(3);
    particleProperties[1].force_laws.emplace_back(&r2_law);
    particleProperties[1].strengths[0].emplace_back(0.1);
    particleProperties[1].strengths[0].emplace_back(0.2);

    particleProperties[1].force_laws.emplace_back(&r4_law);
    particleProperties[1].strengths[1].emplace_back(-0.2);
    particleProperties[1].strengths[1].emplace_back(-0.02);

    particleProperties[1].force_laws.emplace_back(&cross_law);
    particleProperties[1].strengths[2].emplace_back(0.2);
    particleProperties[1].strengths[2].emplace_back(0.0);

    for (int i = 0; i < numParticleTypes; i++) {
        particleColors.emplace_back(particleProperties[i].color);
    }

}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters* cp,
    vector<Vector3D> external_accelerations,
    vector<CollisionObject*>* collision_objects) {
    
    //Particle force calculations-----------------------------
    double delta_t = 1.0f / frames_per_sec / simulation_steps;
    double distSqr;
    Vector3D dir;
    ParticleProperties* curParticleProperties;
    ParticleProperties* otherParticleProperties;
    PointMass* curMass;
    PointMass* otherMass;
    Vector3D forces;
    double collDist;
    Vector3D centroid;

    //loop over all pairs of particles and compute forces
    for (int i = 0; i < point_masses.size(); i += 3) {
        curMass = &point_masses[i];
        curParticleProperties = &particleProperties[curMass->particle_type];
        int cType = curMass->particle_type;
        for (int j = 0; j < point_masses.size(); j += 3) {
            if (i == j) continue;
            forces = Vector3D(0);
            otherMass = &point_masses[j];
            int oType = otherMass->particle_type;
            otherParticleProperties = &particleProperties[otherMass->particle_type];

            dir = otherMass->position - curMass->position;
            distSqr = dir.norm2();
            dir = dir.unit();

            if (distSqr <= 0.000005) continue; //cancel computation if distance too small to avoid explosions

            for (int k = 0; k < otherParticleProperties->force_laws.size(); k++) {
                Vector3D force = -otherParticleProperties->force_laws[k](curMass, otherMass);
                float factor = otherParticleProperties->strengths[k][cType];
                forces += factor * force;
            }
        }

        if (!curMass->pinned) {
            curMass->forces = forces;
            //curMass->forces += 0.3* cross(curMass->position, Vector3D(1.0, 0.0, 0.0));
            curMass->forces += Vector3D(0.0, -0.7, 0.0);
        }
        centroid += curMass->position;
    }

    centroid /= (point_masses.size() / 3);

    for (auto p = point_masses.begin(); p != point_masses.end(); p++) {
        for (auto prim = collision_objects->begin(); prim != collision_objects->end(); prim++) {
            (*prim)->collide(*p);
        }
    }

    Vector3D curPos;
    Vector3D newPos;
    double delta_t_sqr = delta_t * delta_t;
    double dampingFactor = 1 - (cp->damping / 100.0);
    double delta = 0.1;
    Vector3D normal;
    PointMass* curParticle;

    //loop through and update particle positions via verlet integration
    for (int i = 0; i < point_masses.size(); i += 3) {
        curParticle = &point_masses[i];
        curParticleProperties = &particleProperties[curParticle->particle_type];

        curPos = curParticle->position;
        newPos = curPos + dampingFactor * (curPos - curParticle->last_position) + ((curParticle->forces / curParticleProperties->mass) * delta_t_sqr);

        // we cannot move farther than "delta" per timestep
        if ((curPos - newPos).norm2() > delta * delta) {
            newPos = curPos + (newPos - curPos).unit() * delta;
        }

        curParticle->position = newPos;
        curParticle->last_position = curPos;

        normal = (centroid - newPos).unit();
        Vector3D orthogonal = Vector3D(1, 1, -(normal.x + normal.y) / normal.z).unit();
        point_masses[i + 1].position = newPos + orthogonal * curParticleProperties->radius;
        point_masses[i + 2].position = newPos + cross(orthogonal, normal) * curParticleProperties->radius;
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

MeshTriangle* Cloth::getMarchingCubeMesh(int& numTriangles) {
    int size = 80; // number of cells along one side
    double width = 10; // how big the box is
    
    double c = width / size;
    
    int n = size + 1; // how many points? TODO: do we need to do this
    int yz = n * n;        //for little extra speed
    
    double min = -width / 2.0; // centers the entire grid
    
    int numCells = n * n * n;
    ScalarLoc* cells = new ScalarLoc[numCells];
    Vector3D bottomLeft = Vector3D(min);
    
    for (int i = 0; i < numCells; i++) {
        cells[i] = ScalarLoc(c * Vector3D(i / yz, (i / n) % n, i % n) + bottomLeft, 0);
    }
    
    for (int i = 0; i < point_masses.size(); i += 3) {
        Vector3D pos = point_masses[i].position;

        int index = (int)(floor((pos.x - min) / c) * yz + floor((pos.y - min) / c) * n + floor((pos.z - min) / c));
        if (index < 0 || index >= numCells) continue;
        cells[index].value += 0.3;
    }

    /*for (int i = 0; i < point_masses.size(); i += 3) {
        for (int j = 0; j < point_masses.size(); j += 3) {
            if (i == j) continue;
            Vector3D pos = 0.5 * (point_masses[i].position + point_masses[j].position);
            int index = (int)(floor((pos.x - min) / c) * yz + floor((pos.y - min) / c) * n + floor((pos.z - min) / c));
            if (index < 0 || index >= numCells) continue;

            cells[index].value += 0.01;
        }
    }*/
    
    MeshTriangle* triangles = MarchingCubesCross(size, size, size, 0.1f, cells, numTriangles);
    delete[] cells;
    
    return triangles;
}

void Cloth::buildClothMesh() {
    if (point_masses.size() == 0) return;

    ClothMesh* clothMesh = new ClothMesh();
    vector<Triangle*> triangles;

    Vector3D zero = Vector3D(0);
    for (int i = 0; i < point_masses.size(); i += 3) {
        triangles.push_back(new Triangle(&point_masses[i], &point_masses[i + 1], &point_masses[i + 2], zero, zero, zero));
    }

    // Create vector of triangles
    /*for (int y = 0; y < num_height_points - 1; y++) {
      for (int x = 0; x < num_width_points - 1; x++) {
        PointMass *pm = &point_masses[y * num_width_points + x];



        // Get neighboring point masses:
        *                      *
         * pm_A -------- pm_B   *
         *             /        *
         *  |         /   |     *
         *  |        /    |     *
         *  |       /     |     *
         *  |      /      |     *
         *  |     /       |     *
         *  |    /        |     *
         *      /               *
         * pm_C -------- pm_D   *
         *                      *
         *

        */

    // For each triangle in row-order, create 3 edges and 3 internal halfedges
    for (int i = 0; i < triangles.size(); i++) {
        Triangle* t = triangles[i];

        // Allocate new halfedges on heap
        Halfedge* h1 = new Halfedge();
        Halfedge* h2 = new Halfedge();
        Halfedge* h3 = new Halfedge();

        // Allocate new edges on heap
        Edge* e1 = new Edge();
        Edge* e2 = new Edge();
        Edge* e3 = new Edge();

        // Assign a halfedge pointer to the triangle
        t->halfedge = h1;

        // Assign halfedge pointers to point masses
        t->pm1->halfedge = h1;
        t->pm2->halfedge = h2;
        t->pm3->halfedge = h3;

        // Update all halfedge pointers
        h1->edge = e1;
        h1->next = h2;
        h1->pm = t->pm1;
        h1->triangle = t;

        h2->edge = e2;
        h2->next = h3;
        h2->pm = t->pm2;
        h2->triangle = t;

        h3->edge = e3;
        h3->next = h1;
        h3->pm = t->pm3;
        h3->triangle = t;
    }

    clothMesh->triangles = triangles;
    this->clothMesh = clothMesh;
}
