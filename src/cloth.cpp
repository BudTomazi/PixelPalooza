#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

//#include <nanoflann.hpp>

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
    interactionProperties.clear();

    if (clothMesh) {
        delete clothMesh;
    }
}

void Cloth::buildGrid() {
    initializeProperties();

    // TODO (Part 1): Build a grid of masses and springs.
  //    std::vector<std::vector<double> > dann;
  //    size_t num = 10;
  //    dann.resize(num);
  //    for (int i = 0; i < num; i++) {
  //        dann[i].resize(3);
  //
  //        for (size_t d = 0; d < 3; d++)
  //            dann[i][d] = i;
  //    }
  //
  //    KDTreeVectorOfVectorsAdaptor< std::vector<std::vector<double> >, double > mat_index(3 /*dim*/, dann, 10 /* max leaf */ );

    Vector3D right = Vector3D(1, 0, 0);
    Vector3D up = Vector3D(0, 0, 1);

    for (int y = 0; y < num_height_points; y++) {
        for (int x = 0; x < num_width_points; x++) {
            Vector3D pos = Vector3D(0, 1, 0);
            int particleType = rand() % 2;

            pos += Vector3D(2 * (particleType % 2) + rand() % 100 / 100.0 - 0.5, rand() % 100 / 100.0 - 0.5, rand() % 100 / 100.0 - 0.5);

            double radius = particleProperties[particleType].radius;

            point_masses.emplace_back(PointMass(pos, particleType));
            point_masses.emplace_back(PointMass(pos + right * radius, particleType));
            point_masses.emplace_back(PointMass(pos + up * radius, particleType));
        }
    }

    /*
     Vector3D right = Vector3D(0.1, 0, 0);
     Vector3D up = Vector3D(0, 0, 0.1);

     for (int y = 0; y < num_height_points; y++) {
         for (int x = 0; x < num_width_points; x++) {
             Vector3D pos = Vector3D(x * width / num_width_points, 2 * (rand() % 1000) / 1000.0 + 1, y * height / num_height_points);

             point_masses.emplace_back(PointMass(pos, false));
             point_masses.emplace_back(PointMass(pos + right, false));
             point_masses.emplace_back(PointMass(pos + up, false));
         }
     }
     */
}

void Cloth::initializeProperties() {
    // example particles
    particleProperties.emplace_back(ParticleProperties(1, 0.5, Vector3D(1, 0, 0)));
    particleProperties.emplace_back(ParticleProperties(1, 0.5, Vector3D(0, 0.5, 1)));
    particleProperties.emplace_back(ParticleProperties(10, 0.03, Vector3D(0.2, 1, 0)));

    int numParticleTypes = particleProperties.size();
    interactionProperties.resize(numParticleTypes);
    for (int i = 0; i < numParticleTypes; i++) {
        particleColors.emplace_back(particleProperties[i].color);
        interactionProperties[i].resize(numParticleTypes);
    }

    interactionProperties[0][0] = InteractionProperties(0.2, 0.01);
    interactionProperties[0][1] = InteractionProperties(0.0, 0.1);
    interactionProperties[1][0] = InteractionProperties(0.0, 0.1);
    interactionProperties[1][1] = InteractionProperties(0.2, 0.01);
    //    interactionProperties[0][2] = InteractionProperties(0, 0.06);

         // TODO can use less memory by not storing these?

    //    interactionProperties[1][2] = InteractionProperties(0.05, 0.005);
    //
    //    interactionProperties[2][0] = InteractionProperties(0, 0.06);
    //    interactionProperties[2][1] = InteractionProperties(0.05, 0.005);
    //    interactionProperties[2][2] = InteractionProperties(0.1, 0.005);

}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters* cp,
    vector<Vector3D> external_accelerations,
    vector<CollisionObject*>* collision_objects) {
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    double distSqr;
    Vector3D dir;
    InteractionProperties* curInteractionProperties;
    ParticleProperties* curParticleProperties;
    ParticleProperties* otherParticleProperties;
    PointMass* curMass;
    PointMass* otherMass;
    Vector3D forces;
    double collDist;

    Vector3D centroid;
    for (int i = 0; i < point_masses.size(); i += 3) {
        curMass = &point_masses[i];
        forces = Vector3D(0);
        curParticleProperties = &particleProperties[curMass->particle_type];
        
        for (int j = 0; j < point_masses.size(); j += 3) {
            if (i == j) continue;

            otherMass = &point_masses[j];
            otherParticleProperties = &particleProperties[otherMass->particle_type];

            dir = otherMass->position - curMass->position;
            distSqr = dir.norm2();
            dir = dir.unit();

            if (distSqr <= 0.0005) continue;

            curInteractionProperties = &interactionProperties[curMass->particle_type][otherMass->particle_type];
            if (!curInteractionProperties->isActive) continue;

            forces += dir * (curInteractionProperties->attract / distSqr);
            forces -= dir * (curInteractionProperties->repel / (distSqr * distSqr));

            /*collDist = (otherParticleProperties->radius) + (curParticleProperties->radius);
            if (distSqr < (collDist * collDist) && curMass->particle_type == otherMass->particle_type) {
                curMass->particle_type = (curMass->particle_type + 1) % 2;
            }*/
        }

        curMass->forces = forces;

        centroid += curMass->position;
    }

    centroid /= (point_masses.size() / 3);

    Vector3D curPos;
    Vector3D newPos;
    double delta_t_sqr = delta_t * delta_t;
    double dampingFactor = 1 - (cp->damping / 100.0);

    double delta = 0.1;

    Vector3D normal;

    PointMass* curParticle;

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

    // TODO (Part 2): Compute total force acting on each point mass.


    // TODO (Part 2): Use Verlet integration to compute new point mass positions


    // TODO (Part 4): Handle self-collisions.


    // TODO (Part 3): Handle collisions with other primitives.


    // TODO (Part 2): Constrain the changes to be such that the spring does not change
    // in length more than 10% per timestep [Provot 1995].

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

        float u_min = x;
        u_min /= num_width_points - 1;
        float u_max = x + 1;
        u_max /= num_width_points - 1;
        float v_min = y;
        v_min /= num_height_points - 1;
        float v_max = y + 1;
        v_max /= num_height_points - 1;

        PointMass *pm_A = pm                       ;
        PointMass *pm_B = pm                    + 1;
        PointMass *pm_C = pm + num_width_points    ;
        PointMass *pm_D = pm + num_width_points + 1;

        Vector3D uv_A = Vector3D(u_min, v_min, 0);
        Vector3D uv_B = Vector3D(u_max, v_min, 0);
        Vector3D uv_C = Vector3D(u_min, v_max, 0);
        Vector3D uv_D = Vector3D(u_max, v_max, 0);


        // Both triangles defined by vertices in counter-clockwise orientation
        triangles.push_back(new Triangle(pm_A, pm_C, pm_B,
                                         uv_A, uv_C, uv_B));
        triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
                                         uv_B, uv_C, uv_D));
      }
    }*/

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

    // Go back through the cloth mesh and link triangles together using halfedge
    // twin pointers

    // Convenient variables for math
    /*int num_height_tris = (num_height_points - 1) * 2;
    int num_width_tris = (num_width_points - 1) * 2;

    bool topLeft = true;
    for (int i = 0; i < triangles.size(); i++) {
      Triangle *t = triangles[i];

      if (topLeft) {
        // Get left triangle, if it exists
        if (i % num_width_tris != 0) { // Not a left-most triangle
          Triangle *temp = triangles[i - 1];
          t->pm1->halfedge->twin = temp->pm3->halfedge;
        } else {
          t->pm1->halfedge->twin = nullptr;
        }

        // Get triangle above, if it exists
        if (i >= num_width_tris) { // Not a top-most triangle
          Triangle *temp = triangles[i - num_width_tris + 1];
          t->pm3->halfedge->twin = temp->pm2->halfedge;
        } else {
          t->pm3->halfedge->twin = nullptr;
        }

        // Get triangle to bottom right; guaranteed to exist
        Triangle *temp = triangles[i + 1];
        t->pm2->halfedge->twin = temp->pm1->halfedge;
      } else {
        // Get right triangle, if it exists
        if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
          Triangle *temp = triangles[i + 1];
          t->pm3->halfedge->twin = temp->pm1->halfedge;
        } else {
          t->pm3->halfedge->twin = nullptr;
        }

        // Get triangle below, if it exists
        if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
          Triangle *temp = triangles[i + num_width_tris - 1];
          t->pm2->halfedge->twin = temp->pm3->halfedge;
        } else {
          t->pm2->halfedge->twin = nullptr;
        }

        // Get triangle to top left; guaranteed to exist
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm2->halfedge;
      }

      topLeft = !topLeft;
    }*/

    clothMesh->triangles = triangles;
    this->clothMesh = clothMesh;
}
