#ifndef POINTMASS_H
#define POINTMASS_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"
#include <vector>

using namespace CGL;

// Forward declarations
class Halfedge;

struct Particle {
    Particle(Vector3D position, int particle_type)
        : start_position(position), position(position),
        last_position(position), particle_type(particle_type), pinned(false) {}

    Vector3D normal();
    Vector3D velocity(double delta_t) {
        return (position - last_position) / delta_t;
    }

    // static values
    Vector3D start_position;

    // dynamic values
    Vector3D position;
    Vector3D last_position;
    Vector3D forces;
    bool pinned;

    int particle_type;

    // mesh reference
    Halfedge* halfedge;
};

typedef Vector3D (*forceLaw)(Particle* p1, Particle* p2);
struct ParticleProperties {
    ParticleProperties() {}
    ParticleProperties(double mass, double radius, Vector3D color) :
        mass(mass), radius(radius), color(color) {}

    double mass;
    double radius;
    double collRadius;
    Vector3D color;
    bool external_forces;
    bool primitive_collision;
    bool particle_collisions;
    bool pinned;
    std::vector<bool> localized;
    
    int shaderType;

    //vector<f> forces;
    //vector<float> strengths;
    std::vector<forceLaw> force_laws;
    std::vector<std::vector<float>> strengths;
    std::vector<int> collision_transformations;
    // external forces?
    // color
};

#endif /* POINTMASS_H */
