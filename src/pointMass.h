#ifndef POINTMASS_H
#define POINTMASS_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"

using namespace CGL;

// Forward declarations
class Halfedge;

struct PointMass {
    PointMass(Vector3D position, int particle_type)
        : start_position(position), position(position),
        last_position(position), particle_type(particle_type) {}

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

    int particle_type;

    // mesh reference
    Halfedge* halfedge;
};

struct ParticleProperties {
    ParticleProperties(double mass, double radius, Vector3D color) :
        mass(mass), radius(radius), color(color) {}

    double mass;
    double radius;
    Vector3D color;
    // external forces?
    // color
};

struct InteractionProperties {
    InteractionProperties() : attract(0), repel(0), isActive(false) {}

    InteractionProperties(double attract, double repel) :
        attract(attract), repel(repel), isActive(true) {}
    double attract;
    double repel;
    bool isActive;
    // collision function?

};



#endif /* POINTMASS_H */
