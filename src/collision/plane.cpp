#include "iostream"
#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../clothSimulator.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide(PointMass &pm) {
  // TODO (Part 3): Handle collisions with planes.
    double last_side, curr_side;
    last_side = dot(pm.last_position - point, normal);
    curr_side = dot(pm.position - point, normal);
    //abs(curr_side - last_side) > abs(curr_side)

    bool a = last_side > 0;
    bool b = curr_side > 0;


    if (a != b) { // if the signs are different, it means the positions are on opposite sides of the plane
        //if errors, check < vs <= edge case
        Vector3D l0 = pm.last_position;
        Vector3D l = (pm.position - pm.last_position).unit();

        double d = dot((point - l0), normal) / dot(l, normal);

        Vector3D tangent_point = l0 + l * d;

        Vector3D correction = tangent_point - pm.last_position;
        correction = correction - correction.unit() * SURFACE_OFFSET;
        pm.position = pm.last_position + correction * (1 - friction);

    }

}

void Plane::render(GLShader &shader) {
  nanogui::Color color(0.7f, 0.7f, 0.7f, 1.0f);

  Vector3f sPoint(point.x, point.y, point.z);
  Vector3f sNormal(normal.x, normal.y, normal.z);
  Vector3f sParallel(normal.y - normal.z, normal.z - normal.x,
                     normal.x - normal.y);
  sParallel.normalize();
  Vector3f sCross = sNormal.cross(sParallel);

  MatrixXf positions(3, 4);
  MatrixXf normals(3, 4);

  positions.col(0) << sPoint + 2 * (sCross + sParallel);
  positions.col(1) << sPoint + 2 * (sCross - sParallel);
  positions.col(2) << sPoint + 2 * (-sCross + sParallel);
  positions.col(3) << sPoint + 2 * (-sCross - sParallel);

  normals.col(0) << sNormal;
  normals.col(1) << sNormal;
  normals.col(2) << sNormal;
  normals.col(3) << sNormal;

  if (shader.uniform("u_color", false) != -1) {
    shader.setUniform("u_color", color);
  }
  shader.uploadAttrib("in_position", positions);
  if (shader.attrib("in_normal", false) != -1) {
    shader.uploadAttrib("in_normal", normals);
  }

  shader.drawArray(GL_TRIANGLE_STRIP, 0, 4);
}
