#version 330


uniform vec3 u_cam_pos;

uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  //out_color = (vec4(1, 1, 1, 0) + v_normal) / 2;
  //out_color.a = 1;

  vec3 in_ray = u_cam_pos - vec3(v_position);
  vec3 sideways = (in_ray - dot(in_ray, vec3(v_normal)));
  vec3 out_ray = in_ray - 2 * sideways;

  out_color = texture(u_texture_cubemap, out_ray);
  out_color.a = 1.0;

}
