#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
  // texture(u_texture_1, v_uv)
  return texture(u_texture_2, uv).r;
}

void main() {
  // YOUR CODE HERE
  /*float dU = (h(v_uv + vec2(1/u_texture_2_size.y, 0.0)) - h(v_uv)) * u_height_scaling * u_normal_scaling;
  float dV = (h(v_uv + vec2(0.0, 1/u_texture_2_size.x)) - h(v_uv)) * u_height_scaling * u_normal_scaling;
  vec3 no = vec3(-dU, -dV, 1.0);
  mat3 TBN = mat3(vec3(v_tangent), cross(vec3(v_normal), vec3(v_tangent)), vec3(v_normal));

  vec3 nd = TBN * no;


  // (Placeholder code. You will want to replace it.)
  //out_color = vec4(nd, 1.0);
  //out_color.a = 1;

  vec3 l = u_light_pos - vec3(v_position);
  vec3 v = u_cam_pos - vec3(v_position);
  vec3 h = (v + l) / length(v + l);

  vec3 ambient = vec3(0.0);

  vec3 diffuse = u_light_intensity / (length(l) * length(l)) 
  * max(0.0, dot(normalize(vec3(nd)), normalize(l)));

  vec3 specular = u_light_intensity / (length(l) * length(l)) 
  * pow(max(0.0, dot(normalize(vec3(nd)), h)), 40);

  out_color = vec4(ambient + diffuse + 0.8 * specular, 1.0);*/
}

