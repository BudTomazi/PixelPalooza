#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  
  // (Placeholder code. You will want to replace it.)
  out_color = (vec4(1, 1, 1, 0) + v_normal) / 2;
  out_color.a = 1;

  vec3 l = u_light_pos - vec3(v_position);
  vec3 v = u_cam_pos - vec3(v_position);
  vec3 h = (v + l) / length(v + l);

  vec3 ambient = vec3(0.3);

  vec3 diffuse = u_light_intensity / (length(l) * length(l)) 
  * max(0.0, dot(normalize(vec3(v_normal)), normalize(l)));

  vec3 specular = u_light_intensity / (length(l) * length(l)) 
  * pow(max(0.0, dot(normalize(vec3(v_normal)), h)), 40);

  //out_color = vec4(ambient + diffuse + 0.8 * specular, 1.0);
  out_color = vec4(diffuse, 1.0);
}

