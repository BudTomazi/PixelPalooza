#version 330

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec4 v_color;

out vec4 out_color;

void main() {
    out_color = v_normal / 2;
//    out_color = v_color * (2.0 / 3.0) + (v_normal) / 3;
  out_color.a = 1;
}
