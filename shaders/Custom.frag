#version 330

// (Every uniform is available here.)

uniform mat4 u_view_projection;
uniform mat4 u_model;

uniform float u_normal_scaling;
uniform float u_height_scaling;

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

// Feel free to add your own textures. If you need more than 4,
// you will need to modify the skeleton.
uniform sampler2D u_texture_1;
uniform sampler2D u_texture_2;
uniform sampler2D u_texture_3;
uniform sampler2D u_texture_4;

// Environment map! Take a look at GLSL documentation to see how to
// sample from this.
uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec4 v_color;
in float v_shaderType;

out vec4 out_color;

void main() {
  // Your awesome shader here!
    if (v_shaderType < 0.01) {
//        out_color = v_color * (2.0 / 3.0) + (v_normal) / 3;
        float intensity;
        intensity = dot(vec4(u_light_pos.x, u_light_pos.y, u_light_pos.z, 1), normalize(v_normal));
        if (intensity > 0.9) {
            out_color = v_color * 0.83;
        } else if (intensity > 0.5) {
            out_color = v_color * 0.71;
        } else if (intensity > 0.1) {
            out_color = v_color * 0.56;
        } else {
            out_color = v_color * 0.39;
        }
    } else if (v_shaderType < 1.01) {
        // cursed mirror
        vec3 in_ray = u_cam_pos - vec3(v_position);
        vec3 sideways = (in_ray - dot(in_ray, vec3(v_normal)));
        vec3 out_ray = in_ray - 2 * sideways;

        out_color = v_color / 3.0 + texture(u_texture_cubemap, out_ray) * (2.0 / 3.0);
    } else {
        out_color = vec4(0, 0, 0, 1);
    }
//    out_color = v_normal;
  out_color.a = 1;
}
