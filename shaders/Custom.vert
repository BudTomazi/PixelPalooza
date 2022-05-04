#version 330

// Uniform variables are constant throughout the entire shader
// execution. They are also read-only to enable parallelization.
uniform mat4 u_model;
uniform mat4 u_view_projection;

// In a vertex shader, the "in" variables are read-only per-vertex 
// properties. An example of this was shown in the rasterizer project, 
// where each vertex had an associated "color" or "uv" value which we 
// would later interpolate using barycentric coordinates.
in vec4 in_position;
in vec4 in_normal;
in vec4 in_tangent;
in vec4 in_color;
in float in_shaderType;

// In a vertex shader, the "out" variables are per-vertex properties
// that are read/write. These properties allow us to communicate
// information from the vertex shader to the fragment shader.
// That is, in the linked fragment shader, these values become the 
// "in" variables.
out vec4 v_position;
out vec4 v_normal;
out vec4 v_tangent;
out vec4 v_color;
out float v_shaderType;

// Every shader features a "main" function.
// This is typically where we write to the "out" variables that the
// fragment shaders get to use. It is also where "gl_Position" is set,
// which is the final screen-space location of this vertex which the
// GPU's triangle rasterizer takes in.
void main() {
  // Here, we just apply the model's transformation to the various
  // per-vertex properties. That way, when the fragment shader reads
  // them, we already have the position in world-space.
  //vec4 new_pos = round(5.0 * in_position) / 5.0;
  vec4 new_pos = in_position;
  //new_pos.x = new_pos.x + 1;
  v_position = u_model * new_pos;
  //v_position = floor(5.0 * v_position) / 5.0;
  v_normal = normalize(u_model * in_normal);
  v_tangent = normalize(u_model * in_tangent);
    v_color = in_color;
    v_shaderType = in_shaderType;
  
  // The final screen-space location of this vertex which the
  // GPU's triangle rasterizer takes in.
  gl_Position = u_view_projection * u_model * new_pos;
}
