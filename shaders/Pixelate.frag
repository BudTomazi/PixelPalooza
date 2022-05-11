#version 330

out vec4 FragColor;
in vec2 texCoords;

uniform sampler2D screenTexture;
const float pix = 175;

void main() {
    FragColor = texture(screenTexture,
                        vec2(floor(texCoords.x * pix) / pix, floor(texCoords.y * pix) / pix));
}
