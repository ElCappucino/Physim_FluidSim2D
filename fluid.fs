#version 330 core
out vec4 FragColor;
in vec2 TexCoords;

uniform sampler2D smokeTex;
uniform bool showSmoke;

void main()
{
    // Sample the smoke density (stored in Red channel)
    float s = texture(smokeTex, vec2(TexCoords.y, TexCoords.x)).r;
    
    vec3 color = vec3(1.0); // Default White background

    if (showSmoke) {
        color = vec3(s); 
    }

    FragColor = vec4(color, 1.0);
}