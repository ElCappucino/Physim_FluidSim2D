#version 330 core
out vec4 FragColor;
in vec2 TexCoord;

// You can eventually pass these as uniforms from your main script
uniform vec3 obstacleColor = vec3(1.0, 1.0, 1.0); // White
uniform vec3 outlineColor = vec3(0.0, 0.0, 0.0);  // Black
uniform float outlineThickness = 0.01;           // How thick the ring is

void main() {

    vec2 center = vec2(0.5, 0.5);
    float dist = distance(TexCoord, center);

    if (dist > 0.5) {
        discard;
    }

    float edgeSmoothing = 0.005; 
    float alpha = 1.0 - smoothstep(0.5 - edgeSmoothing, 0.5, dist);

    float outlineFactor = smoothstep(0.5 - outlineThickness - edgeSmoothing, 0.5 - outlineThickness, dist);
    
    vec3 finalColor = mix(obstacleColor, outlineColor, outlineFactor);

    FragColor = vec4(finalColor, alpha);
}