#version 330

uniform mat4 Projection;

out vec3 fragPos;

vec4 positions[4] = vec4[](
	vec4(-1, -1, -1, 1),
	vec4( 1, -1, -1, 1), 
	vec4( 1,  1, -1, 1), 
	vec4(-1,  1, -1, 1) 
);

void main(){
	vec4 Position = positions[gl_VertexID];

	fragPos = vec3(inverse(Projection)*Position);
	gl_Position = Position;
} 
