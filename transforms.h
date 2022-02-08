#ifndef TRANSFORMS_H
#define TRANSFORMS_H
#define _USE_MATH_DEFINES
#define M_PI           3.14159265358979323846

#include <cmath>
#include "matrix.h"
#include "vec.h"

inline mat4 loadIdentity(){
    return Id<4>();
}

inline mat4 translate(vec3 u){
	return {
		1, 0, 0, u(0),
		0, 1, 0, u(1),
		0, 0, 1, u(2),
		0, 0, 0,    1
	};
}

inline mat4 translate(float x, float y, float z){
    return translate(vec3{x, y, z});
}

inline mat4 scale(float x, float y, float z){
	return {
		x, 0, 0, 0,
		0, y, 0, 0,
		0, 0, z, 0,
		0, 0, 0, 1
	};
}

inline mat4 rotate_x(float angle){
	float c = cos(angle);
	float s = sin(angle);

	return {
		1, 0,  0, 0,
		0, c, -s, 0,
		0, s,  c, 0,
		0, 0,  0, 1
	};
}

inline mat4 rotate_y(float angle){
	float c = cos(angle);
	float s = sin(angle);

	return {
		 c, 0, s, 0,
		 0, 1, 0, 0,
		-s, 0, c, 0,
		 0, 0, 0, 1
	};
}

inline mat4 rotate_z(float angle){
	float c = cos(angle);
	float s = sin(angle);

	return {
		c, -s, 0, 0,
		s,  c, 0, 0,
		0,  0, 1, 0,
		0,  0, 0, 1
	};
}

inline vec3 rodrigues(vec3 eixo,vec3 uni,float theta){
  //vec3 uni={1,1,1};
  uni=normalize(uni);
  return (dot( eixo,uni))*(1-cos(theta))*uni+
                  cos(theta)*eixo+
                  sin(theta)*cross(uni,eixo);
}

inline mat4 rotate(vec3 n, float theta){
  vec3 eixo1={1,0,0};
  vec3 eixo2={0,1,0};
  vec3 eixo3={0,0,1};

  vec3 rodrigues1= rodrigues(eixo1,n,theta);
  vec3 rodrigues2= rodrigues(eixo2,n,theta);
  vec3 rodrigues3= rodrigues(eixo3,n,theta);

	return {
    rodrigues1[0],rodrigues2[0],rodrigues3[0],0,
    rodrigues1[1],rodrigues2[1],rodrigues3[1],0,
    rodrigues1[2],rodrigues2[2],rodrigues3[2],0,
    0,0,0,1
  };
}

inline mat4 lookAt(vec3 eye, vec3 center, vec3 up){
	vec3 upNormalizado=normalize(up);
  vec3 f=normalize((eye-center));
  f=-f;
  vec3 s=normalize(cross(f,upNormalizado));
  vec3 u=cross(s,f);

  mat4 t1={
    s[0],s[1],s[2],0,
    u[0],u[1],u[2],0,
    -f[0],-f[1],-f[2],0,
    0,0,0,1
  };
  mat4 t2={
    1,0,0,-eye[0],
    0,1,0,-eye[1],
    0,0,1,-eye[2],
    0,0,0,1
  };

	return t1*t2;
}

inline mat4 orthogonal(float l, float r, float b, float t, float n, float f){
	return {
		2/(r-l),      0,     0,      (l+r)/(l-r),
			0,  2/(t-b),     0,      (b+t)/(b-t),
			0,        0, -2/(f-n),   (n+f)/(n-f),
			0,        0,      0,               1
	};
}

inline mat4 frustum(float l, float r, float b, float t, float n, float f){
	mat4 Persp = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0,-(n+f)/n, -f,
		0, 0, -1/n, 0
	};
	mat4 Volume = {
		2/(r-l),      0,     0,      (l+r)/(l-r),
			0,  2/(t-b),     0,      (b+t)/(b-t),
			0,        0, 2/(f-n),    (n+f)/(n-f),
			0,        0,      0,               1
	};
	return Volume*Persp;
}

double GrausParaRadianos(double degrees){
    return degrees * M_PI / 180;
}

inline mat4 perspective(float fovy, float aspect, float Near, float Far){
    float t,b,r,l;
    t=Near*tan(GrausParaRadianos(fovy)/2);
    b=-1*t;
    r=t*aspect;
    l=-r;
    return frustum(l,r,b,t,Near,Far);
}

#endif
