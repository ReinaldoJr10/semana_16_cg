#ifndef RAYTRACING_H
#define RAYTRACING_H

#include "Primitives.h"
#include "vec.h"
#include <algorithm>

inline vec3 reflect(vec3 I, vec3 N){
    return I - 2*dot(N, I)*N;
}

// Menor solucao positiva da equacao ax^2 + bx + c = 0
// Retorna um valor negativo se nao houver solucao positiva
inline float smallerT(float a, float b, float c){
	float delta = b*b - 4*a*c;
	const float EPS = 1e-50;
	if(delta < -EPS)
		return -1.0;

	if(delta < EPS){
		float t = -b/(2*a);
		return t;
	}

    float t1 = (-b-sqrt(delta))/(2*a);
    float t2 = (-b+sqrt(delta))/(2*a);

    if(t1 < 0)
        return t2;

    if(t2 < 0)
        return t1;

    return std::min(t1, t2);
}


struct Sphere{
    float radius;
    vec3 center;
};

struct SphereIntersection{
    float t;
    Sphere sphere;

    bool operator<(const SphereIntersection& other) const{
        return t < other.t;
    }
};

template<class Iterator>
struct TriangleIntersection{
    float t;
    float u;
    float v;
    Iterator it;

    bool operator<(const TriangleIntersection& other) const {
        return t < other.t;
    }
};

inline Triangle<vec3> get_triangle(const Triangle<vec3>& tri){
    return tri;
}

struct Ray{
    vec3 orig;
    vec3 dir;

    vec3 at(float t) const{
        return orig + t*dir;
    }
    
    SphereIntersection sphere_intersection(Sphere sphere) const{ 
        vec3 center = sphere.center;
        vec3 dif = orig - center;
        
        float a = dot(dir, dir);
        float b = 2*dot(dif, dir);
        float c = dot(dif, dif) - sphere.radius*sphere.radius;
        
        float t = smallerT(a, b, c);
        return {t < 1e-50? HUGE_VALF: t, sphere};
    } 

    SphereIntersection min_sphere_intersection(const std::vector<Sphere>& spheres) const{
        SphereIntersection min_intersection;
        min_intersection.t = HUGE_VALF;
        for(Sphere S: spheres)
            min_intersection = std::min(min_intersection, sphere_intersection(S));
        return min_intersection;
    }

    inline bool intersect(const Triangle<vec3>& P, float &t, float &u, float &v) const{ 
        // Algoritmo de Moller-Trumbore
        // Resolver o sistema: [E1 E2 -D][u v t]^T = C0
        // onde E1 = P1 - P0
        //      E2 = P2 - P0
        //      C0 =  C - P0

        vec3 E1 = P[1] - P[0]; 
        vec3 E2 = P[2] - P[0]; 
        vec3 DxE2 = cross(dir, E2); 
        float det = dot(E1, DxE2); 

        // ray and triangle are parallel if det is close to 0
        if (fabs(det) < 1e-15) 
            return false; 

        float invDet = 1/det; 
        vec3 C0 = orig - P[0]; 

        // u =  det(C0, E2, -D) / det(E1, E2, -D) 
        //   = (C0 . (D x E2) ) / (E1 . (D x E2) )
        u = dot(C0, DxE2)*invDet; 
        if (u < 0 || u > 1) 
            return false; 
    
        vec3 C0xE1 = cross(C0, E1); 

        // v =  det(E1, C0, -D) / det(E1, E2, -D) 
        //   = (D . (C0 x E1) ) / (E1 . (D x E2) )
        v = dot(dir, C0xE1)*invDet; 
        if (v < 0 || u + v > 1) 
            return false; 
    
        // t =  det(E1, E2, C0) / det(E1, E2, -D) 
        //   = (E2 . (C0xE1) ) / (E1 . (D x E2) )
        t = dot(E2, C0xE1)*invDet; 
    
        return t > 1e-3; 
    } 
    
    template<class Iterator>
    TriangleIntersection<Iterator> min_tri_intersection(Iterator b, Iterator e) const{
        TriangleIntersection<Iterator> res;
        res.t = HUGE_VALF;
        for(Iterator it = b; it != e; it++){
            float t, u, v;
            if(intersect(get_triangle(*it), t, u, v) && t < res.t)
                res = {t, u, v, it};
        }
        return res;
    }
};

Ray operator*(mat4 M, Ray ray){
    return {
        toVec3(M*toVec4(ray.orig, 1)),
        toVec3(M*toVec4(ray.dir, 0)),
    };
}

class Camera{
    // extrinsic parameters
    vec3 position = vec3{0, 0, 0};
    vec3 direction = vec3{0, 0, -1};
    vec3 up = vec3{0, 1, 0};
    vec3 right = vec3{1, 0, 0};

    // intrinsic parameters
    float near;
    float fovy;
    int w, h;

    public:
    Camera(float near, float fovy, int w, int h) : 
        near{near},fovy{fovy},w{w},h{h}{}

    void lookAt(vec3 pos, vec3 target, vec3 _up){
        position = pos;
        direction = normalize(target - position);
        right = normalize(cross(direction, _up));
        up = cross(right, direction);
    }

    Ray ray(float i, float j) const{
        float ymax = near*tan(fovy/2);
        float xmax = (ymax*w)/h;
        float x = 2*xmax*(i+0.5f)/w - xmax;
        float y = 2*ymax*(j+0.5f)/h - ymax;
        return {position, near*direction + x*right + y*up};
    }
};

vec3 most_distant_point(const std::vector<vec3>& points, vec3 ref){
    return *std::max_element(points.begin(), points.end(), 
        [&](const vec3& p, const vec3& q){ 
            return norm(p-ref) < norm(q-ref);
        }
    );
}

class DummyBoundingVolume{
    public:
    DummyBoundingVolume() = default;
    DummyBoundingVolume(const std::vector<vec3>&){}
    bool intersect(Ray) const{ return true; }
};

class BoundingSphere{
    Sphere bounding_sphere;

    public:
    BoundingSphere() = default;
    BoundingSphere(const std::vector<vec3>& position){
        if(position.size() < 1)
            return;

        vec3 py = most_distant_point(position, position[0]);
        vec3 pz = most_distant_point(position, py);

        bounding_sphere.center = 0.5*(py+pz);
        bounding_sphere.radius = 0.5*norm(py-pz);

        auto cposition = position;
        while(!cposition.empty()){
            auto it = std::find_if(cposition.begin(), cposition.end(), 
                [&](const vec3& p){ 
                    return norm(p-bounding_sphere.center) > bounding_sphere.radius; 
                });

            if(it == cposition.end())
                break;
            
            float A = norm(*it - bounding_sphere.center);
            float t = (A - bounding_sphere.radius)/(2*A);
            bounding_sphere.center = lerp(t, bounding_sphere.center, *it);
            bounding_sphere.radius = norm(*it - bounding_sphere.center);
            cposition.erase(it);
        }
    }

    bool intersect(Ray ray) const{
        SphereIntersection I = ray.sphere_intersection(bounding_sphere);
        return I.t < HUGE_VALF;
    }
};

class BoundingBox{
    bool init = false;
    vec3 pmin, pmax;
    public:
    BoundingBox() = default;
    
    BoundingBox(const std::vector<vec3>& verts){
        for(vec3 v: verts)
            add(v);
    }

    template<class Iterator>
    BoundingBox(Iterator b, Iterator e){
        for(auto it = b; it != e; it++)
            for(vec3 v: get_triangle(*it))
                add(v);
    }

    vec3 mid_point() const{
        return 0.5*(pmin + pmax);
    }

    void add(vec3 p){
        if(!init){
            pmin = pmax = p;
            init = true;
        }else{
            for(int i = 0; i < 3; i++){
                pmin(i) = std::min(pmin(i), p(i));
                pmax(i) = std::max(pmax(i), p(i));
            }
        }
    }

    bool intersect(Ray ray) const{
        for(vec3 v: {pmin, pmax}){
            for(int i=0; i < 3; i++){
                if(fabs(ray.dir[i]) < 1e-5)
                    continue;
                float t = (v[i] - ray.orig[i])/ray.dir[i];
                vec3 p = ray.at(t);
                int j = (i+1)%3;
                int k = (i+2)%3;
                
                if(
                    pmin[j] <= p[j] && p[j] <= pmax[j] && 
                    pmin[k] <= p[k] && p[k] <= pmax[k]
                )
                   return true;
            }
        }
        return false;
    }
};

template<class Tri>
class Octree{
    using Iterator = typename std::vector<Tri>::iterator;
    using Intersection = TriangleIntersection<Iterator>;
    Iterator b, e;
    int level;
    BoundingBox bounding_box;
    std::vector<Octree> children;

    public:
    Octree() = default;

    Octree(Iterator _b, Iterator _e, int level) : 
        b{_b}, e{_e}, level{level}, bounding_box{b, e}
    {
        if(level > 0)
            partition(0, bounding_box.mid_point(), _b, _e);
    }

    Intersection min_tri_intersection(Ray ray) const{
        Intersection res;
        res.t = HUGE_VALF;

        if(!bounding_box.intersect(ray))
           return res;
        
        if(children.empty())
            res = ray.min_tri_intersection(b, e);
        else{
           for(const Octree& child: children)
                res = std::min(res, child.min_tri_intersection(ray));
        }

        return res;
    }
    
    private:
    void partition(int i, vec3 c, Iterator beg, Iterator end){
        Iterator it = std::partition(beg, end, 
            [&](const Tri& tri){ 
                Triangle<vec3> T = get_triangle(tri);
                return T[0](i) < c(i) && 
                       T[1](i) < c(i) && 
                       T[2](i) < c(i) ;
            }
        );
        if(i < 2){ 
            partition(i+1, c, beg, it);
            partition(i+1, c, it, end);
        }else{
            add_child(beg, it);
            add_child(it, end);
        }
    }

    void add_child(Iterator beg, Iterator end){
        if(beg != end)
            children.emplace_back(beg, end, level-1);
    }
};


#endif