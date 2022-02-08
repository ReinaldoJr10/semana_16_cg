#include "Image.h"
#include "raytracing.h"
#include "transforms.h"
#include <random>

//#define USE_BOUNDING_SPHERE
#define USE_BOUNDING_BOX
#define USE_OCTREE
#include "RTMesh.h"

static double random_value(double a, double b){
    static std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(a, b);
    return distribution(generator);
}

// Direção aleatória ponderada pelo cosseno
vec3 cos_random_dir(vec3 n){
    double r1 = random_value(0, 2*M_PI); 
    double r2 = random_value(0, 1); 

    // base ortonormal: n, u, v
    vec3 u = normalize( cross((fabsf(n(0)) > .1 ? vec3{0, 1, 0} : vec3{1, 0, 0}), n) );
    vec3 v = cross(n, u);

    // direção aleatória dentro do hemisfério
    return normalize(sqrt(r2)*(cos(r1)*u + sin(r1)*v) + sqrt(1 - r2)*n);
}

vec3 sky_color(){
    return vec3{1, 1, 1}; // cor do céu
}

void show_progress(int current_value, int max_value){
    float percent = (float)current_value/(float)max_value*100;
    fprintf(stderr, "\r%.1f%%", percent);
}

enum IlluminationModel{
    EMISSION, DIFFUSE, SPECULAR
};

IlluminationModel get_illumination_model(MaterialInfo material){
    if(material.illum == 0)
        return EMISSION;
    
    if(material.illum == 1)
        return DIFFUSE;

    // Choose randomly between diffuse and specular
    return (random_value(0, 1) < 0.5)? DIFFUSE : SPECULAR;
}

struct Scene{
    Image image;
    Camera camera;
    std::vector<RTMesh> meshes;

    void render(int nsamples){
        #pragma omp parallel for schedule(dynamic, 1)  // OpenMP
        for(int y = 0; y < image.height; y++){
            show_progress(y, image.height-1);
            for(int x = 0; x < image.width; x++)
                image(x, y) = color_at(x,y, nsamples);
        }
    }

    Color color_at(int x, int y, int nsamples) const{
        vec3 col = {0, 0, 0};
        for(int i = 0; i < nsamples; i++){
            float rx = x + random_value(-0.1, 0.1);
            float ry = y + random_value(-0.1, 0.1);
            Ray ray = camera.ray(rx, ry);
            col = col + trace_path(ray, 0);
        }
        return toColor(1.0/nsamples*col);
    }

    vec3 trace_path(Ray ray, int depth) const{
        if(depth > 10)
            return vec3{0, 0, 0};

        MatTriIntersection I = min_intersection(ray, meshes);
        if(I.t == HUGE_VALF)
            return sky_color();
        
        vec3 n = normalize(I.normal); 
        sample_textures(I.material, I.texCoords);

        IlluminationModel model = get_illumination_model(I.material);

        if(model == EMISSION)
            return I.material.Kd;

        if(model == DIFFUSE){
            vec3 dir = cos_random_dir(n);
            return I.material.Kd*trace_path(Ray{I.position, dir}, depth+1);
        }

        if(model == SPECULAR){
            vec3 dir = reflect(ray.dir, n);
            return I.material.Ks*trace_path(Ray{I.position, dir}, depth+1);
        }

        return vec3{0, 0, 0};
    }
};

MaterialInfo get_material(int illum, vec3 Kd, vec3 Ks=vec3{0,0,0}){
    MaterialInfo mat = standard_material();
    mat.illum = illum;
    mat.Kd = Kd;
    mat.Ks = Ks;
    return mat;
}

std::vector<RTMesh> get_meshes(){
    std::vector<RTMesh> meshes;
    meshes.emplace_back("modelos/monkey.obj", translate(0, 5.6, -2)*scale(1.4, 1.4, 1.4)*rotate_x(-0.7));
    meshes.emplace_back("modelos/bunny.obj", translate(0, 5.2, 2), get_material(1, {0, 0, 1}));
    meshes.emplace_back("modelos/teapot.obj", translate(6,0,4)*scale(.14,.14,.14)*rotate_x(-M_PI/2), get_material(2, {.8,.8,.8}, {1,1,1}));
    meshes.emplace_back("modelos/wall.obj", scale(20, 20, 20), get_material(1, {0.4, 1.0, 0.4}));
    meshes.emplace_back("modelos/Wood Table/Old Wood Table.obj", translate(0,1.08,0));
    meshes.emplace_back("modelos/pose/pose.obj", translate(-6, 0, 4)*rotate_y(1)*scale(.05, .05, .05));
    meshes.emplace_back("modelos/train-toy-cartoon/train-toy-cartoon.obj", translate(0,0,6)*rotate_y(-2.3)*scale(120, 120, 120));
    meshes.emplace_back("modelos/box.obj", translate(-5, 1, 8), get_material(0, {5,5,0}));
    return meshes;
}

int main(){
    int w = 800;
    int h = 500;

    Scene scene{
        Image{w, h},
        Camera{0.1, 40*M_PI/180, w, h},
        get_meshes()
    };
	scene.camera.lookAt({10, 7, 20}, {0, 4, 0}, {0, 1, 0});
    scene.render(50);
    scene.image.savePNG("output.png");
}