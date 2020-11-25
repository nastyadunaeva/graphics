#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

class Light {
public:
    Vec3f position;
    float intensity;
    Light(const Vec3f &pos, const float &in): position(pos), intensity(in) {};
};

class Material {
public:
    Vec2f albedo;
    Vec3f diffuse_color;
    float specular_exponent;

    Material(): albedo(1,0), diffuse_color(), specular_exponent() {}
    Material(const Vec2f &al, const Vec3f &color, const float &spec): albedo(al), diffuse_color(color), specular_exponent(spec) {}
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*(I*N)*2.f;
}

class Sphere {
public:
    Vec3f center;
    float radius;
    Material material;

    Sphere() {}
    //Sphere (const Vec3f &cen, const float &r) : center(cen), radius(r) {}
    Sphere (const Vec3f &cen, const float &r, const Material &mat) : center(cen), radius(r), material(mat) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres,
                     Vec3f &hit, Vec3f &N, Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    std::vector<Sphere>::const_iterator q = spheres.begin();
    while (q != spheres.end()) {
        float dist_iter;
        if ((*q).ray_intersect(orig, dir, dist_iter) && dist_iter < spheres_dist) {
            spheres_dist = dist_iter;
            hit = orig + dir*dist_iter; //hit, N
            N = (hit - (*q).center).normalize();
            material = (*q).material;
        }
        ++q;
    }
    return spheres_dist < 1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    Vec3f point, N;
    Material material;
    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return Vec3f(0.02, 0.8, 0.96); // цвет фона
    }
    //модель освещения Фонга
    float diffuse_light_intensity = 0;
    float specular_light_intensity = 0;
    std::vector<Light>::const_iterator q = lights.begin();
    while (q != lights.end()) {
        Vec3f light_dir = ((*q).position - point).normalize();
        float light_dist = ((*q).position - point).norm(); //расстояние от точки до света
        Vec3f shadow_orig = point;
        if (light_dir*N < 0) {shadow_orig = shadow_orig - N*1e-3;} // проверяем лежит ли точка в тени источника
        else {shadow_orig = shadow_orig + N*1e-3;} //немного сдвигаем точку в направлении нормали, чтобы не попасть туда же
        Vec3f shadow_pt, shadow_N;
        Material tmp;
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmp) && (shadow_pt-shadow_orig).norm() < light_dist) {
            ++q;
            continue; //если луч точка-источник пересекает объекты сцены, то игнорируем источник
        }

        diffuse_light_intensity += (*q).intensity * std::max(0.f, light_dir*N); //интенсивность зависит от угла между нормалью и направлением света
        specular_light_intensity += powf(std::max(0.f, reflect(light_dir, N)*dir), material.specular_exponent)*(*q).intensity; //отсвет обратно пропорционален углу между направлением взгляда и направлением отраженного свет
        ++q;
    }
    Vec3f v = Vec3f(1.0, 1.0, 1.0)*specular_light_intensity*material.albedo[1];
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + v;
}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {

    const int width = 1024;
    const int height = 768;
    const int fov = M_PI/2;
    std::vector<Vec3f> framebuffer(width*height);

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), dir, spheres, lights);
        }
    }

    std::ofstream ofs;
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        //на этапе отражения - зачем?
        float mx = std::max(framebuffer[i][0], std::max(framebuffer[i][1], framebuffer[i][2]));
        if (mx > 1) {
            framebuffer[i] = framebuffer[i] * (1./mx);
        }
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}



int main() {
    Material ivory(Vec2f(0.6, 0.3), Vec3f(0.4, 0.4, 0.3), 50.0);
    Material red(Vec2f(0.9, 0.1), Vec3f(0.3, 0.1, 0.1), 10.0);
    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));
    spheres.push_back(Sphere(Vec3f(-1, -1.5, -12), 2, red));
    spheres.push_back(Sphere(Vec3f(1.5, -0.5, -18), 3, red));
    spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, ivory));
    std::vector<Light> lights;
    lights.push_back(Light(Vec3f(-20, 20, 20), 1.5));
    lights.push_back(Light(Vec3f(30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f(30, 20, 30), 1.7));
    render(spheres, lights);
    return 0;
}
