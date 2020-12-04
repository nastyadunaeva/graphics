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
    Light(const Vec3f &pos, const float in): position(pos), intensity(in) {};
};

class Material {
public:
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;

    Material(): refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
    Material(const float r, const Vec4f &al, const Vec3f &color, const float spec): refractive_index(r), albedo(al), diffuse_color(color), specular_exponent(spec) {}
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*(I*N)*2.f;
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float n_second, const float n_first=1.f) {
    float cos_inc = - std::max(-1.f, std::min(1.f, I*N));
    if (cos_inc < 0) { //если луч внутри объекта - все наоборот
        return refract(I, -N, n_first, n_second);
    }
    float n_rel = n_first / n_second;
    float rad = 1 - n_rel*n_rel*(1 - cos_inc*cos_inc);
    if (rad < 0) {
        return Vec3f(1,0,0);
    } else {
        return I*n_rel + N*(n_rel*cos_inc - sqrtf(rad));
    }
}

class Sphere {
public:
    Vec3f center;
    float radius;
    Material material;

    Sphere() {}
    Sphere (const Vec3f &cen, const float r, const Material &mat) : center(cen), radius(r), material(mat) {}

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
                     Vec3f &hit, Vec3f &N, Material &material, int &num) {
    float spheres_dist = std::numeric_limits<float>::max();
    std::vector<Sphere>::const_iterator q = spheres.begin();
    int i = 1;
    while (q != spheres.end()) {
        float dist_iter;
        if ((*q).ray_intersect(orig, dir, dist_iter) && dist_iter < spheres_dist) {
            spheres_dist = dist_iter;
            hit = orig + dir*dist_iter; //hit, N
            N = (hit - (*q).center).normalize();
            material = (*q).material;
            num = i;
        }
        ++q;
        i++;
    }
    //return spheres_dist < 1000;

    float checkerboard_dist  =std::numeric_limits<float>::max();
    if (fabs(dir.y) > 1e-3) { //проверяем чтобы не делить на ноль
        float d = -(orig.y+4)/dir.y; // доска находится в плоскости y = -4
        Vec3f pt = orig + dir*d;
        if (d>0 && fabs(pt.x)<20 && pt.z<-5 && pt.z>-30 && d<spheres_dist) {
            checkerboard_dist = d;
            hit = pt; //точка пересечения луча с доской
            N = Vec3f(0,1,0); //нормаль
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(.3, .3, .3) : Vec3f(.98, .42, .02);
        }
    }
    return std::min(spheres_dist, checkerboard_dist)<1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth=0) {
    Vec3f point, N;
    Material material;
    int num;
    if (depth > 2 || !scene_intersect(orig, dir, spheres, point, N, material, num)) { //глубину рекурсии задаем здесь
        return Vec3f(0.1, 0.1, 0.1); // цвет фона
        //return Vec3f(0.1, 0.1, 0.1); // цвет фона
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();

    Vec3f reflect_orig = point;
    if (reflect_dir*N < 0) { //немного сдвигаем точку в направлении нормали, чтобы не попасть туда же
        reflect_orig = reflect_orig - N*1e-3;
    } else {
        reflect_orig = reflect_orig + N*1e-3;
    }
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;

    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, depth + 1);
    Vec3f refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, depth + 1);

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
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmp, num) && (shadow_pt-shadow_orig).norm() < light_dist) {
            ++q;
            continue; //если луч точка-источник пересекает объекты сцены, то игнорируем источник
        }

        diffuse_light_intensity += (*q).intensity * std::max(0.f, light_dir*N); //интенсивность зависит от угла между нормалью и направлением света
        specular_light_intensity += powf(std::max(0.f, reflect(light_dir, N)*dir), material.specular_exponent)*(*q).intensity; //отсвет обратно пропорционален углу между направлением взгляда и направлением отраженного свет
        ++q;
    }
    Vec3f diff_part = material.diffuse_color * diffuse_light_intensity * material.albedo[0];
    Vec3f spec_part = Vec3f(1.0, 1.0, 1.0) * specular_light_intensity * material.albedo[1];
    Vec3f reflect_part = reflect_color * material.albedo[2];
    Vec3f refract_part = refract_color * material.albedo[3];
    if (num == 1) {
        float displacement = (sin(10*point.x)*sin(10*point.y)*sin(10*point.z) + 1.)/2.0;
        return (diff_part + spec_part + reflect_part + refract_part)*displacement;
    } else {
        return diff_part + spec_part + reflect_part + refract_part;
    }

}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {

    const int width = 1024;
    const int height = 768;
    const int fov = M_PI/2;
    std::vector<Vec3f> framebuffer(width*height);

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float dir_x =  (i + 0.5) -  width/2.;
            float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
            float dir_z = -height/(2.*tan(fov/2.));
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
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
    Material      ivory(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3),   50.);
    Material      glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125.);
    Material red(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);
    //Material      black(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.1, 0.1, 0.1),   50.);
    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));
    spheres.push_back(Sphere(Vec3f(-1, -1.5, -12), 2, glass));
    spheres.push_back(Sphere(Vec3f(1.5, -0.5, -18), 3, red));
    spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, mirror));
    //Vec3f li(-4, -1.5, -12);
    //spheres.push_back(Sphere(li, 0.1, glass));

    //spheres.push_back(Sphere(Vec3f(0, 0, -1000), 900, black));
    std::vector<Light> lights;
    lights.push_back(Light(Vec3f(-20, 20, 20), 1.5));
    lights.push_back(Light(Vec3f(30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f(30, 20, 30), 1.7));
    //lights.push_back(Light(li, 10));
    render(spheres, lights);
    return 0;
}
