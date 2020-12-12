#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <sstream>
#include "geometry.h"
#include "model.h"

float fog_density = 0.01;
Model::Model(const char *filename) : verts(), faces() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v[i];
            //Vec3f off(-14, 0, -8.5);
            Vec3f off(0., 0., 2.);
            v = v + off;
            verts.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            Vec3i f;
            int idx, cnt=0;
            iss >> trash;
            while (iss >> idx) {
                idx--; // in wavefront obj all indices start at 1, not zero
                f[cnt++] = idx;
            }
            if (3==cnt) faces.push_back(f);
        }
    }
    std::cerr << "# v# " << verts.size() << " f# "  << faces.size() << std::endl;

    Vec3f min, max;
    get_bbox(min, max);
}

// Moller and Trumbore
bool Model::ray_triangle_intersect(const int &fi, const Vec3f &orig, const Vec3f &dir, float &tnear) {
    Vec3f edge1 = point(vert(fi,1)) - point(vert(fi,0));
    Vec3f edge2 = point(vert(fi,2)) - point(vert(fi,0));
    Vec3f pvec = cross(dir, edge2);
    float det = edge1*pvec;
    if (det<1e-5) return false;

    Vec3f tvec = orig - point(vert(fi,0));
    float u = tvec*pvec;
    if (u < 0 || u > det) return false;

    Vec3f qvec = cross(tvec, edge1);
    float v = dir*qvec;
    if (v < 0 || u + v > det) return false;

    tnear = edge2*qvec * (1./det);
    return tnear>1e-5;
}


int Model::nverts() const {
    return (int)verts.size();
}

int Model::nfaces() const {
    return (int)faces.size();
}

void Model::get_bbox(Vec3f &min, Vec3f &max) {
    min = max = verts[0];
    for (int i=1; i<(int)verts.size(); ++i) {
        for (int j=0; j<3; j++) {
            min[j] = std::min(min[j], verts[i][j]);
            max[j] = std::max(max[j], verts[i][j]);
        }
    }
    std::cerr << "bbox: [" << min << " : " << max << "]" << std::endl;
}

const Vec3f &Model::point(int i) const {
    assert(i>=0 && i<nverts());
    return verts[i];
}

Vec3f &Model::point(int i) {
    assert(i>=0 && i<nverts());
    return verts[i];
}

int Model::vert(int fi, int li) const {
    assert(fi>=0 && fi<nfaces() && li>=0 && li<3);
    return faces[fi][li];
}

std::ostream& operator<<(std::ostream& out, Model &m) {
    for (int i=0; i<m.nverts(); i++) {
        out << "v " << m.point(i) << std::endl;
    }
    for (int i=0; i<m.nfaces(); i++) {
        out << "f ";
        for (int k=0; k<3; k++) {
            out << (m.vert(i,k)+1) << " ";
        }
        out << std::endl;
    }
    return out;
}


Model duck("/home/nastya/graphics/test/duck.obj");
Model triangle("/home/nastya/graphics/test/triangle.obj");
//Model triangle("/home/nastya/graphics/test/lowpolytree.obj");


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
    float bump_map;

    Material(): refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent(), bump_map(0.0) {}
    Material(const float r, const Vec4f &al, const Vec3f &color, const float spec, const float &bm): refractive_index(r), albedo(al), diffuse_color(color), specular_exponent(spec), bump_map(bm) {}
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

class Sphere1 {
public:
    Vec3f center;
    float sphere_radius;

    Sphere1(Vec3f cen, float rad) : center(cen), sphere_radius(rad) {}
    float signed_distance(const Vec3f &p) {
        Vec3f p_copy = p;
        return p_copy.norm() - sphere_radius;
    }
    bool sphere_trace(const Vec3f &orig, const Vec3f &dir, Vec3f &pos) {
        pos = orig;
        for (size_t i=0; i<128; i++) {
            float d = signed_distance(pos);
            if (d < 0) return true;
            pos = pos + dir*std::max(d*0.1f, .01f);
        }
        return false;
    }
};
Sphere1 test(Vec3f(-1, -1.5, -12), 2);

class Sphere {
public:
    Vec3f center;
    float radius;
    Material material;

    Sphere() {}
    Sphere (const Vec3f &cen, const float r, const Material &mat) : center(cen), radius(r), material(mat) {}

    /*float signed_distance(const Vec3f &p) const {
        Vec3f p_copy = p;
        return p_copy.norm() - radius;
    }
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &pos) const {
    pos = orig;
    for (size_t i=0; i<128; i++) {
        float d = signed_distance(pos);
        if (d < 0) return true;
        pos = pos + dir*std::max(d*0.1f, .01f);
    }
    return false;*/

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
            hit = orig + dir*dist_iter;
            N = (hit - (*q).center).normalize();
            material = (*q).material;
            if (material.bump_map > 1e-3) {
                float random_bias_x = (float(rand()) / RAND_MAX * 2 - 1) * material.bump_map;
                float random_bias_y = (float(rand()) / RAND_MAX * 2 - 1) * material.bump_map;
                float random_bias_z = (float(rand()) / RAND_MAX * 2 - 1) * material.bump_map;
                //std::cout << random_bias_x << std::endl << random_bias_y << std::endl << random_bias_z << std::endl;
                Vec3f bias(random_bias_x, random_bias_y, random_bias_z);
                N = (N + bias).normalize();
            }
        }
        ++q;
    }
    //return spheres_dist < 1000;

    /*float spheres_dist1 = std::numeric_limits<float>::max();
    float dist_iter = 0;
    Vec3f hit1;
    if (test.sphere_trace(orig, dir, hit1)) {
        spheres_dist1 = (hit-orig).norm();
        N = (hit - test.center).normalize();
    }*/

    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(dir.y) > 1e-3) { //проверяем чтобы не делить на ноль
        float d = -(orig.y+5.5)/dir.y; // доска находится в плоскости y = -4
        Vec3f pt = orig + dir*d;
        if (d>0 && fabs(pt.x)<40 && pt.z<-5 && pt.z>-200 && d<spheres_dist) {
            checkerboard_dist = d;
            hit = pt; //точка пересечения луча с доской
            N = Vec3f(0,1,0); //нормаль
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(.3, .3, .3) : Vec3f(.06, .04, .0);
            //material = Material(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3),   50.);
        }
    }
    float triangles_dist = std::numeric_limits<float>::max();
    /*for (int i = 0; i < triangle.nfaces(); i++) {
        float dist_iter;
        if (triangle.ray_triangle_intersect(i, orig, dir, dist_iter) && dist_iter < triangles_dist && dist_iter < spheres_dist) {
            triangles_dist = dist_iter;
            hit = orig + dir*dist_iter;
            Vec3f v0 = triangle.point(triangle.vert(i, 0));
            Vec3f v1 = triangle.point(triangle.vert(i, 1));
            Vec3f v2 = triangle.point(triangle.vert(i, 2));
            N = cross(v1-v0, v2-v0).normalize();
            //material = Material(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1),   10.);
            material = Material(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1000., 0.0);
        }
    }*/
    //float triangles_dist = std::numeric_limits<float>::max();
    /*for (int i = 0; i < duck.nfaces(); i++) {
        float dist_iter;
        if (duck.ray_triangle_intersect(i, orig, dir, dist_iter) && dist_iter < triangles_dist && dist_iter < spheres_dist) {
            triangles_dist = dist_iter;
            hit = orig + dir*dist_iter;
            Vec3f v0 = duck.point(duck.vert(i, 0));
            Vec3f v1 = duck.point(duck.vert(i, 1));
            Vec3f v2 = duck.point(duck.vert(i, 2));
            N = cross(v1-v0, v2-v0).normalize();
            //material = Material(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(1.0, 0.8, 0.02), 50., 0.0);
            //material = Material(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(1.0, 0.8, 0.02),   50.);
            //material = Material(1.0, Vec4f(0.8,  0.2, 0.1, 0.0), Vec3f(1.0, 0.8, 0.02),   10.);
            material = Material(1.5, Vec4f(0.3,  1.5, 0.2, 0.5), Vec3f(.24, .21, .09),  125., 0.0);
        }
    }*/

    /*if (num == 1) {
        float displacement = (sin(10*hit.x)*sin(10*hit.y)*sin(10*hit.z))*10;
        spheres_dist = spheres_dist - displacement;
        return std::min(spheres_dist, checkerboard_dist)<1000;
    } else {*/
        //float cur_min = std::min(triangles_dist, spheres_dist);
        return std::min(triangles_dist, std::min(spheres_dist, checkerboard_dist))<1000;
        //return std::min(spheres_dist, checkerboard_dist)<1000;
    //}

}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth=0) {
    Vec3f point, N;
    Material material;
    if (depth > 4 || !scene_intersect(orig, dir, spheres, point, N, material)) { //глубину рекурсии задаем здесь
        return Vec3f(0.98, 0.98, 0.98);
        //return Vec3f(0.12, 0.11, 0.37); // цвет фона
        //return Vec3f(0.1, 0.1, 0.1); // цвет фона
    }
    float z_dist = orig[2]-point[2];

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
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmp) && (shadow_pt-shadow_orig).norm() < light_dist) {
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

    //float fog_intensity = fog_density * z_dist / 60;
    float fog_intensity = fog_density * (1-exp(-z_dist/25));
    Vec3f fog_color(1.0, 1.0, 1.0);
    Vec3f object_color(diff_part + spec_part + reflect_part + refract_part);
    object_color[0] *= (1-fog_intensity);
    object_color[1] *= (1-fog_intensity);
    object_color[2] *= (1-fog_intensity);
    fog_color[0] *=fog_intensity;
    fog_color[1] *=fog_intensity;
    fog_color[2] *=fog_intensity;

    //return diff_part + spec_part + reflect_part + refract_part;
    /*if (num == 1) {
        float displacement = (sin(10*point.x)*sin(10*point.y)*sin(10*point.z) + 1.)/2.0;
        return (diff_part + spec_part + reflect_part + refract_part)*displacement;
    } else {*/
        return object_color + fog_color;
    //}

}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {

    const int width = 1024;
    const int height = 768;
    const int fov = M_PI/2;
    std::vector<Vec3f> framebuffer(width*height);

    int n1 = 2, n2 = 2;
    int num_samples = n1 * n2;

    /*#pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            Vec3f color(0.0, 0.0, 0.0);
            for (int i_sub = 0; i_sub < n1; i_sub++) {
                for (int j_sub = 0; j_sub < n2; j_sub++) {
                    float x_jitter = float(rand()) / RAND_MAX * 2.0 - 1;
                    float y_jitter = float(rand()) / RAND_MAX * 2.0 - 1;
                    float dir_x =  (i + x_jitter + 0.5) -  width/2.;
                    float dir_y = -(j + y_jitter + 0.5) + height/2.;    // this flips the image at the same time
                    float dir_z = -height/(2.*tan(fov/2.));
                    color = color + cast_ray(Vec3f(0,0,3), Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
                }
            }
            color[0] = color[0] / num_samples;
            color[1] = color[1] / num_samples;
            color[2] = color[2] / num_samples;
            framebuffer[i+j*width] = color;
        }
    }*/

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float dir_x =  (i + 0.5) -  width/2.;
            float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
            float dir_z = -height/(2.*tan(fov/2.));
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,3), Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
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
    Material ivory(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3), 50., 0.0);
    Material glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8), 125., 0.0);
    Material glass_bump(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8), 125., 0.02);
    Material red(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.62, 0.07, 0.06), 10., 0.2);
    //Material red(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1), 10., 0.0);
    Material mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425., 0.0);
    Material orange(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.53, 0.25, 0.03), 10., 0.1);
    Material purple(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.25, 0.02, 0.46), 10., 0.1);
    //Material green(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.02, 0.47, 0.05), 10., 0.2);
    std::vector<Sphere> spheres;
    /*spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, orange));
    spheres.push_back(Sphere(Vec3f(-1, -1.5, -12), 2, glass_bump));
    spheres.push_back(Sphere(Vec3f(1.5, -2.0, -30), 2, red));
    spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, mirror));
    spheres.push_back(Sphere(Vec3f(7, -2.5, -14), 1.5, glass));*/
    spheres.push_back(Sphere(Vec3f(1, -1.5, -15), 4, mirror));
    spheres.push_back(Sphere(Vec3f(-1, -4.3, -9.5), 1.2, orange));
    spheres.push_back(Sphere(Vec3f(6, -4, -12), 1.5, glass_bump));
    spheres.push_back(Sphere(Vec3f(1.0, -4.75, -7.7), 0.75, glass));
    //spheres.push_back(Sphere(Vec3f(5, -2.5, -5), 3, purple));
    std::vector<Light> lights;
    lights.push_back(Light(Vec3f(-20, 20, 20), 1.5));
    lights.push_back(Light(Vec3f(30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f(30, 20, 30), 1.7));
    //lights.push_back(Light(Vec3f(0, 0, 0), 2));
    render(spheres, lights);
    return 0;
}
