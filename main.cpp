#include <iostream>
#include <vector>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

using namespace std;

Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& scene,
               int trace_depth);

Color Shade(const std::vector<LightSource>& light_sources,
            const Hittable& hittable_collection,
            const HitRecord& hit_record,
            int trace_depth) {
    // TODO: Add your code here.
    Color color(0.f, 0.f, 0.f);

    color = hit_record.material.ambient;

    for(const auto& light_source : light_sources) {
        glm::vec3 shadow_ray;
        for (int i = 0; i < 3; i++) {
            shadow_ray[i] = light_source.position[i] - hit_record.position[i];
        }
        if (glm::dot(hit_record.normal, shadow_ray) > 0) {
//            if (hittable_collection.Hit(shadow_ray, &hit_record) == false) {
                Color diffuse = (hit_record.material.k_d * hit_record.material.diffuse)
                                * glm::dot(hit_record.normal, shadow_ray);
                glm::vec3 out_direction;
                for (int i = 0; i < 3; i++) {
                    out_direction[i] = hit_record.in_direction[i] * -1.0;
                }
                Color specular = (hit_record.material.k_s * hit_record.material.specular)
                                * glm::pow((glm::dot(hit_record.reflection, out_direction)), hit_record.material.sh);
                color += light_source.intensity * (diffuse + specular);
//            }
        }
    }
    if (trace_depth < kMaxTraceDepth) {
        glm::vec3 reflected_ray;
        for (int i = 0; i < 3; i++) {
            reflected_ray[i] = hit_record.position[i] + hit_record.reflection[i];
        }
//        if (hit_record.material.k_a > 0) {
//            Color r_color = TraceRay(reflected_ray, light_sources, hittable_collection, trace_depth+1);
//            color += r_color * hit_record.material.k_a;
//        }
//        if (hit_record.material.k_d > 0) {
//            Color r_color = TraceRay(reflected_ray, light_sources, hittable_collection, trace_depth+1);
//            color += r_color * hit_record.material.k_d;
//        }
        if (hit_record.material.k_s > 0) {
            Color r_color = TraceRay(reflected_ray, light_sources, hittable_collection, trace_depth+1);
            color += r_color * hit_record.material.k_s;
        }
    }

    for (int i = 0; i < 3; i++) {
        if (color[i] > 1.0) {
            color[i] = 1.0;
        }
    }

    return color;
}

Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& hittable_collection,
               int trace_depth) {
    // TODO: Add your code here.
    HitRecord record;
    Color color(0.0f, 0.0f, 0.0f);

    if (hittable_collection.Hit(ray, &record)) {
        color = Shade(light_sources, hittable_collection, record, trace_depth);
    }
    
    return color;
}

int main() {
    // TODO: Set your workdir (absolute path) here.
//    const std::string work_dir("C:/.../Graphics_PA3_Release/");
    const std::string work_dir("/Users/carl/Documents/Programs/C++/COMP3271/Assignments/3/Graphics_PA3_Release/");


    // Construct scene
    Scene scene(work_dir, "scene/sphere.toml");
    const Camera& camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;

    std::vector<unsigned char> image(width * height * 4, 0);

    float progress = 0.f;

    // Traverse all pixels
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Color color(0.f, 0.f, 0.f);
            int count = 0;
            for (float bias_x = 0.25f; bias_x < 1.f; bias_x += .5f) {
                for (float bias_y = 0.25f; bias_y < 1.f; bias_y += .5f) {
                    Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
                    color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
                    count++;
                }
            }
            color /= float(count);
            int idx = 4 * ((height - y - 1) * width + x);
            for (int i = 0; i < 3; i++) {
                image[idx + i] = (uint8_t) (glm::min(color[i], 1.f - 1e-5f) * 256.f);
            }
            image[idx + 3] = 255;

            float curr_progress = float(x * height + y) / float(height * width);
            if (curr_progress > progress + 0.05f) {
                progress += 0.05f;
                std::cout << "Progress: " << progress << std::endl;
            }
        }
    }

    // Save result as png file
    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    lodepng::save_file(png, work_dir + "output.png");
}
