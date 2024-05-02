#ifndef CAMERA_H
#define CAMERA_H
#include "hittable.h"
#include "rtweekend.h"
#include <fstream>
#include "material.h"
class camera {
public:

    double aspect_ratio = 1.0; 
    int    image_width = 100;
    int samples_per_pixel = 10;
    int max_depth = 10;

    double vfov = 90;
    point3 lookfrom = point3(0, 0, 0);   
    point3 lookat = point3(0, 0, -1);  
    vec3   vup = vec3(0, 1, 0);

    double defocus_angle = 0; 
    double focus_dist = 10;   
    void render(const hittable& world, int num_threads) {
        initialize();

        std::ofstream my_image("image.ppm");
        my_image << "P3\n" << image_width << ' ' << image_height << "\n255\n";

        std::vector<std::thread> threads(num_threads);
        std::vector<std::vector<color>> partial_outputs(num_threads);

        const int rows_per_thread = image_height / num_threads;

        for (int i = 0; i < num_threads; ++i) {
            int start_row = i * rows_per_thread;
            int end_row = (i == num_threads - 1) ? image_height : (i + 1) * rows_per_thread;
            partial_outputs[i].resize(image_width * (end_row - start_row));
            threads[i] = std::thread(&camera::render_chunk, this, start_row, end_row, i, std::ref(partial_outputs[i]), std::cref(world));
        }

        for (auto& t : threads) {
            t.join();
        }

        for (const auto& partial_output : partial_outputs) {
            for (const auto& pixel_color : partial_output) {
                write_color(my_image, pixel_samples_scale * pixel_color);
            }
        }

        std::cout << "Done.\n";
    }

    void render_chunk(int start_row, int end_row, int thread_index, std::vector<color>& output, const hittable& world) {
    for (int j = start_row; j < end_row; j++) {
        std::clog << "\rRemaining scanlines: " << (end_row - j) << std::flush;
        for (int i = 0; i < image_width; i++) {
            color pixel_color(0, 0, 0);
            for (int sample = 0; sample < samples_per_pixel; sample++) {
                ray r = get_ray(i, j);
                pixel_color += ray_color(r, max_depth, world);
            }
            output[(j - start_row) * image_width + i] = pixel_color;
        }
    }
    std::clog << "\rDone.                 \n";
}
   private:
       int    image_height;  
       double pixel_samples_scale;
       point3 center;         
       point3 pixel00_loc;    
       vec3   pixel_delta_u;  
       vec3   pixel_delta_v;  
       vec3 u, v, w;

       vec3   defocus_disk_u;  
       vec3   defocus_disk_v; 
       void initialize()
       {  
           image_height = int(image_width / aspect_ratio);
           image_height = (image_height < 1) ? 1 : image_height;
           pixel_samples_scale = 1.0 / samples_per_pixel;
           center = lookfrom;
           //Camera
          
           auto theta = degrees_to_radians(vfov);
           auto h = tan(theta / 2);
           auto viewport_height = 2*h*focus_dist;
           auto viewport_width = viewport_height * (double(image_width) / image_height);
      

           // Calculating the u,v,w unit basis vectors for  camera coordinate frames.
           w = unit_vector(lookfrom - lookat);
           u = unit_vector(cross(vup,w));
           v = cross(w, u);
           //calculating vector across viewport edges for horizantal and downwards
           auto viewport_u = viewport_width*u;
           auto viewport_v = viewport_height*-v;

           //calculating the horizontal and delta vectors for pixel to pixel
           pixel_delta_u = viewport_u / image_width;
           pixel_delta_v = viewport_v / image_height;

           //claculate the coordinate of upper left pixel
           auto viewport_upper_left = center - (focus_dist*w) - viewport_u / 2 - viewport_v / 2;
            pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

            auto defocus_radius = focus_dist * tan(degrees_to_radians(defocus_angle / 2));
            defocus_disk_u = u * defocus_radius;
            defocus_disk_v = v * defocus_radius;

       }

       ray get_ray(int i, int j) const {
           auto offset = sample_square();
           auto pixel_sample = pixel00_loc + ((i + offset.x()) * pixel_delta_u) + ((j + offset.y()) * pixel_delta_v);

           auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
           auto ray_direction = pixel_sample - ray_origin;

           return ray(ray_origin, ray_direction);
       }

       vec3 sample_square() const {
           return vec3(random_double() - 0.5, random_double() - 0.5,0);
       }
       point3 defocus_disk_sample() const 
       {
          auto p = random_in_unit_disk();
           return center + (p[0] * defocus_disk_u) + (p[1] * defocus_disk_v);
       }
       color ray_color(const ray& r, int depth, const hittable& world) const {
           if (depth <= 0)
               return color(0, 0, 0);

           hit_record rec;

           if (world.hit(r, interval(0.001, infinity), rec)) {
               ray scattered;
               color attenuation;
               if (rec.mat->scatter(r, rec, attenuation, scattered))
                   return attenuation * ray_color(scattered, depth - 1, world);
               return color(0, 0, 0);
           }

           vec3 unit_direction = unit_vector(r.direction());
           auto a = 0.5 * (unit_direction.y() + 1.0);
           return (1.0 - a) * color(1.0, 1.0, 1.0) + a * color(0.5, 0.7, 1.0);
       }
};
#endif