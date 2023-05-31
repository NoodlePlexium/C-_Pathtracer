#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <cstdio>
#include <limits>

#include "Colour.h"
#include "Vector.h"
#include "Camera.h"
#include "Object.h"
#include "Light.h"
#include "Plane.h"
#include "Sphere.h"
#include "Ray.h"



class Renderer {
    Camera camera;
    Plane plane;
    Sphere sphere;
    Light light;

    double Pi;
    Colour** pixels;

public:
    Renderer(const Camera _camera, Plane _plane, Sphere _sphere, Light _light)
    : camera(_camera), plane(_plane), sphere(_sphere), light(_light)
    {
        this->Pi = 3.14159265358979323846;
        InitializePixels();
    }

    void Render(){



        // Iterate over each pixel /////////////////////////////////
        for (int y = 0; y < camera.getHeight(); y++) {
        for (int x = 0; x < camera.getWidth(); x++) {

            // Code runs per Pixel ///////
            pixels[y][x] = PerPixel(x, y);

        // End /////////////////////////////////////////////////////
        }
        }


        // Save Image
        saveBMP("Render.bmp");
    }

    Colour PerPixel(int x, int y){

        double ndcX = (2.0 * (x + 0.5) / camera.getWidth()) - 1.0;
        double ndcY = 1.0 - (2.0 * (y + 0.5) / camera.getHeight());

        double aspectRatio = camera.getWidth() / static_cast<double>(camera.getHeight());

        double fovY = camera.getFOV() * (Pi / 180.0);
        double fovX = 2.0 * std::atan(std::tan(fovY / 2.0) * aspectRatio);
        double angleX = ndcX * fovX / 2.0;
        double angleY = ndcY * fovY / 2.0;

        // Calculate the ray direction for the current pixel
        Vector rayDirection = (camera.Forward() + camera.Right() * std::tan(angleX) + camera.Down() * std::tan(angleY)).normalized();


        // Cast a ray to each object in the scene //////////////////////////
        Ray ray(camera.Position(), rayDirection);
        RayHit planeHit = plane.CalculateRayHit(ray);
        RayHit sphereHit = sphere.CalculateRayHit(ray, light.getPosition());

        // Determine which object is closest ///////////////////////////////
        RayHit hit;
        if (sphereHit.getHitDistance() < planeHit.getHitDistance()) hit = sphereHit;
        else hit = planeHit;



        // Shading ///////////////////////////////////////////////////
        // Calculate the brightness based on the cosine of the angle
        double brightness = std::max(0.0, hit.getHitNormal().dotProduct(light.getDirection()));

        // Return colour
        return (hit.getMeshColour() * brightness);
    }

    void InitializePixels(){
        pixels = new Colour*[camera.getHeight()];
        for (int i = 0; i < camera.getHeight(); i++) {
            pixels[i] = new Colour[camera.getWidth()];
        }
    }

    void saveBMP(const char* filename) {
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            // Handle file opening error
            return;
        }

        int width = camera.getWidth();
        int height = camera.getHeight();
        int dpi = 72;

        int s = width * height;
        int filesize = 54 + 3 * s;

        double factor = 39.375;
        int m = static_cast<int>(factor);
        int ppm = dpi * m;

        unsigned char bmpfileheader[14] = {'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0};
        unsigned char bmpinfoheader[40] = {40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0};

        bmpfileheader[2] = static_cast<unsigned char>(filesize);
        bmpfileheader[3] = static_cast<unsigned char>(filesize >> 8);
        bmpfileheader[4] = static_cast<unsigned char>(filesize >> 16);
        bmpfileheader[5] = static_cast<unsigned char>(filesize >> 24);

        bmpinfoheader[4] = static_cast<unsigned char>(width);
        bmpinfoheader[5] = static_cast<unsigned char>(width >> 8);
        bmpinfoheader[6] = static_cast<unsigned char>(width >> 16);
        bmpinfoheader[7] = static_cast<unsigned char>(width >> 24);

        bmpinfoheader[8] = static_cast<unsigned char>(height);
        bmpinfoheader[9] = static_cast<unsigned char>(height >> 8);
        bmpinfoheader[10] = static_cast<unsigned char>(height >> 16);
        bmpinfoheader[11] = static_cast<unsigned char>(height >> 24);

        bmpinfoheader[21] = static_cast<unsigned char>(3 * s);
        bmpinfoheader[22] = static_cast<unsigned char>(3 * s >> 8);
        bmpinfoheader[23] = static_cast<unsigned char>(3 * s >> 16);
        bmpinfoheader[24] = static_cast<unsigned char>(3 * s >> 24);

        bmpinfoheader[25] = static_cast<unsigned char>(ppm);
        bmpinfoheader[26] = static_cast<unsigned char>(ppm >> 8);
        bmpinfoheader[27] = static_cast<unsigned char>(ppm >> 16);
        bmpinfoheader[28] = static_cast<unsigned char>(ppm >> 24);

        bmpinfoheader[29] = static_cast<unsigned char>(ppm);
        bmpinfoheader[30] = static_cast<unsigned char>(ppm >> 8);
        bmpinfoheader[31] = static_cast<unsigned char>(ppm >> 16);
        bmpinfoheader[32] = static_cast<unsigned char>(ppm >> 24);

        file.write(reinterpret_cast<char*>(bmpfileheader), 14);
        file.write(reinterpret_cast<char*>(bmpinfoheader), 40);

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                Colour pixel = pixels[i][j];
                unsigned char color[3] = {static_cast<unsigned char>(pixel.getBlue() * 255),
                                          static_cast<unsigned char>(pixel.getGreen() * 255),
                                          static_cast<unsigned char>(pixel.getRed() * 255)};
                file.write(reinterpret_cast<char*>(color), 3);
            }
        }

        file.close();
    }
};



int main() {
    Vector pos(0, -0.5, 0);
    Vector dir(0, 0, 1);
    Vector right(1, 0, 0);
    Vector down(0, -1, 0);


    // Camera setup
    // Set camera rotation
    Vector angles(0, 0, 0);  // Rotation angles in degrees

    // Apply rotation
    dir = dir.rotate(angles);
    right = right.rotate(angles);

    int width = 500;
    int height = 350;
    double fov = 70;

    Camera cam(pos, dir, right, down, width, height, fov);

    // Create a plane
    Plane plane(Vector(0, -2, 5), Vector(0, 0, 0), Vector(5, 5, 5));
    plane.setColour(Colour(0.1, 0.8, 0.3, 1));

    // Create a sphere
    Sphere sphere(Vector(0, -1.5, 5), 0.5, Colour(1, 0.2, 0.2, 1));

    // Create a light
    Light light(Vector(2, 10, -5), Vector(-1, -1, 1), Colour(1,1,1,1));


    Renderer renderer(cam, plane, sphere, light);
    renderer.Render();

    return 0;
}
