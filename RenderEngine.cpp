#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <cstdio>
#include <limits>
#include <vector>
#include <random>

#include "Colour.h"
#include "Vector.h"
#include "Camera.h"
#include "Object.h"
#include "Light.h"
#include "Plane.h"
#include "Sphere.h"
#include "Ray.h"
#include "Helper.h"



class Renderer {
    Camera camera;
    std::vector<Object*> objects;
    Light light;

    double Pi;
    Colour** pixels;

public:
    Renderer(const Camera _camera, std::vector<Object*> _objects, Light _light)
    : camera(_camera), objects(_objects), light(_light)
    {
        this->Pi = 3.14159265358979323846;
        InitializePixels();
    }

    void Render(){



        // Iterate over each pixel /////////////////////////////////
        for (int y = 0; y < camera.getHeight(); y++) {
        for (int x = 0; x < camera.getWidth(); x++) {

        // Code runs per Pixel ///////
        pixels[y][x] = PerPixel(x,y);
        }
        }


        // Save Image
        saveBMP("Render4.bmp");
    }

    Colour PerPixel(int x, int y){

    	// Shading ///////////////////////////////////////////////////
        Colour skyColour(0.55, 0.69, 0.84);
        Colour directIllumination; // Set to sky colour by default
        Colour indirectIllumination; // Illumination from bounce light

        int samples = 10;
		for (int s = 0; s < samples; s++){

			// Initialize sample colours /////////////////
           	Colour sampleDirect = skyColour; 
           	Colour sampleIndirect;

           	// First ray from camera ///////////////////////////
	        Ray ray(camera.Position(), DirFromPixel(x, y));
	        RayHit hit = CastRay(ray);

	        // Ray hits scene object
	        if (hit.DidHit()){
	            double brightness = std::max(0.0, hit.getHitNormal().dotProduct(light.getDirection())); // Based on angle to light source
	            sampleDirect = hit.getMeshColour() * light.getColour() * brightness;
	            Ray shadowRay(hit.getHitPosition()-hit.getHitNormal()*0.00001, light.getDirection()*-1); // Check if point is in shadow
	        	RayHit shadowHit = CastRay(shadowRay);
	        	if (shadowHit.DidHit()){sampleDirect *= Colour(0,0,0,1);} 
	        }


	        // Perform ray tracing for reflections
	        if (hit.DidHit()){

	        	Vector reflectionVector = hit.getHitNormal() * 2 * ray.getRayDirection().dotProduct(hit.getHitNormal());
    			Vector diffuseDir = Helper().RandomHemisphereDirection(hit.getHitNormal());
    			Vector specularDir = ray.getRayDirection() - reflectionVector;

    			ray.setOrigin(hit.getHitPosition() - hit.getHitNormal() * 0.00001);
    			ray.setDirection(diffuseDir * hit.getMeshMaterial().getRoughness() +  specularDir * (1-hit.getMeshMaterial().getRoughness()));

				int bounces = 4;
	          	for (int i=0; i < bounces; i++){

	                // Cast a reflection ray
	                hit = CastRay(ray);

	                if (!hit.DidHit()){ // Break the loop nothing is hit (sky)
	                    //sampleIndirect.blendAdd(skyColour / (bounces+1));
	                    break;  
	                } 

	               	// Calculate the indirect colour
	                double angleIntensity = std::max(0.0, hit.getHitNormal().dotProduct(light.getDirection())); // Based on angle to light source
	                Colour indirectColour = hit.getMeshColour() * angleIntensity;

	                // Shadow Check
	                Ray shadowRay(hit.getHitPosition()-hit.getHitNormal()*0.00001, light.getDirection()*-1);
	        		RayHit shadowHit = CastRay(shadowRay);
	        		if (shadowHit.DidHit()){indirectColour = Colour(0,0,0,1);} 

	                // Blend the colour
	                sampleIndirect.blendAdd(indirectColour / (bounces));

	                // Set origin and direction for next bounce
	                Vector reflectionVector = hit.getHitNormal() * 2 * ray.getRayDirection().dotProduct(hit.getHitNormal());
    				Vector diffuseDir = Helper().RandomHemisphereDirection(hit.getHitNormal());
    				Vector specularDir = ray.getRayDirection() - reflectionVector;

    				ray.setOrigin(hit.getHitPosition() - hit.getHitNormal() * 0.00001);
    				ray.setDirection(diffuseDir * hit.getMeshMaterial().getRoughness() +  specularDir * (1-hit.getMeshMaterial().getRoughness()));
	            }
	        }

	        directIllumination.blendAdd(sampleDirect / samples);
	        indirectIllumination.blendAdd(sampleIndirect / samples);
        }

        return directIllumination + indirectIllumination;
    }

    Vector DirFromPixel(int x, int y){
    	double ndcX = (2.0 * (x + 0.5) / camera.getWidth()) - 1.0;
        double ndcY = 1.0 - (2.0 * (y + 0.5) / camera.getHeight());

        double aspectRatio = camera.getWidth() / static_cast<double>(camera.getHeight());

        double fovY = camera.getFOV() * (Pi / 180.0);
        double fovX = 2.0 * std::atan(std::tan(fovY / 2.0) * aspectRatio);
        double angleX = ndcX * fovX / 2.0;
        double angleY = ndcY * fovY / 2.0;

        // Calculate the ray direction for the current pixel
        return (camera.Forward() + camera.Right() * std::tan(angleX) + camera.Down() * std::tan(angleY)).normalized();
    }


    RayHit CastRay(Ray ray){
        RayHit hit;

        // Cast ray for every object in the scene and get the closest
        for (int j=0; j < objects.size(); j++){
            RayHit newHit = objects[j]->CalculateRayHit(ray);
            if (newHit.getHitDistance() < hit.getHitDistance()) {
                hit = newHit;
            }
        }

        return hit;
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
    Vector pos(0, 0, 0);
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
    int height = 500;
    double fov = 70;

    // Create a camera
    Camera cam(pos, dir, right, down, width, height, fov);


    // Create Scene objects /////////////////////////////////////////////////////////////
    std::vector<Object*> objects;

    Material mat1 = Material(0.1, 0, 0, Colour(0.8, 0.8, 0.8, 1));

    Material groundMat = Material(0, 0, 0, Colour(0.9, 0.9, 0.9, 1));
    Material BWMat = Material(0.5, 0, 0, Colour(0.2, 0.8, 0.2, 1));
    Material LWMat = Material(0.1, 0, 0, Colour(0.8, 0.2, 0.2, 1));
    Material RWMat = Material(0.1, 0, 0, Colour(0.2, 0.2, 0.8, 1));



    Sphere* sphere = new Sphere(Vector(0, -1.5, 5), Vector(0, 0, 0), Vector(1, 1, 1), 0.5, mat1); // Create a sphere
    Sphere* sphere2 = new Sphere(Vector(-1.1, -1.5, 6), Vector(0, 0, 0), Vector(0.5, 0.5, 0.5), 0.5, mat1); // Create a sphere
    Sphere* sphere3 = new Sphere(Vector(1.2, -1.5, 5), Vector(0, 0, 0), Vector(0.5, 0.5, 0.5), 0.5, mat1); // Create a sphere

    Plane* ground =  new Plane(Vector(0, -2, 5), Vector(0, 0, 0), Vector(5, 5, 5), groundMat); // Create a plane
    Plane* Right_Wall =  new Plane(Vector(2.5, -2, 5), Vector(0, 0, 90), Vector(5, 5, 5), RWMat); // Create a plane
    Plane* Left_Wall =  new Plane(Vector(-2.5, -2, 5), Vector(0, 0, -90), Vector(5, 5, 5), LWMat); // Create a plane
    Plane* Back_Wall =  new Plane(Vector(0, -2, 7.5), Vector(90, 180, 0), Vector(5, 5, 5), BWMat); // Create a plane

    objects.push_back(sphere);
    objects.push_back(sphere2);
    objects.push_back(sphere3);
    objects.push_back(ground);
    objects.push_back(Right_Wall);
    objects.push_back(Left_Wall);
    objects.push_back(Back_Wall);


    // Create a light
    Light light(Vector(2, 10, -5), Vector(0.5, -1, 1), Colour(1,1,1,1));


    Renderer renderer(cam, objects, light);
    renderer.Render();

    return 0;
}

