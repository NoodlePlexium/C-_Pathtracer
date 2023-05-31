#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <cstdio>

class Vector {
    double x, y, z;

public:
    Vector(double _x = 0.0, double _y = 0.0, double _z = 0.0) : x(_x), y(_y), z(_z) {}

    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }

    double magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vector normalized() const {
        double mag = magnitude();
        return Vector(x / mag, y / mag, z / mag);
    }

    Vector negative() const {
        return Vector(-x, -y, -z);
    }

    double dotProduct(const Vector& v) const {
        return x * v.getX() + y * v.getY() + z * v.getZ();
    }

    Vector crossProduct(const Vector& v) const {
        return Vector(
            y * v.getZ() - z * v.getY(),
            z * v.getX() - x * v.getZ(),
            x * v.getY() - y * v.getX()
        );
    }

    Vector operator+(const Vector& v) const {
        return Vector(x + v.getX(), y + v.getY(), z + v.getZ());
    }

    Vector operator-(const Vector& v) const {
        return Vector(x - v.getX(), y - v.getY(), z - v.getZ());
    }

    Vector operator*(double scalar) const {
        return Vector(x * scalar, y * scalar, z * scalar);
    }

    Vector rotate(const Vector& angles) const {

        // PI constant
        double Pi = 3.14159265358979323846;

        // Convert rotation angles from degrees to radians
        double angleX = angles.getX() * (Pi / 180.0);
        double angleY = angles.getY() * (Pi  / 180.0);
        double angleZ = angles.getZ() * (Pi / 180.0);

        // Apply rotations around X, Y, and Z axes
        Vector rotatedVector = *this;

        // Rotation around X-axis
        double cosX = std::cos(angleX);
        double sinX = std::sin(angleX);
        double rotatedY = rotatedVector.getY() * cosX - rotatedVector.getZ() * sinX;
        double rotatedZ = rotatedVector.getY() * sinX + rotatedVector.getZ() * cosX;
        rotatedVector = Vector(rotatedVector.getX(), rotatedY, rotatedZ);

        // Rotation around Y-axis
        double cosY = std::cos(angleY);
        double sinY = std::sin(angleY);
        double rotatedZ2 = rotatedVector.getZ() * cosY - rotatedVector.getX() * sinY;
        double rotatedX = rotatedVector.getZ() * sinY + rotatedVector.getX() * cosY;
        rotatedVector = Vector(rotatedX, rotatedVector.getY(), rotatedZ2);

        // Rotation around Z-axis
        double cosZ = std::cos(angleZ);
        double sinZ = std::sin(angleZ);
        double rotatedX2 = rotatedVector.getX() * cosZ - rotatedVector.getY() * sinZ;
        double rotatedY2 = rotatedVector.getX() * sinZ + rotatedVector.getY() * cosZ;
        rotatedVector = Vector(rotatedX2, rotatedY2, rotatedVector.getZ());

        return rotatedVector;
    }
};

class Colour {
    double red, green, blue, alpha;

public:
    Colour(double _red=0.0, double _green=0.0, double _blue=0.0, double _alpha=0.0)
        : red(_red), green(_green), blue(_blue), alpha(_alpha) {}

    double getRed() const { return red; }
    double getGreen() const { return green; }
    double getBlue() const { return blue; }
    double getAlpha() const { return alpha; }

    void setColourRed(double redValue) { red = redValue; }
    void setColourGreen(double greenValue) { green = greenValue; }
    void setColourBlue(double blueValue) { blue = blueValue; }
    void setColourAlpha(double alphaValue) { alpha = alphaValue; }

    // Add two colors
    Colour add(const Colour other) const {
        double newRed = red + other.red;
        double newGreen = green + other.green;
        double newBlue = blue + other.blue;
        double newAlpha = alpha + other.alpha;
        return Colour(newRed, newGreen, newBlue, newAlpha);
    }

    // Multiply the color by a scalar value
    Colour multiply(double scalar) const {
        double newRed = red * scalar;
        double newGreen = green * scalar;
        double newBlue = blue * scalar;
        double newAlpha = alpha * scalar;
        return Colour(newRed, newGreen, newBlue, newAlpha);
    }
};

class Camera {
    Vector campos, camdir, camright, camdown;
    int width, height;
    double fov;

public:
    Camera(const Vector& _pos, const Vector& _dir, const Vector& _right, 
        const Vector& _down, int _width, int _height, double _fov)
        : campos(_pos), camdir(_dir), camright(_right), camdown(_down),
          width(_width), height(_height), fov(_fov) {}

    Vector Position() const { return campos; }
    Vector Forward() const { return camdir; }
    Vector Right() const { return camright; }
    Vector Down() const { return camdown; }

    int getWidth() const { return width; }
    int getHeight() const { return height; }
    double getFOV() const { return fov; }
};

class Ray {
    Vector origin, direction;

public:
    Ray(const Vector& _origin, const Vector& _direction) : origin(_origin), direction(_direction) {}

    Vector getRayOrigin() const { return origin; }
    Vector getRayDirection() const { return direction; }
};

class RayHit {
    Ray ray;
    Vector hitPosition;
    Vector hitNormal;
    Colour meshColour;
    bool didHit;

public:
    RayHit(Ray _ray = Ray(Vector(0,0,0), Vector(0,0,0)), Vector _hitPosition = Vector(0,0,0), Vector _hitNormal = Vector(0,0,0), Colour _meshColour = Colour(0,0,0,0), bool _didHit = false)
    : ray(_ray), hitPosition(_hitPosition), hitNormal(_hitNormal), meshColour(_meshColour), didHit(_didHit) {}

    Ray getRay() {return ray;}
    Vector getHitPosition() {return hitPosition;}
    Vector getHitNormal() {return hitNormal;}
    Colour getMeshColour() {return meshColour;}
    bool DidHit() {return didHit;}
    void setMeshColour(Colour colour) {meshColour = colour;}
};


struct RGBType{
    double r;
    double g;
    double b;
};

void saveBMP(const char* filename, int width, int height, int dpi, Colour** pixels) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        // Handle file opening error
        return;
    }

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

class Light {
    Vector position;
    Colour colour;

public:
    Light(const Vector& _position, const Colour& _colour)
        : position(_position), colour(_colour) {}

    Vector getPosition() const { return position; }
    Colour getColour() const { return colour; }
};

class Plane {
    // Plane Transform
    Vector position;
    Vector rotation;
    Vector scale;
    Vector normal;
    Colour colour;

    // Vertices
    Vector vertices[4];


public:
    Plane(const Vector _position, const Vector _rotation, const Vector _scale) 
        : position(_position), rotation(_rotation), scale(_scale) {
        this->colour = Colour(0.8, 0.8, 0.8, 1);    
        CalculateVertices();
        CalculateNormal();
    }

    const Vector Vertex(int index) const {return vertices[index];}
    const Vector Position() const {return position;}
    const Vector Rotation() const {return rotation;}
    const Vector Scale() const {return scale;}
    const Vector Normal() const {return normal;}

    RayHit CalculateRayLighting(const Ray ray)
    {
        // Ray intersects with face 1 or 2
        RayHit hit1 = RayTriangleIntersect(ray, vertices[1], vertices[0], vertices[2]);  
        RayHit hit2 = RayTriangleIntersect(ray, vertices[2], vertices[0], vertices[3]);

        RayHit hit;
        if (hit1.DidHit()) {hit = hit1;} // ray collided with face 1
        if (hit2.DidHit()) {hit = hit2;} // ray collided with face 2

        // Calculate the colour
        if (hit.DidHit()){

        // Get the light direction from the ray hit point
        Vector lightDirection = Vector(3, 8, 20) - hit.getHitPosition();

        // Calculate the dot product of the surface normal and light direction
        double cosTheta = hit.getHitNormal().dotProduct(lightDirection);

        // Calculate the brightness based on the cosine of the angle
        double brightness = std::max(0.0, cosTheta);

        // Multiply the mesh colour by the brightness
        Colour meshColour = hit.getMeshColour();
        meshColour = meshColour.multiply(brightness);

        // Set the modified mesh colour back to the hit object
        hit.setMeshColour(meshColour);

        }

        return hit;
    }


private:
    void CalculateVertices() {
        // Calculate the vertices based on position, rotation, and scale
        // Assuming the initial vertices are defined in local space
        Vector vertex1(-0.5, 0, -0.5);
        Vector vertex2(-0.5, 0, 0.5);
        Vector vertex3(0.5, 0, 0.5);
        Vector vertex4(0.5, 0, -0.5);

        // Apply scale
        vertex1 = Vector(vertex1.getX() * scale.getX(), vertex1.getY() * scale.getY(), vertex1.getZ() * scale.getZ());
        vertex2 = Vector(vertex2.getX() * scale.getX(), vertex2.getY() * scale.getY(), vertex2.getZ() * scale.getZ());
        vertex3 = Vector(vertex3.getX() * scale.getX(), vertex3.getY() * scale.getY(), vertex3.getZ() * scale.getZ());
        vertex4 = Vector(vertex4.getX() * scale.getX(), vertex4.getY() * scale.getY(), vertex4.getZ() * scale.getZ());

        // Apply rotation
        vertex1 = vertex1.rotate(rotation);
        vertex2 = vertex2.rotate(rotation);
        vertex3 = vertex3.rotate(rotation);
        vertex4 = vertex4.rotate(rotation);

        // Apply position
        vertex1 = vertex1 + position;
        vertex2 = vertex2 + position;
        vertex3 = vertex3 + position;
        vertex4 = vertex4 + position;

        vertices[0] = vertex1;
        vertices[1] = vertex2;
        vertices[2] = vertex3;
        vertices[3] = vertex4;
    }

    void CalculateNormal() {
        // Calculate the normal direction based on the vertices
        Vector edge1 = vertices[1] - vertices[0];
        Vector edge2 = vertices[2] - vertices[0];

        normal = edge1.crossProduct(edge2).normalized();
    }

    RayHit RayTriangleIntersect(const Ray& ray, const Vector& vertex1, const Vector& vertex2, const Vector& vertex3)
    {
        Vector edge1 = vertex2 - vertex1;
        Vector edge2 = vertex3 - vertex1;

        Vector p = ray.getRayDirection().crossProduct(edge2);
        double determinant = edge1.dotProduct(p);

        // Check if the ray is parallel to the triangle (no intersection)
        if (determinant == 0.0)
            return RayHit(ray);

        double invDeterminant = 1.0 / determinant;

        Vector t = ray.getRayOrigin() - vertex1;
        double u = t.dotProduct(p) * invDeterminant;

        // Check if the intersection point is outside the triangle
        if (u < 0.0 || u > 1.0)
            return RayHit(ray);

        Vector q = t.crossProduct(edge1);
        double v = ray.getRayDirection().dotProduct(q) * invDeterminant;

        // Check if the intersection point is outside the triangle
        if (v < 0.0 || u + v > 1.0)
            return RayHit(ray);

        double tValue = edge2.dotProduct(q) * invDeterminant;

        // Check if the intersection point is behind the ray's origin
        if (tValue < 0.0)
            return RayHit(ray);

        // Calculate the intersection position
        Vector intersectionPos = ray.getRayOrigin() + ray.getRayDirection() * tValue;

        // Calculate the face normal
        Vector faceNormal = edge1.crossProduct(edge2).normalized();

        // The ray intersects the triangle
        return RayHit(ray, intersectionPos, faceNormal, colour, true);
    }
};


void Render(Camera camera, Plane plane) {
    std::cout << "Rendering Scene..." << std::endl;

    // Constant for PI
    double Pi = 3.14159265358979323846;

    // Create a 2D array of booleans
    bool** grid = new bool*[camera.getHeight()];
    for (int i = 0; i < camera.getHeight(); i++) {
        grid[i] = new bool[camera.getWidth()];
    }

    // Create a 2D array of colours
    Colour** colourArray = new Colour*[camera.getHeight()];
    for (int i = 0; i < camera.getHeight(); i++) {
        colourArray[i] = new Colour[camera.getWidth()];
    }


    // run code for every pixel in camera resolution
    for (int y = 0; y < camera.getHeight(); y++) {
        for (int x = 0; x < camera.getWidth(); x++) {
            double ndcX = (2.0 * (x + 0.5) / camera.getWidth()) - 1.0;
            double ndcY = 1.0 - (2.0 * (y + 0.5) / camera.getHeight());

            double aspectRatio = camera.getWidth() / static_cast<double>(camera.getHeight());

            double fovY = camera.getFOV() * (Pi / 180.0);
            double fovX = 2.0 * std::atan(std::tan(fovY / 2.0) * aspectRatio);
            double angleX = ndcX * fovX / 2.0;
            double angleY = ndcY * fovY / 2.0;

            // Calculate the ray direction for the current pixel
            Vector rayDirection = (camera.Forward() + camera.Right() * std::tan(angleX) + camera.Down() * std::tan(angleY)).normalized();

            // Cast the ray with the calculated direction
            Ray ray(camera.Position(), rayDirection);

            RayHit hit = plane.CalculateRayLighting(ray);
            colourArray[y][x] = hit.getMeshColour();
            //grid[y][x] = hitPlane;
        }
    }

    saveBMP("first_render.bmp", camera.getWidth(), camera.getHeight(), 72, colourArray);

    // // Render the grid
    // for (int i = 0; i < camera.getHeight(); i++) {
    //     for (int j = 0; j < camera.getWidth(); j++) {
    //         if (grid[i][j]) {
    //             std::cout << '#';  // Print # for true value
    //         } else {
    //             std::cout << '.';  // Print '.' for false value
    //         }
    //     }
    //     std::cout << std::endl;  // Move to the next line
    // }
}


 

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
    int height = 400;
    double fov = 70;

    Camera cam(pos, dir, right, down, width, height, fov);

    // Create a plane
    Plane plane(Vector(0, 0, 1), Vector(-90, 0, 45), Vector(1, 1, 1));

    // Create a light
    Light light(Vector(2, 6, 4), Colour(1,1,1,1));


    Render(cam, plane);

    return 0;
}
