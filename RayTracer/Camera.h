#ifndef CAMERA_H
#define CAMERA_H

#include "Vector.h"

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

#endif
