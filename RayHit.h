#ifndef RAYHIT_H
#define RAYHIT_H

#include "Vector.h"
#include "Colour.h"
#include "Ray.h"

class RayHit {
    Ray ray;
    Vector hitPosition;
    Vector hitNormal;
    double hitDistance;
    Colour meshColour;
    bool didHit;

public:
    RayHit(Ray _ray = Ray(Vector(0,0,0), Vector(0,0,0)), Vector _hitPosition = Vector(0,0,0), double _hitDistance = std::numeric_limits<double>::infinity(), Vector _hitNormal = Vector(0,0,0), Colour _meshColour = Colour(0,0,0,0), bool _didHit = false)
    : ray(_ray), hitPosition(_hitPosition), hitDistance(_hitDistance), hitNormal(_hitNormal), meshColour(_meshColour), didHit(_didHit) {}

    Ray getRay() {return ray;}
    Vector getHitPosition() {return hitPosition;}
    Vector getHitNormal() {return hitNormal;}
    Colour getMeshColour() {return meshColour;}
    double getHitDistance() {return hitDistance;}
    bool DidHit() {return didHit;}
    void setMeshColour(Colour colour) {meshColour = colour;}
};

#endif