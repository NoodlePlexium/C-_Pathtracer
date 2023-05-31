#ifndef SPHERE_H
#define SPHERE_H

#include "Vector.h"

class Sphere {
    Vector center;
    double radius;
    Colour color;

public:
    Sphere(const Vector& _center, double _radius, const Colour& _color)
        : center(_center), radius(_radius), color(_color) {}

    RayHit CalculateRayHit(const Ray& ray, Vector lightPos) {
        RayHit hit;

        Vector oc = ray.getRayOrigin() - center;
        double a = ray.getRayDirection().dotProduct(ray.getRayDirection());
        double b = 2.0 * oc.dotProduct(ray.getRayDirection());
        double c = oc.dotProduct(oc) - radius * radius;
        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0) {
            return hit;
        }

        double sqrtDiscriminant = std::sqrt(discriminant);
        double denominator = 2.0 * a;

        // Find the closest intersection distance
        double t1 = (-b - sqrtDiscriminant) / denominator;
        double t2 = (-b + sqrtDiscriminant) / denominator;

        // Ray hits
        if (t1 > 0 || t2 > 0) {
            double t = (t1 > 0) ? t1 : t2;

            // Calculate the intersection position
            Vector hitPos = ray.getRayOrigin() + ray.getRayDirection() * t;
            Vector normal = (center - hitPos).normalized();

            return RayHit(ray, hitPos, t, normal, color, true);
        }
        return hit;
    }

};

#endif