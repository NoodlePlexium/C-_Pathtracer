#ifndef SPHERE_H
#define SPHERE_H

#include "Vector.h"
#include "Material.h"

class Sphere : public Object {
    double radius;
    Material material;

public:
    Sphere(const Vector& position, const Vector& rotation, const Vector& scale, double _radius, const Material& _material)
        : Object(position, rotation, scale), radius(_radius), material(_material) {}

    RayHit CalculateRayHit(const Ray& ray) const override {
        RayHit hit;

        Vector oc = ray.getRayOrigin() - position;
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
            Vector normal = (position - hitPos).normalized();
            return RayHit(ray, hitPos, t, normal, material, true);
        }
        return hit;
    }

};

#endif
