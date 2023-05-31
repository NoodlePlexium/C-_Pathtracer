#ifndef PLANE_H
#define PLANE_H

#include "Vector.h"
#include "RayHit.h"

class Plane : public Object {
    Vector normal;
    Colour colour;
    Vector vertices[4];

public:
    Plane(const Vector& position, const Vector& rotation, const Vector& scale)
        : Object(position, rotation, scale), colour(0.8, 0.4, 0.2, 1) {
        CalculateVertices();
        CalculateNormal();
    }

    const Vector& Vertex(int index) const { return vertices[index]; }
    const Vector& Normal() const { return normal; }

    RayHit CalculateRayHit(const Ray& ray) {
        // Ray intersects with face 1 or 2
        RayHit hit1 = RayTriangleIntersect(ray, vertices[1], vertices[0], vertices[2]);
        RayHit hit2 = RayTriangleIntersect(ray, vertices[2], vertices[0], vertices[3]);

        RayHit hit;
        if (hit1.DidHit()) hit = hit1; // ray collided with face 1
        if (hit2.DidHit()) hit = hit2; // ray collided with face 2

        return hit;
    }

    void setColour(Colour _colour){
        colour = _colour;
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
        Vector edge1 = vertices[1] - vertices[0];
        Vector edge2 = vertices[2] - vertices[0];

        normal = edge1.crossProduct(edge2).normalized();
    }

    RayHit RayTriangleIntersect(const Ray& ray, const Vector& vertex1, const Vector& vertex2, const Vector& vertex3) {
        Vector edge1 = vertex2 - vertex1;
        Vector edge2 = vertex3 - vertex1;

        Vector p = ray.getRayDirection().crossProduct(edge2);
        double determinant = edge1.dotProduct(p);

        if (determinant == 0.0)
            return RayHit(ray);

        double invDeterminant = 1.0 / determinant;

        Vector origin = ray.getRayOrigin() - vertex1;
        double u = origin.dotProduct(p) * invDeterminant;

        if (u < 0.0 || u > 1.0)
            return RayHit(ray);

        Vector q = origin.crossProduct(edge1);
        double v = ray.getRayDirection().dotProduct(q) * invDeterminant;

        if (v < 0.0 || u + v > 1.0)
            return RayHit(ray);

        double dist = edge2.dotProduct(q) * invDeterminant;

        if (dist < 0.0)
            return RayHit(ray);

        Vector intersectionPos = ray.getRayOrigin() + ray.getRayDirection() * dist;
        Vector faceNormal = edge1.crossProduct(edge2).normalized();

        return RayHit(ray, intersectionPos, dist, faceNormal, colour, true);
    }
};

#endif