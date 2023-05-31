#ifndef OBJECT_H
#define OBJECT_H

#include "Vector.h"

class Object {
protected:
    Vector position;
    Vector rotation;
    Vector scale;

public:
    Object(const Vector& _position = Vector(), const Vector& _rotation = Vector(), const Vector& _scale = Vector(1.0, 1.0, 1.0))
        : position(_position), rotation(_rotation), scale(_scale) {}

    // Getters and setters for position, rotation, and scale
    const Vector& getPosition() const { return position; }
    const Vector& getRotation() const { return rotation; }
    const Vector& getScale() const { return scale; }

    void setPosition(const Vector& _position) { position = _position; }
    void setRotation(const Vector& _rotation) { rotation = _rotation; }
    void setScale(const Vector& _scale) { scale = _scale; }
};

#endif  
