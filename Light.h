#ifndef LIGHT_H
#define LIGHT_H

#include "Object.h"
#include "Colour.h"
#include "Vector.h"

class Light : public Object {
    Colour colour;
    Vector direction;

public:
    Light(const Vector& _position, const Vector& _direction, const Colour& _colour)
        : Object(_position), direction(_direction.normalized()), colour(_colour) {}

    Colour getColour() const { return colour; }
    Vector getDirection() const {return direction; }
};

#endif
