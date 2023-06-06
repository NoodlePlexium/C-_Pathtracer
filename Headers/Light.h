#ifndef LIGHT_H
#define LIGHT_H

#include "Colour.h"
#include "Vector.h"

class Light{
    Colour colour;
    Vector direction;
    Vector position;
    double intensity = 1.6;

public:
    Light(const Vector& _position, const Vector& _direction, const Colour& _colour)
        : position(_position), direction(_direction.normalized()), colour(_colour) {}

    Colour getColour() const { return colour; }
    Vector getDirection() const {return direction; }
    const Vector& getPosition() const { return position; }
    double getIntensity() const { return intensity; }
};

#endif
