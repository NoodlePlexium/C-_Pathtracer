#ifndef COLOUR_H
#define COLOUR_H


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
    Colour operator*(double scalar) const {
        double newRed = red * scalar;
        double newGreen = green * scalar;
        double newBlue = blue * scalar;
        double newAlpha = alpha * scalar;
        return Colour(newRed, newGreen, newBlue, newAlpha);
    }
};

#endif