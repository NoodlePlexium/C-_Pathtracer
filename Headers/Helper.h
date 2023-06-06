#ifndef HELPER_H
#define HELPER_H

#include "Vector.h"
#include <random>

class Helper{

public:

	float Random(float min, float max) {
        return min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max - min)));
    }

    float RandomValueNormalDistribution(){
    	float theta = 2 * 3.1415926 * Random(0, 1);
    	float rho = std::sqrt(-2 * std::log(Random(0, 1)));
    	return rho * std::cos(theta);
    }

    Vector RandomDirection(){
    	float x = RandomValueNormalDistribution();
        float y = RandomValueNormalDistribution();
        float z = RandomValueNormalDistribution();
        Vector randomDir(x,y,z);
        return randomDir.normalized();
    }

   	Vector RandomHemisphereDirection(const Vector& normal) {

        float x = RandomValueNormalDistribution();
        float y = RandomValueNormalDistribution();
        float z = RandomValueNormalDistribution();

        Vector randomDirection = RandomDirection();

        if (randomDirection.dotProduct(normal) < 0.0f){
        	randomDirection = randomDirection * -1.0f;
        }

        return randomDirection;
    }

};

#endif
