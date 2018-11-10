#ifndef SPRING_H
#define SPRING_H

#include "mass.h"
#include "common.hpp"

class Mass;

class Spring
{
public:
    float k;         // spring constant
    float kd;        // damping constant
    float rl;
    Mass *m1;
    Mass *m2;

    Spring(float springConstant, float restLength, Mass *mass1, Mass *mass2);
    vec3 getForce(Mass *refMass);
};

#endif