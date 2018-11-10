#ifndef MASS_HPP
#define MASS_HPP

#include "common.hpp"
#include "mass.h"

Mass::Mass(float m, float x, float y, float z)
{
    mass = m;
    position = vec3(x, y, z);
    velocity = vec3(0.0, 0.0, 0.0);
}
void Mass::setPosition(float x, float y, float z)
{
    position = vec3(x,y,z);
}

void Mass::addSpring(Spring *s)
{
    springs.push_back(s);
}

vec3 Mass::calculateForces(vec3 gravity,vec3 wind)  // calculate all forces acting upon the mass
{
    vec3 fg = gravity * mass;
    vec3 fs = vec3(0.0, 0.0, 0.0);

    for (int i=0; i<springs.size(); i++) 
        fs += springs[i]->getForce(this);

    if (wind.norm()>2)
    {
        wind = wind*2/wind.norm();
    }

    vec3 force = fg + fs + wind;
    return force;
}

#endif
