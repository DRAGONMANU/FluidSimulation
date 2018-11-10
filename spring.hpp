#ifndef SPRING_HPP
#define SPRING_HPP

#include "common.hpp"
#include "spring.h"

Spring::Spring(float springConstant, float restLength, Mass *mass1, Mass *mass2)
{
	k = springConstant;
	kd = 0.1;
	rl = restLength;
	m1 = mass1;
	m2 = mass2;
}

vec3 Spring::getForce(Mass *refMass)
{
	if (!(refMass == m1 or refMass == m2))
    	return vec3(0, 0, 0);

    float dist = (m1->position - m2->position).norm();
    float scalar = k*(dist - rl);
    vec3 dir = (m2->position - m1->position).normalized();

    float s1 = m1->velocity.dot(dir);
    float s2 = m2->velocity.dot(dir);
    float dampingScalar = -kd * (s1 + s2);

    if (refMass == m1)
        return (scalar + dampingScalar) * dir;
	else
        return (-scalar + dampingScalar) * dir;

}
#endif
