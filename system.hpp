#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include "spring.hpp"
#include "mass.hpp"
#include <math.h>

using namespace std;

class MassSpringSystem
{
public:
	vector<Mass*> masses;
    vector<Spring*> springs;
	MassSpringSystem()
	{
	}

	Mass* addMass(float mass, float x, float y, float z)
	{
		Mass *m = new Mass(mass, x, y, z);
	    masses.push_back(m);
	    return m;
	}

    Spring* addSpring(float springConstant, float restLength, Mass *mass1, Mass *mass2)
    {
    	Spring *s = new Spring(springConstant, restLength, mass1, mass2);
	    mass1->addSpring(s);
	    mass2->addSpring(s);
	    springs.push_back(s);
	    return s;
    }

    void update(float dt,int choice, vec3 gravity, vec3 wind, float AT, float AN)
    {
    	
    	if(gravity.norm()==0) // simple spring
    	{
    		vec3 windForce;
    		for (int i = 0; i < masses.size(); ++i)
	    	{
	    		if (masses[i]->mass==0) 
		        		continue;

	        	if(choice==0) // Symplectic Euler
			    {
			        vec3 acc = masses[i]->calculateForces(gravity,windForce) / masses[i]->mass;
			        masses[i]->velocity = masses[i]->velocity + acc * dt;
			        masses[i]->position = masses[i]->position + masses[i]->velocity * dt;
			    }
			    else if (choice==1) // Midpoint
			    {
			        vec3 init_pos = masses[i]->position;
			        vec3 init_vel = masses[i]->velocity;
			        masses[i]->velocity = masses[i]->velocity + masses[i]->calculateForces(gravity,windForce) / masses[i]->mass * dt/2;
			        masses[i]->position = masses[i]->position + masses[i]->velocity * dt/2;        
			        vec3 v_half = masses[i]->velocity;
			        masses[i]->velocity = init_vel + masses[i]->calculateForces(gravity,windForce) / masses[i]->mass * dt;
			        masses[i]->position = init_pos + v_half * dt;
			    }
			    else if (choice==2) // RK4
			    {
			        vec3 k1,k2,k3,k4,l1,l2,l3,l4;
			        vec3 init_pos = masses[i]->position;
			        vec3 init_vel = masses[i]->velocity;

			        l1 = masses[i]->calculateForces(gravity,windForce) / masses[i]->mass;
			        k1 = masses[i]->velocity;

			        masses[i]->position = init_pos + k1 * dt/2;
			        masses[i]->velocity = init_vel + l1 * dt/2;
			        l2 = masses[i]->calculateForces(gravity,windForce) / masses[i]->mass;
			        k2 = masses[i]->velocity;

			        masses[i]->position = init_pos + k2 * dt/2;
			        masses[i]->velocity = init_vel + l2 * dt/2;
			        l3 = masses[i]->calculateForces(gravity,windForce) / masses[i]->mass;
			        k3 = masses[i]->velocity;

			        masses[i]->position = init_pos + k3 * dt;
			        masses[i]->velocity = init_vel + l3 * dt;
			        l4 = masses[i]->calculateForces(gravity,windForce) / masses[i]->mass;
			        k4 = masses[i]->velocity;

			        masses[i]->velocity = init_vel + (1/6.0*(l1)+1/3.0*(l2)+1/3.0*(l3)+1/6.0*(l4)) * dt;
			        masses[i]->position = init_pos + (1/6.0*(k1)+1/3.0*(k2)+1/3.0*(k3)+1/6.0*(k4)) * dt;   
			    }
			    else
			        exit(0);
	    	}
    	}
		else
		{
	    	int side = (int)sqrt(masses.size());
	    	for (int i = 0; i < side; ++i)
	    	{
		    	for(int j = 0;j < side; ++j)
		    	{
		    		vec3 windForce;    	
		    		if (masses[i*side+j]->mass==0) 
		        		continue;

				    if (wind.norm()!=0)
			    	{
			    		vec3 normal;
			    		float area;
			    		if(i>0 && j>0)
			    		{
			    			normal = (masses[(i-1)*side+j]->position-masses[i*side+j]->position).cross(masses[i*side+j-1]->position-masses[i*side+j]->position);
			    			area = normal.norm();
			    			normal = normal/normal.norm();
			    		}
			    		vec3 vn = normal.dot(wind)*normal;
			    		vec3 vt = wind - vn;
			    		vec3 fn = AN*area*wind.dot(vn)*normal;
			    		vec3 ft = AT*area*vt;
			    		windForce = ft + fn;
			    	}

				    if(choice==0) // Symplectic Euler
				    {
				        vec3 acc = masses[i*side+j]->calculateForces(gravity,windForce) / masses[i*side+j]->mass;
				        masses[i*side+j]->velocity = masses[i*side+j]->velocity + acc * dt;
				        masses[i*side+j]->position = masses[i*side+j]->position + masses[i*side+j]->velocity * dt;
				    }
				    else if (choice==1) // Midpoint
				    {
				        vec3 init_pos = masses[i*side+j]->position;
				        vec3 init_vel = masses[i*side+j]->velocity;
				        masses[i*side+j]->velocity = masses[i*side+j]->velocity + masses[i*side+j]->calculateForces(gravity,windForce) / masses[i*side+j]->mass * dt/2;
				        masses[i*side+j]->position = masses[i*side+j]->position + masses[i*side+j]->velocity * dt/2;        
				        vec3 v_half = masses[i*side+j]->velocity;
				        masses[i*side+j]->velocity = init_vel + masses[i*side+j]->calculateForces(gravity,windForce) / masses[i*side+j]->mass * dt;
				        masses[i*side+j]->position = init_pos + v_half * dt;
				    }
				    else if (choice==2) // RK4
				    {
				        vec3 k1,k2,k3,k4,l1,l2,l3,l4;
				        vec3 init_pos = masses[i*side+j]->position;
				        vec3 init_vel = masses[i*side+j]->velocity;

				        l1 = masses[i*side+j]->calculateForces(gravity,windForce) / masses[i*side+j]->mass;
				        k1 = masses[i*side+j]->velocity;

				        masses[i*side+j]->position = init_pos + k1 * dt/2;
				        masses[i*side+j]->velocity = init_vel + l1 * dt/2;
				        l2 = masses[i*side+j]->calculateForces(gravity,windForce) / masses[i*side+j]->mass;
				        k2 = masses[i*side+j]->velocity;

				        masses[i*side+j]->position = init_pos + k2 * dt/2;
				        masses[i*side+j]->velocity = init_vel + l2 * dt/2;
				        l3 = masses[i*side+j]->calculateForces(gravity,windForce) / masses[i*side+j]->mass;
				        k3 = masses[i*side+j]->velocity;

				        masses[i*side+j]->position = init_pos + k3 * dt;
				        masses[i*side+j]->velocity = init_vel + l3 * dt;
				        l4 = masses[i*side+j]->calculateForces(gravity,windForce) / masses[i*side+j]->mass;
				        k4 = masses[i*side+j]->velocity;

				        masses[i*side+j]->velocity = init_vel + (1/6.0*(l1)+1/3.0*(l2)+1/3.0*(l3)+1/6.0*(l4)) * dt;
				        masses[i*side+j]->position = init_pos + (1/6.0*(k1)+1/3.0*(k2)+1/3.0*(k3)+1/6.0*(k4)) * dt;   
				    }
				    else
				        exit(0);
		    	}
		    }
		}
	}
};

#endif
