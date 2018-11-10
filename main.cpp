#include "camera.hpp"
#include "draw.hpp"
#include "gui.hpp"
#include "lighting.hpp"
#include "text.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <cmath>
#include <unordered_map>

using namespace std;

Window window;
Camera camera;
Lighting lighting;
Text text;

float dt = 0.001;
float t = 0;
bool paused = false;

/////////////////////////////////////////////////

float rand01()
{
    return (float)rand() * (1.f / RAND_MAX);
}

struct Particle;
struct Neighbor
{
    Particle* j;
    float q, q2;
};

struct Particle
{
    float r;
    float g;
    float b;

    vec2 pos;
    vec2 pos_old;
    vec2 vel;
    vec2 force;
    float mass;
    float rho;
    float rho_near;
    float press;
    float press_near;
    float sigma;
    float beta;
    std::vector< Neighbor > neighbors;
};

std::vector< Particle > particles;



//PARAMETERS    

const float G = .02f * .25;           // Gravitational Constant for our simulation
const float spacing = 1.f;            // Spacing of particles
const float k = spacing / 1000.0f;    // Far pressure weight
const float k_near = k * 10;          // Near pressure weight
const float rest_density = 3;         // Rest Density
const float r = spacing * 1.25f;      // Radius of Support
const float rsq = r * r;              // ... squared for performance stuff
const float SIM_W = 50;               // The size of the world
const float bottom = 0;               // The floor of the world

void init( const unsigned int N )
{
    float w = SIM_W / 4;
    for( float y = bottom + 1; y <= 10000; y += r * 0.5f )
    {
        for(float x = -w; x <= w; x += r * 0.5f )
        {
            if( particles.size() > N )
            {
                break;
            }

            Particle p;
            p.pos = vec2(x, y);
            p.pos_old = p.pos;// + 0.001f * glm::vec2(rand01(), rand01());
            p.force = vec2(0,0);
            p.sigma = 3.f;
            p.beta = 4.f;
            particles.push_back(p);
        }
    }
}

template< typename T >
class SpatialIndex
{
public:
    typedef std::vector< T* > NeighborList;

    SpatialIndex
        (
        const unsigned int numBuckets,  // number of hash buckets
        const float cellSize,           // grid cell size
        const bool twoDeeNeighborhood   // true == 3x3 neighborhood, false == 3x3x3
        )
        : mHashMap( numBuckets )
        , mInvCellSize( 1.0f / cellSize )
    {
        // initialize neighbor offsets
        for( int i = -1; i <= 1; i++ )
            for( int j = -1; j <= 1; j++ )
                if( twoDeeNeighborhood )
                    mOffsets.push_back( ivec3( i, j, 0 ) );
                else
                    for( int k = -1; k <= 1; k++ )
                        mOffsets.push_back( ivec3( i, j, k ) );
    }

    void Insert( const vec3& pos, T* thing )
    {
        mHashMap[ Discretize( pos, mInvCellSize ) ].push_back( thing );
    }

    void Neighbors( const vec3& pos, NeighborList& ret ) const
    {
        const ivec3 ipos = Discretize( pos, mInvCellSize );
        for( const auto& offset : mOffsets )
        {
            typename HashMap::const_iterator it = mHashMap.find( offset + ipos );
            if( it != mHashMap.end() )
            {
                ret.insert( ret.end(), it->second.begin(), it->second.end() );
            }
        }
    }

    void Clear()
    {
        mHashMap.clear();
    }

private:
    // "Optimized Spatial Hashing for Collision Detection of Deformable Objects"
    // Teschner, Heidelberger, et al.
    // returns a hash between 0 and 2^32-1
    struct TeschnerHash : std::unary_function< ivec3, std::size_t >
    {
        std::size_t operator()( ivec3 const& pos ) const
        {
            const unsigned int p1 = 73856093;
            const unsigned int p2 = 19349663;
            const unsigned int p3 = 83492791;
            return size_t( ( pos[0] * p1 ) ^ ( pos[1] * p2 ) ^ ( pos[2] * p3 ) );
        };
    };

    // returns the indexes of the cell pos is in, assuming a cellSize grid
    // invCellSize is the inverse of the desired cell size
    static inline ivec3 Discretize( const vec3& pos, const float invCellSize )
    {
        return ivec3( floor( pos[0] * invCellSize ),floor( pos[1] * invCellSize ),floor( pos[2] * invCellSize ) );
    }

    typedef std::unordered_map< ivec3, NeighborList, TeschnerHash > HashMap;
    HashMap mHashMap;

    std::vector< ivec3 > mOffsets;

    const float mInvCellSize;
};

typedef SpatialIndex< Particle > IndexType;
IndexType indexsp( 4093, r, false );

void drawStuff() 
{
    
    // if(OPTION==0)
    // {
    //     setPointSize(10);
    //     drawPoint(mss.masses.at(0)->position);
    //     drawPoint(mss.masses.at(1)->position);
    //     drawLine(mss.masses.at(0)->position,mss.masses.at(1)->position);
    // }
    // else if (OPTION==1)
    // {
    //     setPointSize(10);
    //     drawPoint(p0);
    //     drawPoint(p1);
    //     setColor(vec3(0.8,0.2,0.2));
    //     for (int i=1; i<SIZE; i++) 
    //     {
    //         for (int j = 1; j < SIZE; ++j)
    //         {
    //             drawQuad(mss.masses.at(SIZE*i+j)->position,mss.masses.at(SIZE*i+j-1)->position,mss.masses.at(SIZE*(i-1)+j-1)->position,mss.masses.at(SIZE*(i-1)+j)->position);
    //         }
    //     }
    // }
    // else
    // {
    //     setColor(vec3(0.0,0.0,0.0));
    //     setLineWidth(4);
    //     drawLine(vec3(0.0,0.0,0.0),vec3(0,-SIZE*STRUCT_LEN*4,0));
    //     setColor(vec3(0.2,0.2,0.8));
    //     for (int i=1; i<SIZE; i++) 
    //     {
    //         for (int j = 1; j < SIZE; ++j)
    //         {
    //             drawQuad(mss.masses.at(SIZE*i+j)->position,mss.masses.at(SIZE*i+j-1)->position,mss.masses.at(SIZE*(i-1)+j-1)->position,mss.masses.at(SIZE*(i-1)+j)->position);
    //         }
    //     }
    // }
}

void drawWorld() 
{
    camera.apply(window);
    lighting.apply();
    clear(vec3(0.9,0.9,0.9));
    drawStuff();
    setColor(vec3(0,0,0));
    text.draw("WASD and LShift/LCtrl to move camera", -0.9, 0.90);
    text.draw("Mouse to rotate view", -0.9, 0.85);
    text.draw("Space to play/pause animation", -0.9, 0.80);
}

void update(float dt) 
{
    t += dt;
    for( int i = 0; i < (int)particles.size(); ++i )
    {
        // Apply the currently accumulated forces
        particles[i].pos += particles[i].force;

        // Restart the forces with gravity only. We'll add the rest later.
        particles[i].force = vec2( 0.0f, -::G );

        // Calculate the velocity for later.
        particles[i].vel = particles[i].pos - particles[i].pos_old;

        // If the velocity is really high, we're going to cheat and cap it.
        // This will not damp all motion. It's not physically-based at all. Just
        // a little bit of a hack.
        const float max_vel = 2.0f;
        const float vel_mag = particles[i].vel.dot(particles[i].vel);
        // If the velocity is greater than the max velocity, then cut it in half.
        if( vel_mag > max_vel * max_vel )
        {
            particles[i].vel *= .5f;
        }

        // Normal verlet stuff
        particles[i].pos_old = particles[i].pos;
        particles[i].pos += particles[i].vel;

        // If the Particle is outside the bounds of the world, then
        // Make a little spring force to push it back in.
        if( particles[i].pos[0] < -SIM_W ) particles[i].force[0] -= ( particles[i].pos[0] - -SIM_W ) / 8;
        if( particles[i].pos[0] >  SIM_W ) particles[i].force[0] -= ( particles[i].pos[0] - SIM_W ) / 8;
        if( particles[i].pos[1] < bottom ) particles[i].force[1] -= ( particles[i].pos[1] - bottom ) / 8;
        // if( particles[i].pos[1] > SIM_W * 50 ) particles[i].force[1] -= ( particles[i].pos[1] - SIM_W * 50 ) / 8;

        // Reset the nessecary items.
        particles[i].rho = 0;
        particles[i].rho_near = 0;
        particles[i].neighbors.clear();
    }

    // update spatial index
    indexsp.Clear();
    for( auto& particle : particles )
    {
        indexsp.Insert( vec3( particle.pos[0],particle.pos[1], 0.0f ), &particle );
    }

    // DENSITY
    // Calculate the density by basically making a weighted sum
    // of the distances of neighboring particles within the radius of support (r)
    for( int i = 0; i < (int)particles.size(); ++i )
    {
        particles[i].rho = 0;
        particles[i].rho_near = 0;

        // We will sum up the 'near' and 'far' densities.
        float d = 0;
        float dn = 0;

        IndexType::NeighborList neigh;
        neigh.reserve( 64 );
        indexsp.Neighbors( vec3( particles[i].pos[0],particles[i].pos[1], 0.0f ), neigh );
        for( int j = 0; j < (int)neigh.size(); ++j )
        {
            if( neigh[j] == &particles[i] )
            {
                // do not calculate an interaction for a Particle with itself!
                continue;
            }

            // The vector seperating the two particles
            const vec2 rij = neigh[j]->pos - particles[i].pos;

            // Along with the squared distance between
            const float rij_len2 = rij.dot(rij);

            // If they're within the radius of support ...
            if( rij_len2 < rsq )
            {
                // Get the actual distance from the squared distance.
                float rij_len = sqrt( rij_len2 );

                // And calculated the weighted distance values
                const float q = 1 - ( rij_len / r );
                const float q2 = q * q;
                const float q3 = q2 * q;

                d += q2;
                dn += q3;

                // Set up the Neighbor list for faster access later.
                Neighbor n;
                n.j = neigh[j];
                n.q = q;
                n.q2 = q2;
                particles[i].neighbors.push_back(n);
            }
        }

        particles[i].rho += d;
        particles[i].rho_near += dn;
    }

    // PRESSURE
    // Make the simple pressure calculation from the equation of state.
    for( int i = 0; i < (int)particles.size(); ++i )
    {
        particles[i].press = k * ( particles[i].rho - rest_density );
        particles[i].press_near = k_near * particles[i].rho_near;
    }

    // PRESSURE FORCE
    // We will force particles in or out from their neighbors
    // based on their difference from the rest density.
    for( int i = 0; i < (int)particles.size(); ++i )
    {
        // For each of the neighbors
        vec2 dX = vec2( 0,0 );
        for( const Neighbor& n : particles[i].neighbors )
        {
            // The vector from Particle i to Particle j
            const vec2 rij = (*n.j).pos - particles[i].pos;

            // calculate the force from the pressures calculated above
            const float dm
                = n.q * ( particles[i].press + (*n.j).press )
                + n.q2 * ( particles[i].press_near + (*n.j).press_near );

            // Get the direction of the force
            const vec2 D = ( rij.normalized() ) * dm;
            dX += D;
        }

        particles[i].force -= dX;
    }

    // VISCOSITY
    // This simulation actually may look okay if you don't compute
    // the viscosity section. The effects of numerical damping and
    // surface tension will give a smooth appearance on their own.
    // Try it.
    for( int i = 0; i < (int)particles.size(); ++i )
    {        
        particles[i].r = 0.20f ;
        particles[i].g = 0.50f ;
        particles[i].b = 1.f ;

        // For each of that particles neighbors
        for( const Neighbor& n : particles[i].neighbors )
        {
            const vec2 rij = (*n.j).pos - particles[i].pos;
            const float l = rij.norm();
            const float q = l / r;

            const vec2 rijn = ( rij / l );
            // Get the projection of the velocities onto the vector between them.
            const float u = (particles[i].vel - (*n.j).vel).dot(rijn );
            if( u > 0 )
            {
                // Calculate the viscosity impulse between the two particles
                // based on the quadratic function of projected length.
                const vec2 I
                    = ( 1 - q )
                    * ( (*n.j).sigma * u + (*n.j).beta * u * u )
                    * rijn;

                // Apply the impulses on the current particle
                particles[i].vel -= I * 0.5f;
            }
        }
    } 
}

void keyPressed(int key) 
{
    // See http://www.glfw.org/docs/latest/group__keys.html for key codes
    if (key == GLFW_KEY_SPACE)
        paused = !paused;
    if (key == GLFW_KEY_Q)
        exit(0);
}

int main(int argc, char **argv) 
{
    window.create("Animation", 1024, 768);
    window.onKeyPress(keyPressed);
    camera.lookAt(vec3(1,1.5,5), vec3(0,0,0));
    lighting.createDefault();
    text.initialize();

    init( 2048 );   

    while (!window.shouldClose()) 
    {
        camera.processInput(window);
        if (!paused)
            update(dt);
        window.prepareDisplay();
        drawWorld();
        window.updateDisplay();
        window.waitForNextFrame(dt);
    }
}
