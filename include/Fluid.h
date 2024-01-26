#ifndef FLUID_H_
#define FLUID_H_
#include <ngl/Vec3.h>
#include <cstdlib>

enum class ParticleState : bool {Alive,Dead};
struct Particle
{
    Particle()=default;
    ngl::Vec3 pos;
    ngl::Vec3 dir;
    ngl::Vec3 colour = {0.0f, 0.35f, 0.85f};
    ngl::Vec3 velocity;
    ngl::Vec3 acceleration;

    float restDensity = .01f;
//    float size=1.0f;
    float density = 0.5f;
    float mass = static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX);
    float stiffness = 1.0f;
    int life=100;
    ParticleState alive = ParticleState::Dead;
};

struct Fluid
{
    std::vector<Particle> particles;
};

#endif