#ifndef FLUID_H_
#define FLUID_H_
#include <ngl/Vec3.h>
#include <cstdlib>

struct Particle
{
    Particle()=default;
    ngl::Vec3 pos;
    ngl::Vec3 dir;
    ngl::Vec3 colour;
    ngl::Vec3 velocity;
    ngl::Vec3 acceleration;

    float restDensity = .01f;
    float size=0.5f;
    float density = 1.0f;
    float pressure = 0.5f;
    float mass = static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX);
    float stiffness = 0.2f;
    int life=200;
};

struct Fluid
{
    std::vector<Particle> particles;
};

#endif