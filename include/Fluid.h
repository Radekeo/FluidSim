#ifndef FLUID_H
#define FLUID_H

#include "Vec3.h"

struct FluidParticle
{
    FluidParticle() = default;
    Vec3 pos;
    Vec3 dir;
    Vec3 colour;
    float size = 1.0f;
    int life = 100;
};
#endif
