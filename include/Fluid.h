#ifndef FLUID_H
#define FLUID_H

#include "Vec3.h"

struct Fluid
{
    Fluid() = default;
    Vec3 pos;
    Vec3 colour;
    Vec3 dir;

    float size;
    float viscosity;
    //float velocity;
    float density;
};
#endif
