#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include "Fluid.h"

class FluidSimulator
{
public:
    FluidSimulator(size_t _numParticles);
    size_t numParticles() const;
    void Initialize();
    void Draw() const;
    void Update();
private:

};

class SPHSimulator
{
    void Initialize(); //smoothed particle hydrodynamics
};

#endif
