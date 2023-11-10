#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include "Fluid.h"

class Simulator
{
public:
    Simulator(Vec3 _pos);
    void draw() const;
    void update();
private:
    void rebound(Fluid &_f); //to handle collision with wall and surfaces etc

};

#endif
