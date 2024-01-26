//
// Created by s5636764 on 23/01/24.
//

#ifndef KERNEL_H_
#define KERNEL_H_

#include "Fluid.h"

struct Kernel
{
public :
    Kernel() = default;
    float smoothingKernel(float _rSquared);
    ngl::Vec3 smoothingKernelGrad(ngl::Vec3 _r);
    float smoothingKernelLaplacian(Particle &_p1, Particle &_p2);
    ngl::Vec3 pressureGrad(Particle &_p1, Particle &_p2);

private:
    float m_h = 1000.0f; //smoothing bandwidth


};
#endif
