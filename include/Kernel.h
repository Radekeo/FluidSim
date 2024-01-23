//
// Created by s5636764 on 23/01/24.
//

#ifndef KERNEL_H_
#define KERNEL_H_

#include "Fluid.h"
class Kernel
{
public :
    Kernel();
    bool nonZeroDist(Particle _p1, Particle _p2);
    float smoothingKernel(float _r);
    ngl::Vec3 smoothingKernelGrad(const ngl::Vec3 &_r);
    float smoothingKernelLaplacian(const ngl::Vec3 &_r);
private:
    float m_h = 10.0f;
};
#endif
