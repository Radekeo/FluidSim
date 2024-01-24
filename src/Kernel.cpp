//
// Created by s5636764 on 23/01/24.
//
#include "Kernel.h"
#include <cmath>

float Kernel::smoothingKernel(float _rSquared)
{
    // @smoothingKernel Equation
    // Author: Andrew Gibiansky
    // https://andrew.gibiansky.com/blog/physics/computational-fluid-dynamics/
    // Wg(r,h)=(315/(64π(h^9))) * ((h^2)−(r^2))^3.

    float coefficient = 315.0f/ (64.0f * M_PI * std::pow(m_h, 9));
    float kernel = coefficient * std::pow((std::pow(m_h, 2) - _rSquared), 3);

    return kernel;
}

float Kernel::smoothingKernelLaplacian(Particle &_p1, Particle &_p2)
{
    // @Viscocity Kernel Equation
    // Author: Andrew Gibiansky
    // https://andrew.gibiansky.com/blog/physics/computational-fluid-dynamics/
    // Wg_laplacian(r,h) = −((r^3)/2h^3)+(r^2/h^2)+(h/2r)−1.
    ngl::Vec3 dist = _p1.pos - _p2.pos;
    float rSquared = dist.lengthSquared();
    float r = std::sqrt(rSquared);

    float coefficient1 = 1.0f/(2.0f * std::pow(m_h, 3));
    float coefficient2 = 1.0f/std::pow(m_h, 2);
    float coefficient3 = m_h;
    return (-(coefficient1 * (std::pow(r, 3))) + (coefficient2 * std::pow(r,2)) + (coefficient3/(2*r))) - 1;
}

ngl::Vec3 Kernel::pressureGrad(Particle &_p1, Particle &_p2)
{
    // @Pressure Kernel Equation
    // Author: Andrew Gibiansky
    // https://andrew.gibiansky.com/blog/physics/computational-fluid-dynamics/
    // Wp(r,h)=15πh6(h−r)3
    ngl::Vec3 dist = _p1.pos - _p2.pos;
    float rSquared = dist.lengthSquared();
    float r = std::sqrt(rSquared);

    float coefficient = 15.0f /(M_PI * std::pow(m_h, 6));
    float k = coefficient * std::pow((m_h -r),3);
    ngl::Vec3 kernel = {k,k,k};

    return kernel;
}

ngl::Vec3 Kernel::smoothingKernelGrad(const ngl::Vec3 _r)
{
    // Gradient of the smoothing kernel
    // Wg_grad(r,h) = d(Wg(r,h))/dr
    float coefficient = -945.0f / (32.0f * M_PI * std::pow(m_h, 9));
    float q = std::sqrt(_r.lengthSquared()) / m_h;
    return coefficient * _r * std::pow((1.0f - q * q), 2);
}

//const bool Kernel::nonZeroDist(Particle &_p1, Particle &_p2)
//{
////    float normsq = pow((_p1.pos.m_x- _p2.pos.m_x ), 2)+ pow((_p1.pos.m_y- _p2.pos.m_y ), 2) + pow((_p1.pos.m_z- _p2.pos.m_z), 2);
//    float normSq =  (_p1.pos - _p2.pos).lengthSquared();
//    float sq = std::pow(m_h,2);
//    if(normSq > sq)
//    {
//        return false;
//    }
//    return true;
//}
