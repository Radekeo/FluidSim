#include "Simulator.h"
#include "Kernel.h"
#include <ngl/Random.h>
//#include <iostream>
//#include <fstream>
#include <algorithm>
#include <ngl/VAOFactory.h>

#include <cmath>

Simulator::Simulator(ngl::Vec3 _pos, size_t _numParticles) : m_pos{_pos}
{
    m_fluid.particles.resize(_numParticles);
    for(auto &p : m_fluid.particles)
    {
        resetParticle(p);
    }

    birthParticles();

    m_vao = ngl::VAOFactory::createVAO(ngl::simpleVAO,GL_POINTS);

}

void Simulator::resetParticle(Particle &_p)
{
    _p.pos=m_pos;
    // Set initial direction to fall downward
   _p.dir = ngl::Vec3(0.0f, -1.0f, 0.0f);
//   _p.dir = randomVectorOnSphere(m_spread);

//  _p.dir.m_y = std::abs(_p.dir.m_y);
    _p.colour=ngl::Random::getRandomColour3();
//    _p.life = 10 + static_cast<int>(ngl::Random::randomPositiveNumber(2000.0f));
    _p.size = 0.1f;
    _p.density = 0.0f;  // Initialize density
    _p.pressure = 0.0f; // Initialize pressure

    // Randomize position
    float random_x = ngl::Random::randomNumber(0.5f);
    _p.pos.m_x += random_x;
    float y_interval = -0.1f;
    _p.pos.m_y += y_interval;
}

void Simulator::birthParticles()
{
    auto births = 0 + static_cast<int>(ngl::Random::randomPositiveNumber(m_numPerFrame));
    for(size_t i=0; i<births; ++i)
    {
        for(auto &p : m_fluid.particles)
        {
//            if(p.life <= 0)
//            {
            resetParticle(p);
            break;
//            }
        }
    }
}

void Simulator::update(float _delta)
{
    // std::cout<<"update\n";
    float dt=_delta;
    /* Find accelerations */
    for(auto &p : m_fluid.particles){
        calculateAcceleration(p);
    }

//    getPressureDensity(); // update pressure and density


    for (auto &p : m_fluid.particles)
    {
        // acceleration
        p.velocity = p.acceleration + gravity * _delta;
        // Apply gravity
//        p.dir += gravity * dt * 0.5f;
        /* Apply velocities */

    for(auto &p : m_fluid.particles){
            p.pos.m_x += p.velocity.m_x * dt;
            p.pos.m_y += p.velocity.m_x * dt;
            p.pos.m_z += p.velocity.m_x * dt;
        }

        /* Check for collisions and reverse velocity vectors if needed */
//        for(Particle p : particles){
//            applyCollisions(p);
//        }

//        // Apply pressure force
//        ngl::Vec3 pressureForce;
//        for (auto &p2 : m_fluid.particles)
//        {
//            ngl::Vec3 r = p2.pos - p.pos;
//            float rSquared = r.lengthSquared();
//
//            if (rSquared > 0.0f && rSquared < 0.2f)
//            {
//                pressureForce += (p.pressure + p2.pressure) / (2 * p2.density) * smoothingKernelGrad(r) * p2.mass;
//            }
//        }
//        p.dir += pressureForce;
//
//        // Apply viscosity force
//        ngl::Vec3 viscosityForce = applyViscosity(p);
//        p.dir += viscosityForce;
//
//        // Update particle position
//        p.pos += p.dir * 0.5f;
//        p.life -= 1;
//        p.size += 0.01f;
//        p.size = std::clamp(p.size, 0.0f, 2.0f);
//
//        // Clamp particles so they don't fall below the grid
        p.pos.m_y = std::max(p.pos.m_y, 0.0f);
    }
}
void Simulator::calculateAcceleration(Particle &_p)
{
    /* Calculate things which are needed to calculate force */
    precalculateDensities();

    /* Calculate individual force terms */
    ngl::Vec3 pressureForce = calculatePressure(_p);
    ngl::Vec3 viscosityForce = calculateViscosity(_p);
//    ngl::Vec3 externalForce = calculateExternal(_p);
//    ngl::Vec3 surfaceForce = calculateSurfaceTension(_p);

//    ngl::Vec3 totalForce = pressureForce + viscosityForce + externalForce + surfaceForce;
    ngl::Vec3 totalForce = pressureForce + viscosityForce;

    ngl::Vec3 density = ngl::Vec3(_p.density, _p.density, _p.density);
    _p.acceleration = totalForce/density;

//    float accelerationX = totalForceX / _p.density;
//    float accelerationY = totalForceY / _p.density;

//    /* Apply accelerations */
//    p.ax = accelerationX;
//    p.ay = accelerationY;

}

void Simulator::precalculateDensities()
{
    for(auto &p1 : m_fluid.particles)
    {
        float sum = 0.0f;
        for(auto &p2 : m_fluid.particles)
        {
            float product = 1.0f;
            product *= p2.mass;
            ngl::Vec3 r = p2.pos - p1.pos;
            float rSquared = r.lengthSquared();
            product *= smoothingKernel(rSquared);
            sum += product;
        }
        p1.density = sum;
    }
}

ngl::Vec3 Simulator::calculatePressure(Particle &_p)
{
    //fpressurei=−∇p(r⃗ )=−∑mjpi+pj2ρj∇W(|r⃗ i−r⃗ j|,h).
    ngl::Vec3 sum;
    float c = 100000.0f;
    float pressure1 = c * (_p.density - _p.restDensity);
    for (auto &p2: m_fluid.particles) {
        if (p2.pos != _p.pos) {
            float pressure2 = c * (p2.density - _p.restDensity);
            float product = 1.0f;
            product *= p2.mass;
            product *= pressure1 + pressure2;
            product *= 0.5f * p2.density;

            /* Account for direction of pressure! */
            if (nonZeroDist(_p, p2)) {
                ngl::Vec3 factor = pressureGrad(_p, p2);

                double dist = std::sqrt((_p.pos - p2.pos).lengthSquared());
                float fromQx = (_p.pos.m_x - p2.pos.m_x) / dist;
                float fromQy = (_p.pos.m_y - p2.pos.m_y) / dist;
                float fromQz = (_p.pos.m_z - p2.pos.m_z) / dist;
                ngl::Vec3 fromQ = {fromQx, fromQy, fromQz};

                /* Calculate directional derivative */
                ngl::Vec3 directional = factor * fromQ;
                sum += product * directional;
            }
            ngl::Vec3 pressure = -sum;
            return pressure;
        }
    }
}

void Simulator::draw() const
{
    glPointSize(4.0);
    m_vao->bind();
    m_vao->setData(ngl::AbstractVAO::VertexData(m_fluid.particles.size()*sizeof(Particle),m_fluid.particles[0].pos.m_x));
    m_vao->setVertexAttributePointer(0,3,GL_FLOAT,sizeof(Particle),0);
    m_vao->setVertexAttributePointer(1,3,GL_FLOAT,sizeof(Particle),6);

    m_vao->setNumIndices(m_fluid.particles.size());
    m_vao->draw();
    m_vao->unbind();


    //std::cout<<"draw \n";
//    std::cout<<"**********************\n";
//    for(auto p : m_particles)
//    {
//        std::cout<<p.pos.x<<' '<<p.pos.y<<' '<<p.pos.z<<'\n';
//    }
}

ngl::Vec3 Simulator::randomVectorOnSphere(float _radius)
{
    auto phi = ngl::Random::randomPositiveNumber(M_PI * 2.0f);
    auto costheta = ngl::Random::randomNumber();
    float u = ngl::Random::randomPositiveNumber();
    float theta = acosf(costheta);
    float r = _radius * std::cbrt(u);
    return ngl::Vec3(r*sinf(theta) * cosf(phi),
                     r*sinf(theta) * sinf(phi),
                     r*cosf(theta)
    );
}

void Simulator::getPressureDensity()
{
    float restDensity = 0.01f;
    for (auto &p1 : m_fluid.particles)
    {
        p1.density = 0.0f;

        for (auto &p2 : m_fluid.particles)
        {
            ngl::Vec3 r = p2.pos - p1.pos;
            float rSquared = r.lengthSquared();

            if (rSquared > 0.0f && rSquared < 1.0)
            {
                p1.density += p1.mass * smoothingKernel(rSquared);
            }
        }

        p1.pressure = p1.stiffness * (p1.density - restDensity);
    }
}


float Simulator::smoothingKernel(float _rSquared)
{
    // @smoothingKernel Equation
    // Author: Andrew Gibiansky
    // https://andrew.gibiansky.com/blog/physics/computational-fluid-dynamics/
    // Wg(r,h)=(315/(64π(h^9))) * ((h^2)−(r^2))^3.

    //smoothing bandwith = m_h
    float coefficient = 315/ (64.0f * M_PI * std::pow(m_h, 9));
    auto kernel = coefficient * std::pow((std::pow(m_h, 2) - _rSquared), 3);

    // Viscosity Kernel
    // Wv(r,h)=−(r^3/2h^3)+(r^2/h^2)+(h/2r)−1.

    return kernel;

}

ngl::Vec3 Simulator::pressureGrad(Particle &_p1, Particle &_p2)
{
    //Wp(r,h)=15πh6(h−r)3
    ngl::Vec3 dist = _p1.pos - _p2.pos;
    float rSquared = dist.lengthSquared();
    float r = std::sqrt(rSquared);

    float coefficient = 15.0f /(M_PI * std::pow(m_h, 6));
    auto kernel = coefficient * std::pow((m_h -r),3);
    ngl::Vec3 m_kernel = {kernel, kernel, kernel};

    return m_kernel;
}

ngl::Vec3 Simulator::smoothingKernelGrad(const ngl::Vec3 &_r)
{
    // Gradient of the smoothing kernel
    // Wg_grad(r,h) = d(Wg(r,h))/dr
    float coefficient = -945.0f / (32.0f * M_PI * std::pow(m_h, 9));
    float q = std::sqrt(_r.lengthSquared()) / m_h;
    return coefficient * _r * std::pow((1.0f - q * q), 2);
}

ngl::Vec3 Simulator::calculateViscosity(Particle &_p)
{
    // fviscosityi=μ∇2v(r⃗ i)=μ∑mjvjρj∇2W(|r⃗ i−r⃗ j|,h).
    float sumX = 0.0f;
    float sumY = 0.0f;
    float sumZ = 0.0f;
    for(auto &p1 : m_fluid.particles)
    {
        float product = 1.0f;
        product *= p1.mass;
        product *= 1/p1.density;

        ngl::Vec3 difference = ngl::Vec3(p1.velocity - _p.velocity);
        double factor = smoothingKernelLaplacian(&_p,&p1);
        sumX += product * difference.m_x * factor;
        sumY += product * difference.m_y * factor;
        sumZ += product * difference.m_z * factor;
    }

    ngl::Vec3 viscosityForce = ngl::Vec3(sumX, sumY,sumZ);
    return viscosityForce;
}


const bool Simulator::nonZeroDist(Particle &_p1, Particle &_p2)
{
//    float normsq = pow((_p1.pos.m_x- _p2.pos.m_x ), 2)+ pow((_p1.pos.m_y- _p2.pos.m_y ), 2) + pow((_p1.pos.m_z- _p2.pos.m_z), 2);
    float normSq =  std::pow(_p1.pos - _p2.pos, 2);
    float sq = m_h * m_h;
    if(normSq > sq)
    {
        return false;
    }
    return true;
}
float Simulator::smoothingKernelLaplacian(Particle &_p1, Particle &_p2)
{
    // Viscocity Kernel
    // Wg_laplacian(r,h) = d^2(Wg(r,h))/dr^2
    // −r3/2h3+r2/h2+h/2r−1.
    ngl::Vec3 dist = _p1.pos - _p2.pos;
    float rSquared = dist.lengthSquared();
    float r = std::sqrt(rSquared);

    double coefficient1 = 1.0f/(2.0f * std::pow(m_h, 3));
    double coefficient2 = 1.0f/std::pow(m_h, 2);
    float coefficient3 = m_h;
    return (-(coefficient1 * (std::pow(r, 3))) + (coefficient2 * std::pow(r,2)) + (coefficient3/(2*r))) - 1;
}

//ngl::Vec3 Simulator::applyViscosity(Particle &_p)
//{
//    ngl::Vec3 viscosityForce;
//
//    // Loop over neighbors and accumulate viscosity forces
//    for (auto &neighbor : m_fluid.particles)
//    {
//        ngl::Vec3 r = neighbor.pos - _p.pos;
//        float rSquared = r.lengthSquared();
//
//        if (nonZeroDist(_p, neighbor))
//        {
//            viscosityForce += (neighbor.dir - _p.dir) * smoothingKernelLaplacian(r) * _p.mass;
//        }
//    }

//    // Scale by viscosity coefficient
//    viscosityForce *= 0.8f;
//
//    return viscosityForce;
//}


size_t Simulator::numParticles()const
{
    return m_fluid.particles.size();
}

ngl::Vec3 Simulator::getPosition() const
{
    return m_pos;
}

