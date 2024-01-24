#include "Simulator.h"
#include <ngl/Random.h>
#include <iostream>
//#include <fstream>
#include <algorithm>
#include <ngl/VAOFactory.h>

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
//   _p.dir = ngl::Vec3(0.0f, -1.0f, 0.0f);
   _p.dir = randomVectorOnSphere(m_spread);

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
            if(p.life <= 0)
            {
                resetParticle(p);
                break;
            }
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
//
//
    for (auto &p : m_fluid.particles)
    {
        // apply acceleration
        p.velocity = p.acceleration * _delta;
        // Apply gravity
//        p.dir += gravity * dt * 0.5f;
        // Reset Acceleration to Zero
        p.acceleration.m_x = 0.0f;
        p.acceleration.m_y = 0.0f;
        p.acceleration.m_z = 0.0f;
    }

    for(auto &p : m_fluid.particles) {
        p.pos += p.velocity * dt;
        p.pos.m_y = std::max(p.pos.m_y, 0.0f);
    }
    // Check for collisions
    for(auto &p : m_fluid.particles) {
        applyCollisions(p);
    }
}

//void Simulator::update(float _delta)
//{
//    // std::cout<<"update\n";
//    float dt=_delta;
////    ngl::Vec3 gravity(0.0,-9.81f,0.0);
//
//    precalculatePressureDensity(); // update pressure and density
//    for (auto &p : m_fluid.particles)
//    {
//        // Apply gravity
//        p.dir += m_gravity * dt * 0.5f;
//
//        // Apply pressure force
//        ngl::Vec3 pressureForce = calculatePressure(p);
//        p.dir += pressureForce;
//
//        // Apply viscosity force
//        ngl::Vec3 viscosityForce = applyViscosity(p);
//        p.dir += viscosityForce;
//
//        // Apply Surface Tension
//        ngl::Vec3 surfaceForce = calculateSurfaceTension(p);
//        p.dir += surfaceForce;
//
//        // Apply External Force
//        ngl::Vec3 externalForce = calculateExternalForce(p);
//        p.dir += externalForce;
//
//        // Update particle position
//        p.pos += p.dir * 0.5f;
//        p.life -= 1;
//
//        // Clamp particles so they don't fall below the grid
//        p.pos.m_y = std::max(p.pos.m_y, 0.0f);
//    }
//}

void Simulator::calculateAcceleration(Particle &_p)
{
    precalculatePressureDensity();
    /* Calculate individual force terms */
    ngl::Vec3 pressureForce = calculatePressure(_p);
    ngl::Vec3 viscosityForce = calculateViscosity(_p);
    ngl::Vec3 externalForce = calculateExternalForce(_p);
    ngl::Vec3 surfaceForce = calculateSurfaceTension(_p);

    ngl::Vec3 totalForce = pressureForce + viscosityForce + externalForce + surfaceForce;
    _p.acceleration = totalForce / _p.density;
}

void Simulator::precalculatePressureDensity()
{
    for (auto &p1 : m_fluid.particles) {
        p1.density = 0.0f;

        for (auto &p2: m_fluid.particles) {
            ngl::Vec3 dist = p2.pos - p1.pos;
            float rSquared = dist.lengthSquared();
            float r = std::sqrt(rSquared);

            if (r>0) {
                p1.density += p1.mass * m_kernel.smoothingKernel(rSquared);
            }
        }
        // Update pressure for each particle
        p1.pressure = p1.stiffness * (p1.density - p1.restDensity);
    }
}

ngl::Vec3 Simulator::applyViscosity(Particle &_p)
{
    ngl::Vec3 viscosityForce;

    // Loop over neighbors and accumulate viscosity forces
    for (auto &neighbor : m_fluid.particles)
    {
        ngl::Vec3 dist = neighbor.pos - _p.pos;
        float rSquared = dist.lengthSquared();
        float  r = std::sqrt(rSquared);

        if (r > 0)
        {
            viscosityForce += (neighbor.dir - _p.dir) * m_kernel.smoothingKernelLaplacian(_p, neighbor) * _p.mass;
        }
    }

    // Scale by viscosity coefficient
    viscosityForce *= 0.8f;

    return viscosityForce;
}
ngl::Vec3 Simulator::calculateViscosity(Particle &_p)
{
    // fviscosityi=μ∇2v(r⃗ i)=μ∑mjvjρj∇2W(|r⃗ i−r⃗ j|,h).
    ngl::Vec3 sum;
    for(auto &p1 : m_fluid.particles)
    {
        float product = 1.0f;
        product *= p1.mass;
        product *= 1/p1.density;

        ngl::Vec3 difference = p1.velocity - _p.velocity;
        float factor = m_kernel.smoothingKernelLaplacian(_p,p1);
        sum += product * difference* factor;
    }

    ngl::Vec3 viscosityForce = sum;
    return viscosityForce;
}

ngl::Vec3 Simulator::calculateExternalForce(Particle &_p)
{
    return m_gravity * _p.mass;
}

ngl::Vec3 Simulator::calculateSurfaceTension(Particle &_p)
{
    ngl::Vec3 sum = {0.0f, 0.0f, 0.0f};

    for(auto &p1 : m_fluid.particles)
    {
        float product = 1.0f;
        product *= p1.mass;
        product *= 1/p1.density;

        ngl::Vec3 factor = m_kernel.smoothingKernelGrad(_p.pos - p1.pos);
        sum += product * factor;
    }

    double magN = std::sqrt(std::pow(sum.m_x, 2) + std::pow(sum.m_y, 2) + std::pow(sum.m_z, 2));
    float nThreshold = 0.25f;
    if(magN >= nThreshold){
        /* Calculate laplacian of color */
        float colorLaplacian = 0.0f;
        for(auto &p1 : m_fluid.particles){
            float product = 1;
            product *= p1.mass;
            product *= 1/p1.density;

            double factor = m_kernel.smoothingKernelLaplacian(_p, p1);
            colorLaplacian += product * factor;
        }

        /* Calculate surface tension force */
        float surfaceForce = 1.0f;
        surfaceForce *= - m_surfaceTension;
        surfaceForce *= colorLaplacian;
        surfaceForce /= magN;
        return surfaceForce * sum;
    } else {
        ngl::Vec3 zeroVec = {0.0f, 0.0f, 0.0f};
        return zeroVec;
    }
}

ngl::Vec3 Simulator::calculatePressure(Particle &_p)
{
    ngl::Vec3 pressureForce;
    //p=kρ
    // stiffness represents constant k
    for (auto &p2: m_fluid.particles)
    {
        ngl::Vec3 dist = p2.pos - _p.pos;
        float rSquared = dist.lengthSquared();
        if (rSquared > 0.0f && rSquared < 1.0)
        {
            // @Pressure Equation
            // Author: Andrew Gibiansky
            // https://andrew.gibiansky.com/blog/physics/computational-fluid-dynamics/
            //fpressurei=−∇p(r⃗ )=−∑(mj((pi+pj)/2ρj))∇W(|r⃗ i−r⃗ j|,h).
            pressureForce += (((_p.pressure + p2.pressure) / (2 * p2.density)) * p2.mass) * m_kernel.pressureGrad(_p, p2);
        }
    }
    return pressureForce;
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
    for(auto &p : m_fluid.particles)
    {
        std::cout<<p.pos.m_x<<' '<<p.pos.m_y<<' '<<p.pos.m_z<<'\n';
    }
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

void Simulator::applyCollisions(Particle &_p) {
    float dampingFactor = 0.8f;
    if (_p.pos.m_y <= 40) {
        _p.pos.m_y = 40.0f;
        _p.velocity.m_y = -_p.velocity.m_y * dampingFactor;
    }
    if (_p.pos.m_x < 50 ) {
        _p.pos.m_x = 50.0f;
        _p.velocity.m_x = -_p.velocity.m_x * dampingFactor;
    }
    if (_p.pos.m_x > 350) {
        _p.pos.m_x = 350;
        _p.velocity.m_x = -_p.velocity.m_x  * dampingFactor;
    }
    if (_p.pos.m_z < 50 ) {
        _p.pos.m_z = 50.0f;
        _p.velocity.m_z = -_p.velocity.m_z * dampingFactor;
    }
    if (_p.pos.m_x > 350) {
        _p.pos.m_x = 350;
        _p.velocity.m_z = -_p.velocity.m_z  * dampingFactor;
    }
}


size_t Simulator::numParticles()const
{
    return m_fluid.particles.size();
}

ngl::Vec3 Simulator::getPosition() const
{
    return m_pos;
}

