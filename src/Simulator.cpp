#include "Simulator.h"
#include <ngl/Random.h>
#include <iostream>
#include <fstream>
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
    // Set initial direction to point downward
//    _p.dir = ngl::Vec3(0.0f, -1.0f, 0.0f);
    _p.dir = randomVectorOnSphere(m_spread);

//  _p.dir.m_y = std::abs(_p.dir.m_y);
    _p.colour=ngl::Random::getRandomColour3();
    _p.life = 100 + static_cast<int>(ngl::Random::randomPositiveNumber(1000.0f));
    _p.size = 0.01f;
    _p.density = 0.0f;  // Initialize density
    _p.pressure = 0.0f; // Initialize pressure

    // Randomize position
    float random_x = ngl::Random::randomNumber(0.5f);
    _p.pos.m_x += random_x;


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
    ngl::Vec3 gravity(0.0,-9.81f,0.0);

    getPressureDensity(); // update pressure and density


    for (auto &p : m_fluid.particles)
    {
        // Apply gravity
        p.dir += gravity * dt * 0.5f;

        // Apply pressure force
        ngl::Vec3 pressureForce;
        for (auto &p2 : m_fluid.particles)
        {
            ngl::Vec3 r = p2.pos - p.pos;
            float rSquared = r.lengthSquared();

            if (rSquared > 0.0f && rSquared < 1.0)
            {
                pressureForce += (p.pressure + p2.pressure) / (2 * p2.density) * smoothingKernelGrad(r) * p2.mass;
            }
        }
        p.dir += pressureForce;

        // Apply viscosity force (you need to implement this function)
        ngl::Vec3 viscosityForce = applyViscosity(p);
        p.dir += viscosityForce;

        // Update particle position
        p.pos += p.dir * 0.5f;
        p.life -= 1;
        p.size += 0.01f;
        p.size = std::clamp(p.size, 0.0f, 2.0f);

        // Clamp particles so they don't fall below the grid
        p.pos.m_y = std::max(p.pos.m_y, 0.0f);
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

    float h = 1.0f; //smoothing bandwith
    float coefficient = 315/ (64.0f * M_PI * std::pow(h, 9));
    auto kernel = coefficient * std::pow((std::pow(h, 2) - _rSquared), 3);

    return kernel;

}

ngl::Vec3 Simulator::smoothingKernelGrad(const ngl::Vec3 &_r)
{
    // Gradient of the smoothing kernel
    // Wg_grad(r,h) = d(Wg(r,h))/dr
    float h = 1.0f; // smoothing bandwidth
    float coefficient = -945.0f / (32.0f * M_PI * std::pow(h, 9));
    float q = std::sqrt(_r.lengthSquared()) / h;
    return coefficient * _r * std::pow((1.0f - q * q), 2);
}

float Simulator::smoothingKernelLaplacian(const ngl::Vec3 &_r)
{
    // Laplacian of the smoothing kernel
    // Wg_laplacian(r,h) = d^2(Wg(r,h))/dr^2
    float h = 1.0f; // smoothing bandwidth
    float coefficient = -945.0f / (32.0f * M_PI * std::pow(h, 9));
    float q = std::sqrt(_r.lengthSquared()) / h;
    return coefficient * (3.0f * q * q - 1.0f) * std::pow((1.0f - q * q), 2);
}

ngl::Vec3 Simulator::applyViscosity(const Particle &_p)
{
    ngl::Vec3 viscosityForce;

    // Loop over neighbors and accumulate viscosity forces
    for (const auto &neighbor : m_fluid.particles)
    {
        ngl::Vec3 r = neighbor.pos - _p.pos;
        float rSquared = r.lengthSquared();

        if (rSquared > 0.0f && rSquared < 1.0)
        {
            viscosityForce += (neighbor.dir - _p.dir) * smoothingKernelLaplacian(r) * _p.mass;
        }
    }

    // Scale by viscosity coefficient
    viscosityForce *= 1.0;

    return viscosityForce;
}


size_t Simulator::numParticles()const
{
    return m_fluid.particles.size();
}

ngl::Vec3 Simulator::getPosition() const
{
    return m_pos;
}