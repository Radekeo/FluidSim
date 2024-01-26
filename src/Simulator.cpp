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
    float y_interval = -0.1f;
    // Set initial direction to fall downward
   _p.dir = ngl::Vec3(0.0f, -1.0f, 0.0f);
   _p.life = 200 + static_cast<int>(ngl::Random::randomPositiveNumber(1000.0f));

//   _p.dir = randomVectorOnSphere(1.0f);
//   _p.dir = randomFlow(m_spread);

    // Randomize x position
    float random_x = ngl::Random::randomNumber(10.0f);
    float diff1 = abs(_p.pos.m_x - 10);
    float diff2 = abs(_p.pos.m_x - (-10));

    if (diff1 <= diff2)
    {
        _p.pos.m_x += random_x;
    }
    else
    {
        _p.pos.m_x += -random_x;
    }
    _p.pos.m_y += y_interval;
    m_pos.m_y = _p.pos.m_y;
}

void Simulator::birthParticles()
{
    auto dont = static_cast<int>(ngl::Random::randomPositiveNumber(100.0f));
    if(dont <= 10)
        return;
    auto births = 0 + static_cast<int>(ngl::Random::randomPositiveNumber(m_numPerFrame));
    for(size_t i=0; i<births; ++i)
    {
        for(auto &p : m_fluid.particles)
        {
            if(p.alive == ParticleState::Dead)
            {
                resetParticle(p);
                p.alive=ParticleState::Alive;
                break;
            }
        }
    }
}

void Simulator::addParticles()
{
    for (auto &p : m_fluid.particles)
    {
        if (p.alive == ParticleState::Dead)
        {
            resetParticle(p);
            p.alive = ParticleState::Alive;
        }
    }
}


void Simulator::update(float _delta)
{
    // std::cout<<"update\n";
    float dt=_delta;
//    ngl::Vec3 gravity(0.0,-9.81f,0.0);
    size_t numAlive = std::count_if(std::begin(m_fluid.particles),std::end(m_fluid.particles),
                                    [](auto p)
                                    {
                                        return p.alive == ParticleState::Alive;
                                    });
    if(numAlive < m_maxAlive)
    {
        birthParticles();
    }

    precalculateDensity(); // update density for each particle
    for (auto &p : m_fluid.particles)
    {
        p.dir = randomVectorOnSphere(m_spread);
        if (p.alive == ParticleState::Alive)
        {
            // Check if particle crosses boundary on X-axis
            if (p.pos.m_x <= -40.0f || p.pos.m_x >= 40.0f)
            {
                float diff1 = abs(p.pos.m_x - 40);
                float diff2 = abs(p.pos.m_x - (-40));
                // Reset particle position to stay within the boundary
                float new_x = std::max(0.0f, ngl::Random::randomNumber(40));
                if (diff1 >= diff2)
                {
                    p.pos.m_x = new_x;
                }
                else
                {
                    p.pos.m_x = -new_x;
                }
            }
        }
        if(p.alive == ParticleState::Dead)
            continue;

        // Apply Force Terms
        //Pressure Force
        ngl::Vec3 pressureForce = calculatePressure(p);
        // Apply viscosity force
        ngl::Vec3 viscosityForce = applyViscosity(p);
        // Apply External Force
        ngl::Vec3 externalForce = calculateExternalForce(p);
        // Apply Surface Tension
        ngl::Vec3 surfaceForce = calculateSurfaceTension(p);

        ngl::Vec3 totalForce = pressureForce + viscosityForce + surfaceForce + externalForce;
        p.dir += totalForce/p.density;

        // Update particle position
        p.pos += p.dir * 0.5f;
        p.life -= 1;

//        Clamp particles so they don't fall below the grid
        p.pos.m_y = std::max(p.pos.m_y, -3.0f);

        if(p.pos.m_x <= -50.0f || p.pos.m_x > 50.0f || p.life <=0.0f)
        {
            resetParticle(p);
        }
    }
    addParticles();

    for (auto &p : m_fluid.particles)
    {
        applyCollisions(p);
    }

   m_spread += 0.1;
}

void Simulator::precalculateDensity()
{
    for (auto &p1 : m_fluid.particles) {

        for (auto &p2: m_fluid.particles) {
            ngl::Vec3 dist = p2.pos - p1.pos;
            float rSquared = dist.lengthSquared();
            float r = std::sqrt(rSquared);

            if (r>0 && r < 1) {
                p1.density += p1.mass * m_kernel.smoothingKernel(rSquared);
            }
        }
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

        if (r > 0 && r < 1)
        {
            viscosityForce += (neighbor.dir - _p.dir) * m_kernel.smoothingKernelLaplacian(_p, neighbor) * _p.mass;
        }
    }

    // Scale by viscosity coefficient
    viscosityForce *= m_viscosity;

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

        float r = std::sqrt((_p.pos- p1.pos).lengthSquared());
        if (r > 0 && r < 1)
        {
            ngl::Vec3 factor = m_kernel.smoothingKernelGrad(_p.pos - p1.pos);
            sum += product * factor;
        }

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
        ngl::Vec3 zeroVec;
        return zeroVec;
    }
}

ngl::Vec3 Simulator::calculatePressure(Particle &_p)
{
    ngl::Vec3 pressureForce;
    //p=kρ
    // stiffness represents constant k
    float pressure1 = _p.stiffness * (_p.density - _p.restDensity);
    for (auto &p2: m_fluid.particles)
    {
        float pressure2 = p2.stiffness * (p2.density - p2.restDensity);
        ngl::Vec3 dist = p2.pos - _p.pos;
        float rSquared = dist.lengthSquared();
        float r = std::sqrt(rSquared);
        if (r > 0.0f && r < 1.0f)
        {
            // @Pressure Equation
            // Author: Andrew Gibiansky
            // https://andrew.gibiansky.com/blog/physics/computational-fluid-dynamics/
            //fpressurei=−∇p(r⃗ )=−∑(mj((pi+pj)/2ρj))∇W(|r⃗ i−r⃗ j|,h).
            pressureForce += (((pressure1 + pressure2) / (2 * p2.density)) * p2.mass) * m_kernel.pressureGrad(_p, p2);
        }
    }
    return pressureForce;
}

void Simulator::draw() const
{
    glPointSize(8.0);
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
        std::cout<<m_viscosity<<'\n';
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
ngl::Vec3 Simulator::randomFlow(float _radius)
{
    float phi, theta, r;
    float u = ngl::Random::randomPositiveNumber();

    // Adjust behavior based on radius
    if (_radius <= 1.0f) {
        // Flowing out of the tap
    phi = ngl::Random::randomPositiveNumber(M_PI * 2.0f);
    theta = acosf(ngl::Random::randomNumber());
    r = _radius * std::cbrt(u);
    } else if (_radius <= 3.0f) {
        // Transition to wider space
        phi = ngl::Random::randomPositiveNumber(M_PI * 4.0f);  // Adjust angle range
        theta = acosf(ngl::Random::randomNumber());
        r = _radius * std::cbrt(u);
    } else {
        // Freely flowing
        phi = ngl::Random::randomPositiveNumber(M_PI * 2.0f);
        theta = acosf(ngl::Random::randomNumber());
        r = _radius * std::cbrt(u);
    }

    // Convert spherical coordinates to Cartesian coordinates
    return ngl::Vec3(r * sinf(theta) * cosf(phi),
                     r * sinf(theta) * sinf(phi),
                     r * cosf(theta)
    );
}
void Simulator::applyCollisions(Particle &_p) {
    float dampingFactor = 1.2f;

    // Check for collisions with boundaries
    if (_p.pos.m_y <= -3.1f) {
        // Collision with ground
        _p.pos.m_y = -3.0f; // Set particle position to ground level
        _p.velocity.m_y = -_p.velocity.m_y * dampingFactor; // Reverse and dampen y velocity
    }
    if (abs(_p.pos.m_x) > 40.0f) {
        // Collision with x-axis boundary
        _p.pos.m_x = 10.0f * (_p.pos.m_x > 0 ? 1 : -1); // Clamp particle position within boundary
        _p.velocity.m_x = -_p.velocity.m_x * dampingFactor; // Reverse and dampen x velocity
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
