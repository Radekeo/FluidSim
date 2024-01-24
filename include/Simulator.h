#ifndef SIMULATOR_H_
#define SIMULATOR_H_
#include <cstdlib>
#include <vector>
#include "Fluid.h"
#include "Kernel.h"
#include <string_view>

#include <ngl/AbstractVAO.h>
#include <memory>

class Simulator
{
public :
    Simulator(ngl::Vec3 _pos, size_t _numParticles);
    size_t numParticles() const;
    ngl::Vec3 getPosition() const;

    void draw() const;
    void update(float _delta);

private :
    [[nodiscard]] ngl::Vec3 randomVectorOnSphere(float _radius = 1.0f);

    void resetParticle(Particle &_p);
    void birthParticles();

    void precalculatePressureDensity();
    void calculateAcceleration(Particle &_p);
    ngl::Vec3 calculatePressure(Particle &_p);
    ngl::Vec3 calculateViscosity(Particle &_p);
    ngl::Vec3 calculateExternalForce(Particle &_p);
    ngl::Vec3 calculateSurfaceTension(Particle &_p);
    void applyCollisions(Particle &_p);

    ngl::Vec3 applyViscosity(Particle &_p);

    Fluid m_fluid;
    ngl::Vec3 m_pos;
    ngl::Vec3 m_emitDir={0.0f,1.0f,0.0f};
    ngl::Vec3 m_gravity= {0.0f,-9.8f,0.0f};
    float m_spread=0.5f;
    // max alive at one time
    size_t m_maxAlive = 4000;
    // max birthed at one time
    size_t m_numPerFrame =100;
    std::unique_ptr<ngl::AbstractVAO> m_vao;

    float m_viscosity = 10.0f;
    float m_surfaceTension = 1000.0f;
    Kernel m_kernel;
};

#endif