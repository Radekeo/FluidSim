#ifndef SIMULATOR_H_
#define SIMULATOR_H_
#include <cstdlib>
#include <vector>
#include "Fluid.h"
#include <string_view>

#include <ngl/AbstractVAO.h>
#include <memory>

class Simulator
{
public :
    Simulator(ngl::Vec3 _pos, size_t _numParticles);
    size_t numParticles() const;
    ngl::Vec3 getPosition() const;
    float smoothingKernel(float _r);
    float smoothingKernelGrad(const ngl::Vec3 &_r);
    float smoothingKernelLaplacian(Particle &_p1, Particle &_p2);

    void draw() const;
    void update(float _delta);
    void getPressureDensity();
    const bool nonZeroDist(Particle &_p1, Particle &_p2);
    ngl::Vec3 pressureGrad(Particle &_p1, Particle &_p2);
    ngl::Vec3 calculatePressure(Particle &_p);
    ngl::Vec3 calculateViscosity(Particle &_p);


private :
    [[nodiscard]] ngl::Vec3 randomVectorOnSphere(float _radius = 1.0f);

    void resetParticle(Particle &_p);
    void calculateAcceleration(Particle &_p);
    void birthParticles();


    void precalculateDensities();


//    ngl::Vec3 applyViscosity(Particle &_p);
    ngl::Vec3 calculateViscosity(Particle &_p);
    ngl::Vec3 gravity= {0.0f,-9.8f,0.0f};
    Fluid m_fluid;
    ngl::Vec3 m_pos;
    ngl::Vec3 m_emitDir={0.0f,1.0f,0.0f};
    float m_spread=0.5f;
    // max alive at one time
    size_t m_maxAlive = 4000;
    // max birthed at one time
    size_t m_numPerFrame =500;
    std::unique_ptr<ngl::AbstractVAO> m_vao;

    float m_viscosity = 10.0f;
    float m_surfaceTension = 1000.0f;
    float m_h = 10.0f;
};


#endif