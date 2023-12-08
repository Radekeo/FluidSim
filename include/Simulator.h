#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <cstdlib>
#include <vector>
#include "Fluid.h"
#include <string_view>0

#include <ngl/AbstractVAO.h>
#include <memory>

class Simulator
{
public :
    Simulator(ngl::Vec3 _pos,size_t _numParticles);
    size_t numParticles() const;
    ngl::Vec3 getPosition() const;
    void draw() const;
    void update(float _delta);
    void getPressureDensity();
    float smoothingKernel(float _r); 
private :
    [[nodiscard]] ngl::Vec3 randomVectorOnSphere(float _radius = 2.0f);
    void resetParticle(Particle &_p);
    void birthParticles();
    std::vector<Particle> m_particles;
    ngl::Vec3 m_pos;
    ngl::Vec3 m_emitDir={0.5f,0.5f,0.0f};
    float m_spread=5.0f;
    // max alive at one time
    size_t m_maxAlive =5000;
    // max birthed at one time
    size_t m_numPerFrame =200;
    std::unique_ptr<ngl::AbstractVAO> m_vao;
};


#endif //SIMULATOR_H_
