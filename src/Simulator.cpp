#include "Simulator.h"
#include <ngl/Random.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <ngl/VAOFactory.h>

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

// @ Pressure & Density Calculation
// Modified from ChatGPT result
// Query; Write a function to get pressure and density for particle based fluid animation using SPH.
//          Some scode snippets were also shared as part of the query
// Result:
// // void Emitter::calculateDensityAndPressure()
// {
//     for (auto &p1 : m_particles)
//     {
//         if (p1.alive == ParticleState::Dead)
//             continue;
//         p1.density = 0.0f;
//         for (auto &p2 : m_particles)
//         {
//             if (p2.alive == ParticleState::Dead)
//                 continue;
//             ngl::Vec3 r = p2.pos - p1.pos;
//             float rSquared = r.lengthSquared();
//             if (rSquared > 0.0f && rSquared < hSquared) // Use an appropriate smoothing radius, h
//             {
//                 p1.density += particleMass * poly6Kernel(rSquared);
//             }
//         }
//         p1.pressure = stiffness * (p1.density - restDensity);
//     }
// }

void Simulator::getPressureDensity()
{
  float restDensity = 0.01f;
  for (auto &p1 : m_particles)
    {
        if (p1.alive == ParticleState::Dead)
            continue;

        p1.density = 0.0f;

        for (auto &p2 : m_particles)
        {
            if (p2.alive == ParticleState::Dead)
                continue;

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

void Simulator::resetParticle(Particle &_p)
{
  _p.pos=m_pos;
  _p.dir = ngl::Vec3(0.0f, 1.0f, 0.0f);
  // _p.dir.m_y = std::abs(_p.dir.m_y);
  _p.colour=ngl::Random::getRandomColour3();
  _p.life = 200 + static_cast<int>(ngl::Random::randomPositiveNumber(1000.0f));
  _p.size = 0.01f;
  _p.density = 0.0f;  // Initialize density
  _p.pressure = 0.0f; // Initialize pressure
}

void Simulator::birthParticles()
{
  auto dont = static_cast<int>(ngl::Random::randomPositiveNumber(100.0f));
  if(dont <=80)
    return;
  auto births = 0 + static_cast<int>(ngl::Random::randomPositiveNumber(m_numPerFrame));
  for(size_t i=0; i<births; ++i)
  {
    for(auto &p : m_particles)
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

Simulator::Simulator(ngl::Vec3 _pos,size_t _numParticles) : m_pos{_pos}
{
    m_particles.resize(_numParticles);
    for(auto &p : m_particles)
    {
      resetParticle(p);
    }

  birthParticles();

  m_vao = ngl::VAOFactory::createVAO(ngl::simpleVAO,GL_POINTS);

}

size_t Simulator::numParticles()const
{
    return m_particles.size();
}

ngl::Vec3 Simulator::getPosition() const
{
    return m_pos;
}


void Simulator::draw() const
{
  glPointSize(4.0);
  m_vao->bind();
  m_vao->setData(ngl::AbstractVAO::VertexData(m_particles.size()*sizeof(Particle),m_particles[0].pos.m_x));
  m_vao->setVertexAttributePointer(0,3,GL_FLOAT,sizeof(Particle),0);
  m_vao->setVertexAttributePointer(1,3,GL_FLOAT,sizeof(Particle),6);

  m_vao->setNumIndices(m_particles.size());
  m_vao->draw();
  m_vao->unbind();
}

void Simulator::update(float _delta)
{
   // std::cout<<"update\n";
    float dt=_delta;
    ngl::Vec3 gravity(0.0,-9.81f,0.0);
    // find how many particles alive
    size_t numAlive = std::count_if(std::begin(m_particles),std::end(m_particles),
    [](auto p)
    {
      return p.alive == ParticleState::Alive;
    });
    if(numAlive < m_maxAlive)
    {
      birthParticles();
    }
    // add new particles if not enough
    for(auto &p : m_particles)
    {
        if(p.alive == ParticleState::Dead)
          continue;
        p.dir += gravity * dt * 0.5f;
        p.pos += p.dir * 0.5f;
        p.life -=1;
        p.size+=0.01f;
        p.size=std::clamp(p.size,0.0f,2.0f);
        if(p.pos.m_y <= 0.0f || p.life <=0)
        {
          resetParticle(p);
        }
    }
}
