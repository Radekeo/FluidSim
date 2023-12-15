#include <gtest/gtest.h>
#include "Fluid.h"
//#include "Simulator.h"

//@Fluid Tests
//Modified from ASE Labs

TEST(Particle, ctor)
{
    Particle p;

    EXPECT_FLOAT_EQ(p.pos.m_x, 0.0f);
    EXPECT_FLOAT_EQ(p.pos.m_y, 0.0f);
    EXPECT_FLOAT_EQ(p.pos.m_z, 0.0f);

    EXPECT_FLOAT_EQ(p.dir.m_x, 0.0f);
    EXPECT_FLOAT_EQ(p.dir.m_y, 0.0f);
    EXPECT_FLOAT_EQ(p.dir.m_z, 0.0f);
}


//@Fluid Tests
// Modified from chatGPT
// Prompt "what would a test for this fluid struct look like?:
// struct Fluid
//{
//    std::vector<Particle> fluidParticles;
//};"
// Result
// TEST_CASE("Fluid struct tests", "[Fluid]") {
//    // Test case 1: Default constructor
//    SECTION("Default constructor") {
//        Fluid fluid;
//        REQUIRE(fluid.fluidParticles.empty());
//    }
//    // Test case 2: Adding particles
//    SECTION("Adding particles") {
//        Fluid fluid;
//        Particle particle1(/* initialize particle here */);
//        Particle particle2(/* initialize another particle here */);
//        fluid.fluidParticles.push_back(particle1);
//        fluid.fluidParticles.push_back(particle2);
//        REQUIRE(fluid.fluidParticles.size() == 2);
//        REQUIRE(fluid.fluidParticles[0] == particle1);
//        REQUIRE(fluid.fluidParticles[1] == particle2);
//    }
//    // Add more test cases as needed...
//    // Test case n: Clearing particles
//    SECTION("Clearing particles") {
//        Fluid fluid;
//        Particle particle(/* initialize particle here */);
//        fluid.fluidParticles.push_back(particle);
//
//        fluid.fluidParticles.clear();
//        REQUIRE(fluid.fluidParticles.empty());
//    }

TEST(Fluid, ctor)
{
    Fluid f;

    // Adding particles
    Particle p1;
    Particle p2;
    Particle p3;

    f.particles.push_back(p1);
    f.particles.push_back(p2);
    f.particles.push_back(p3);

    EXPECT_EQ(f.particles.size(), 3);
}
//
//TEST(Simulator, ctor)
//{
//   Simulator s(ngl::Vec3(0.5f, 0.0f, 0.0f),100);
//   EXPECT_EQ(s.numParticles(), 20);
//   auto pos = s.getPosition();
//   EXPECT_FLOAT_EQ(pos.x, 0.5f);
//   EXPECT_FLOAT_EQ(pos.y, 0.0f);
//   EXPECT_FLOAT_EQ(pos.z, 0.0f);
//}
//
//TEST(Simulator, Smoothing)
//{
//
//}
