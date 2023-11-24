#include <gtest/gtest.h>
#include "Fluid.h"
#include "Vec3.h"


TEST(Fluid,ctor)
{
    FluidParticle fp;
    EXPECT_FLOAT_EQ(fp.pos.x, 0.0f);
    EXPECT_FLOAT_EQ(fp.pos.y, 0.0f);
    EXPECT_FLOAT_EQ(fp.pos.z, 0.0f);

    EXPECT_FLOAT_EQ(fp.dir.x, 0.0f);
    EXPECT_FLOAT_EQ(fp.dir.y, 0.0f);
    EXPECT_FLOAT_EQ(fp.dir.z, 0.0f);

    EXPECT_FLOAT_EQ(fp.colour.x, 0.0f);
    EXPECT_FLOAT_EQ(fp.colour.y, 0.0f);
    EXPECT_FLOAT_EQ(fp.colour.z, 0.0f);

    EXPECT_FLOAT_EQ(fp.size, 1.0f);
    EXPECT_EQ(fp.life, 100);
}

TEST(Vec3, ctor)
{
    Vec3 v;
    EXPECT_FLOAT_EQ(v.x, 0.0f);
    EXPECT_FLOAT_EQ(v.y, 0.0f);
    EXPECT_FLOAT_EQ(v.z, 0.0f);
}

TEST(Simulator, ctor)
{
    // TODO
}


TEST(Domain, ctor)
{
    // TODO
}

