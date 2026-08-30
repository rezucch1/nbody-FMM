/**
 * @file test_tensor.cpp
 * @brief GoogleTest unit tests for Tensor mathematical vector class.
 */

#include <gtest/gtest.h>
#include "math_utils/tensor.hpp"

TEST(TensorTest, DimensionConstructor) {
    Tensor t(3);
    EXPECT_EQ(t.dim, 3);
    EXPECT_DOUBLE_EQ(t[0], 0.0);
    EXPECT_DOUBLE_EQ(t[1], 0.0);
    EXPECT_DOUBLE_EQ(t[2], 0.0);
}

TEST(TensorTest, InitializerListConstructor) {
    Tensor t{1.0, 2.0, 3.0};
    EXPECT_EQ(t.dim, 3);
    EXPECT_DOUBLE_EQ(t[0], 1.0);
    EXPECT_DOUBLE_EQ(t[1], 2.0);
    EXPECT_DOUBLE_EQ(t[2], 3.0);
}

TEST(TensorTest, CopyAndMoveConstructors) {
    Tensor original{1.0, 2.0, 3.0};
    Tensor copy_tensor(original);
    EXPECT_TRUE(copy_tensor == original);

    Tensor moved_tensor(std::move(copy_tensor));
    EXPECT_TRUE(moved_tensor == original);
}

TEST(TensorTest, IndexingAndEquality) {
    Tensor t1{4.0, 5.0, 6.0};
    Tensor t2{4.0, 5.0, 6.0};
    Tensor t3{4.0, 5.0, 6.1};

    EXPECT_TRUE(t1 == t2);
    EXPECT_FALSE(t1 == t3);

    t1[1] = 10.0;
    EXPECT_DOUBLE_EQ(t1[1], 10.0);
}

TEST(TensorTest, VectorAdditionAndSubtraction) {
    Tensor a{1.0, 2.0, 3.0};
    Tensor b{4.0, 5.0, 6.0};

    Tensor sum = a + b;
    EXPECT_DOUBLE_EQ(sum[0], 5.0);
    EXPECT_DOUBLE_EQ(sum[1], 7.0);
    EXPECT_DOUBLE_EQ(sum[2], 9.0);

    Tensor diff = b - a;
    EXPECT_DOUBLE_EQ(diff[0], 3.0);
    EXPECT_DOUBLE_EQ(diff[1], 3.0);
    EXPECT_DOUBLE_EQ(diff[2], 3.0);

    Tensor neg = -a;
    EXPECT_DOUBLE_EQ(neg[0], -1.0);
    EXPECT_DOUBLE_EQ(neg[1], -2.0);
    EXPECT_DOUBLE_EQ(neg[2], -3.0);
}

TEST(TensorTest, InPlaceOperators) {
    Tensor a{1.0, 2.0, 3.0};
    Tensor b{4.0, 5.0, 6.0};

    a += b;
    EXPECT_DOUBLE_EQ(a[0], 5.0);
    EXPECT_DOUBLE_EQ(a[1], 7.0);
    EXPECT_DOUBLE_EQ(a[2], 9.0);

    a -= b;
    EXPECT_DOUBLE_EQ(a[0], 1.0);
    EXPECT_DOUBLE_EQ(a[1], 2.0);
    EXPECT_DOUBLE_EQ(a[2], 3.0);
}

TEST(TensorTest, ScalarScalingAndDivision) {
    Tensor a{1.0, -2.0, 3.0};

    Tensor r_scaled = a * 2.5;
    EXPECT_DOUBLE_EQ(r_scaled[0], 2.5);
    EXPECT_DOUBLE_EQ(r_scaled[1], -5.0);
    EXPECT_DOUBLE_EQ(r_scaled[2], 7.5);

    Tensor l_scaled = 2.5 * a;
    EXPECT_TRUE(l_scaled == r_scaled);

    Tensor div = r_scaled / 2.5;
    EXPECT_TRUE(div == a);
}

TEST(TensorTest, Norms) {
    Tensor a{3.0, 4.0, 0.0};
    EXPECT_DOUBLE_EQ(a.squared_norm(), 25.0);
    EXPECT_DOUBLE_EQ(a.norm(), 5.0);

    Tensor b{1.0, 2.0, 2.0};
    EXPECT_DOUBLE_EQ(b.squared_norm(), 9.0);
    EXPECT_DOUBLE_EQ(b.norm(), 3.0);
}
