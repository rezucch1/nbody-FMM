#pragma once
#include <memory>
#include <initializer_list>
#include <ostream>

/**
 * @file tensor.hpp
 * @brief Multi-dimensional vector/tensor representation and linear algebra operations.
 */

/**
 * @class Tensor
 * @brief Represents a dynamic N-dimensional vector with basic vector space operations.
 *
 * Provides operations such as vector addition, scalar multiplication, norm calculation,
 * and element indexing.
 */
class Tensor {
  public:
    /**
     * @brief Default constructor. Constructs an uninitialized or zero-dimensional Tensor.
     */
    Tensor();

    /**
     * @brief Constructs a Tensor initialized with elements from an initializer list.
     * @param init Initializer list of double values.
     */
    Tensor(std::initializer_list<double> init);

    /**
     * @brief Constructs a zero-initialized Tensor of specified dimension.
     * @param dim Spatial dimension of the tensor.
     */
    Tensor(size_t dim);

    /**
     * @brief Copy constructor.
     * @param other Tensor to copy from.
     */
    Tensor(const Tensor &other) noexcept;

    /**
     * @brief Move constructor.
     * @param other Tensor to move from.
     */
    Tensor(Tensor &&other) noexcept;

    /**
     * @brief Copy assignment operator.
     * @param other Tensor to copy.
     * @return Reference to this Tensor.
     */
    Tensor &operator=(const Tensor &other);

    /**
     * @brief Move assignment operator.
     * @param other Tensor to move.
     * @return Reference to this Tensor.
     */
    Tensor &operator=(const Tensor &&other);

    /**
     * @brief Access element at index.
     * @param idx Index of element (0 to dim-1).
     * @return Reference to element.
     */
    double &operator[](size_t idx) const;

    /**
     * @brief Vector addition (const reference).
     * @param other Tensor to add.
     * @return Resulting Tensor sum.
     */
    Tensor operator+(const Tensor &other) const;

    /**
     * @brief Vector addition (rvalue reference).
     * @param other Tensor to add.
     * @return Resulting Tensor sum.
     */
    Tensor operator+(const Tensor &&other) const;

    /**
     * @brief Unary negation operator.
     * @return Negated Tensor.
     */
    Tensor operator-() const;

    /**
     * @brief Vector subtraction (const reference).
     * @param other Tensor to subtract.
     * @return Resulting Tensor difference.
     */
    Tensor operator-(const Tensor &other) const;

    /**
     * @brief Vector subtraction (rvalue reference).
     * @param other Tensor to subtract.
     * @return Resulting Tensor difference.
     */
    Tensor operator-(const Tensor &&other) const;

    /**
     * @brief In-place vector addition (const reference).
     * @param other Tensor to add.
     * @return Reference to this Tensor.
     */
    Tensor &operator+=(const Tensor &other);

    /**
     * @brief In-place vector addition (rvalue reference).
     * @param other Tensor to add.
     * @return Reference to this Tensor.
     */
    Tensor &operator+=(const Tensor &&other);

    /**
     * @brief In-place vector subtraction (const reference).
     * @param other Tensor to subtract.
     * @return Reference to this Tensor.
     */
    Tensor &operator-=(const Tensor &other);

    /**
     * @brief In-place vector subtraction (rvalue reference).
     * @param other Tensor to subtract.
     * @return Reference to this Tensor.
     */
    Tensor &operator-=(const Tensor &&other);

    /**
     * @brief In-place scalar multiplication.
     * @param scalar Multiplier factor.
     * @return Reference to this Tensor.
     */
    Tensor &operator*=(double scalar);

    /**
     * @brief Scalar multiplication.
     * @param scalar Multiplier factor.
     * @return Scaled Tensor.
     */
    Tensor operator*(double scalar) const;

    /**
     * @brief Scalar division.
     * @param scalar Division factor.
     * @return Scaled Tensor.
     */
    Tensor operator/(double scalar) const;

    /**
     * @brief In-place scalar division.
     * @param scalar Division factor.
     * @return Reference to this Tensor.
     */
    Tensor &operator/=(double scalar);

    /**
     * @brief Element-wise equality comparison.
     * @param other Tensor to compare with.
     * @return True if dimensions and components match.
     */
    bool operator==(const Tensor &other) const;

    /**
     * @brief Stream insertion operator for Tensor formatting.
     * @param os Output stream.
     * @param tensor Tensor to print.
     * @return Output stream reference.
     */
    friend std::ostream &operator<<(std::ostream &os, const Tensor &tensor);

    /**
     * @brief Computes squared Euclidean norm \f$ \|x\|^2 = \sum_{i=1}^{d} x_i^2 \f$.
     * @return Squared norm value.
     */
    double squared_norm() const;

    /**
     * @brief Computes Euclidean norm \f$ \|x\| = \sqrt{\sum_{i=1}^{d} x_i^2} \f$.
     * @return Euclidean norm value.
     */
    double norm() const;

    size_t dim; /**< Dimension of the tensor vector space. */

  protected:
    std::unique_ptr<double[]> data; /**< Dynamic array storing vector components. */
};

/**
 * @brief Left scalar multiplication for Tensor.
 * @param scalar Scalar value.
 * @param vect Target Tensor vector.
 * @return Scaled Tensor.
 */
Tensor operator*(double scalar, const Tensor &vect);

/**
 * @brief Left scalar multiplication for rvalue Tensor.
 * @param scalar Scalar value.
 * @param vect Target rvalue Tensor.
 * @return Scaled Tensor.
 */
Tensor operator*(double scalar, const Tensor &&vect);

/**
 * @brief Stream insertion operator for rvalue Tensor formatting.
 * @param os Output stream.
 * @param tensor Rvalue Tensor.
 * @return Output stream reference.
 */
std::ostream &operator<<(std::ostream &os, const Tensor &&tensor);
