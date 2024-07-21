/**
* @Function: GNSS Wheel Velocity error used for loosely integration with other sensors
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#pragma once

#pragma diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
// Eigen 3.2.7 uses std::binder1st and std::binder2nd which are deprecated since c++11
// Fix is in 3.3 devel (http://eigen.tuxfamily.org/bz/show_bug.cgi?id=872).
#include <ceres/ceres.h>
#include <Eigen/Core>
#pragma diagnostic pop

#include "gici/estimate/error_interface.h"
#include "gici/gnss/geodetic_coordinate.h"
#include "gici/gnss/gnss_types.h"

namespace gici {

// wheel velocity error with slip factor
// The candidate parameter setups are:
// Group 1: P1. body velocity in body frame, P2. slip factor (1)
// Group 2: P1. body pose in ENU (7), P2. body velocity in body frame, P3. relative 
//              position from body to receiver in body frame (3), P4. slip factor (1)
template<int... Ns>
class WheelVelocityError :
    public ceres::SizedCostFunction<
    3 /* number of residuals */,
    Ns ... /* parameter blocks */>,
    public ErrorInterface
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// \brief The base class type.
  typedef ceres::SizedCostFunction<3, Ns ...> base_t;

  /// \brief Number of residuals (3).
  static const int kNumResiduals = 3;

  /// \brief The information matrix type (3x3).
  typedef Eigen::Matrix<double, 3, 3> information_t;

  /// \brief The covariance matrix type (same as information).
  typedef Eigen::Matrix<double, 3, 3> covariance_t;

  /// \brief Default constructor.
  WheelVelocityError();

  /// \brief Construct with measurement and information matrix
  /// @param[in] velocity_measurement The wheel speed measurement.
  /// @param[in] information The information matrix.
  WheelVelocityError(const Eigen::Vector3d& velocity_measurement,
                   const Eigen::Matrix3d& information);

  /// \brief Construct with measurement and information matrix
  /// @param[in] velocity_measurement The wheel speed measurement.
  /// @param[in] information The information matrix.
  /// @param[in] angular_velocity The angular velocity.
  WheelVelocityError(const Eigen::Vector3d& velocity_measurement,
                   const Eigen::Matrix3d& information,
                   const Eigen::Vector3d& angular_velocity);

  /// \brief Trivial destructor.
  virtual ~WheelVelocityError() {}

  // setters
  /// \brief Set the measurement.
  /// @param[in] velocity_measurement The measurement.
  void setMeasurement(const Eigen::Vector3d& velocity_measurement)
  {
    measurement_ = velocity_measurement;
  }

  /// \brief Set the information.
  /// @param[in] information The information (weight) matrix.
  void setInformation(const Eigen::Matrix3d& information);

  // Set coordinate for ENU to ECEF conversion
  void setCoordinate(const GeoCoordinatePtr& coordinate) {
    coordinate_ = coordinate;
  }

  // getters
  /// \brief Get the measurement.
  /// @return The measurement vector.
  const Eigen::Vector3d& measurement() const { return measurement_; }

  // error term and Jacobian implementation
  /**
    * @brief This evaluates the error term and additionally computes the Jacobians.
    * @param parameters Pointer to the parameters (see ceres)
    * @param residuals Pointer to the residual vector (see ceres)
    * @param jacobians Pointer to the Jacobians (see ceres)
    * @return success of the evaluation.
    */
  virtual bool Evaluate(double const* const * parameters, double* residuals,
                        double** jacobians) const;

  /**
   * @brief This evaluates the error term and additionally computes
   *        the Jacobians in the minimal internal representation.
   * @param parameters Pointer to the parameters (see ceres)
   * @param residuals Pointer to the residual vector (see ceres)
   * @param jacobians Pointer to the Jacobians (see ceres)
   * @param jacobians_minimal Pointer to the minimal Jacobians (equivalent to jacobians).
   * @return Success of the evaluation.
   */
  bool EvaluateWithMinimalJacobians(double const* const * parameters,
                                    double* residuals, double** jacobians,
                                    double** jacobians_minimal) const;

  // sizes
  /// \brief Residual dimension.
  size_t residualDim() const { return kNumResiduals; }

  /// \brief Number of parameter blocks.
  size_t parameterBlocks() const { return dims_.kNumParameterBlocks; }

  /// \brief Dimension of an individual parameter block.
  size_t parameterBlockDim(size_t parameter_block_idx) const
  {
    return dims_.GetDim(parameter_block_idx);
  }

  /// @brief Residual block type as string
  virtual ErrorType typeInfo() const
  {
    return ErrorType::kWheelVelocityError;
  }

  // Convert normalized residual to raw residual
  virtual void deNormalizeResidual(double *residuals) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 1>> Residual(residuals);
    Residual = square_root_information_inverse_ * Residual;
  }

protected:
  Eigen::Vector3d measurement_; ///< The measurement.
  Eigen::Vector3d angular_velocity_;

  // weighting related
  Eigen::Matrix3d information_; ///< The 3x3 information matrix.
  Eigen::Matrix3d square_root_information_; ///< The 3x3 square root information matrix.
  Eigen::Matrix3d square_root_information_inverse_;
  Eigen::Matrix3d covariance_; ///< The 3x3 covariance matrix.

  // Parameter dimensions
  ceres::internal::StaticParameterDims<Ns...> dims_;

  // Geodetic coordinate
  GeoCoordinatePtr coordinate_;

  // parameter types
  bool is_estimate_body_;
  int parameter_block_group_;
};

// Explicitly instantiate template classes
template class WheelVelocityError<3, 1>;  // Group 1
template class WheelVelocityError<7, 3, 3, 1>;  // Group 2

}  // namespace gici
