/**
* @Function: Wheel Velocity error block for ceres backend
*
* @Author  : [Your Name]
* @Email   : [Your Email]
*
* Copyright (C) 2023 by [Your Name], All rights reserved.
**/
#include "gici/gnss/wheel_velocity_error.h"

#include "gici/utility/transform.h"
#include "gici/gnss/gnss_common.h"

namespace gici {

// This evaluates the error term and additionally computes the Jacobians.
template<int... Ns>
bool WheelVelocityError<Ns ...>::Evaluate(double const* const * parameters,
                                          double* residuals, double** jacobians) const
{
  return EvaluateWithMinimalJacobians(parameters, residuals, jacobians, nullptr);
}

// This evaluates the error term and additionally computes the Jacobians
// in the minimal internal representation.
template<int... Ns>
bool WheelVelocityError<Ns ...>::EvaluateWithMinimalJacobians(
    double const* const * parameters, double* residuals, double** jacobians,
    double** jacobians_minimal) const
{
  Eigen::Vector3d t_WR_ECEF, v_WR, angular_velocity;
  
  // Get the current state
  t_WR_ECEF = Eigen::Map<const Eigen::Vector3d>(parameters[0]);
  v_WR = Eigen::Map<const Eigen::Vector3d>(parameters[1]);
  angular_velocity = Eigen::Map<const Eigen::Vector3d>(parameters[2]);

  // Compute velocity in ECEF frame
  Eigen::Vector3d v_WR_ECEF = coordinate_->rotate(v_WR, GeoType::ENU, GeoType::ECEF);

  // Compute error
  Eigen::Vector3d error = measurement_ - v_WR_ECEF;

  // weigh it
  Eigen::Map<Eigen::Vector3d> weighted_error(residuals);
  weighted_error = square_root_information_ * error;

  // compute Jacobians
  if (jacobians != nullptr)
  {
    // Velocity in ENU
    Eigen::Matrix<double, 3, 3> J_v_ENU = coordinate_->rotationMatrix(GeoType::ENU, GeoType::ECEF);

    // Velocity in ECEF
    Eigen::Matrix<double, 3, 3> J_v_ECEF = Eigen::Matrix<double, 3, 3>::Identity();

    // Jacobians w.r.t. ENU velocity
    if (jacobians[1] != nullptr) {
      Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> J1(jacobians[1]);
      J1 = square_root_information_ * J_v_ENU;

      if (jacobians_minimal != nullptr && jacobians_minimal[1] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> J1_minimal(jacobians_minimal[1]);
        J1_minimal = J1;
      }
    }
  }

  return true;
}

}  // namespace gici

// Explicitly instantiate template classes
template class gici::WheelVelocityError<3, 1>;

