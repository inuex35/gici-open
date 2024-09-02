/**
* @Function: Double-differenced phaserange residual block for ceres backend
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

// Double-differenced phaserange error
// The candidate parameter setups are:
// Group 1: P1. receiver position in ECEF (3), P2. ambiguity (1), 
//          P3. ambiguity of base satellite (1)
// Group 2: P1. body pose in ENU (7), P2. relative position from body to receiver
//          in body frame (3), P3. ambiguity (1), P4. ambiguity of base satellite (1)
// Group 3: Group 1 + P4. troposphere delay at rov (1), P5. troposphere delay at ref (1) 
//          P6. ionosphere delay (1), P7. ionosphere delay of base satellite (1)
// Group 4: Group 2 + P5. troposphere delay at rov (1), P6. troposphere delay at ref (1), 
//          P7. ionosphere delay (1), P8. ionosphere delay of base satellite (1)
template<int... Ns>
class PhaserangeErrorDD :
    public ceres::SizedCostFunction<
    1 /* number of residuals */,
    Ns ... /* parameter blocks */>,
    public ErrorInterface
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// \brief The base class type.
  typedef ceres::SizedCostFunction<1, Ns ...> base_t;

  /// \brief Number of residuals (1).
  static const int kNumResiduals = 1;

  /// \brief The information matrix type (1x1).
  typedef Eigen::Matrix<double, 1, 1> information_t;

  /// \brief The covariance matrix type (same as information).
  typedef Eigen::Matrix<double, 1, 1> covariance_t;

  /// \brief Default constructor.
  PhaserangeErrorDD();

  /// \brief Construct with measurement and information matrix
  /// @param[in] measurement_rov The measurement at rover.
  /// @param[in] measurement_ref The measurement at reference.
  /// @param[in] measurement_rov_base The measurement of base satellite at rover.
  /// @param[in] measurement_ref_base The measurement of base satellite at reference.
  /// @param[in] error_parameter To compute GNSS information matrix.
  PhaserangeErrorDD(const GnssMeasurement& measurement_rov,
                    const GnssMeasurement& measurement_ref,
                    const GnssMeasurementIndex index_rov,
                    const GnssMeasurementIndex index_ref,
                    const GnssMeasurementIndex index_rov_base,
                    const GnssMeasurementIndex index_ref_base,
                    const GnssErrorParameter& error_parameter,
                    bool compass);

  /// \brief Trivial destructor.
  virtual ~PhaserangeErrorDD() {}

  // setters
  /// \brief Set the measurement.
  /// @param[in] measurement The measurement.
  void setMeasurement(const GnssMeasurement& measurement_rov,
                      const GnssMeasurement& measurement_ref)
  {
    measurement_rov_ = measurement_rov;
    measurement_ref_ = measurement_ref;
  }

  /// \brief Set the information.
  /// @param[in] information The information (weight) matrix.
  void setInformation(const GnssErrorParameter& error_parameter);

  // Set coordinate for ENU to ECEF convertion
  void setCoordinate(const GeoCoordinatePtr& coordinate) {
    coordinate_ = coordinate;
  }

  // error term and Jacobian implementation
  /**
    * @brief This evaluates the error term and additionally computes the Jacobians.
    * @param parameters Pointer to the parameters (see ceres)
    * @param residuals Pointer to the residual vector (see ceres)
    * @param jacobians Pointer to the Jacobians (see ceres)
    * @return success of th evaluation.
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
    return ErrorType::kPhaserangeErrorDD;
  }

  // Convert normalized residual to raw residual
  virtual void deNormalizeResidual(double *residuals) const
  {
    Eigen::Map<Eigen::Matrix<double, 1, 1>> Residual(residuals);
    Residual = square_root_information_inverse_ * Residual;
  }

  // Get GNSS index
  inline GnssMeasurementDDIndexPair getGnssMeasurementIndex() { 
    return GnssMeasurementDDIndexPair(
      GnssMeasurementIndex(satellite_rov_.prn, observation_rov_.raw_code),
      GnssMeasurementIndex(satellite_ref_.prn, observation_ref_.raw_code),
      GnssMeasurementIndex(satellite_rov_base_.prn, observation_rov_base_.raw_code),
      GnssMeasurementIndex(satellite_ref_base_.prn, observation_ref_base_.raw_code));
  }

protected:
  GnssMeasurement measurement_rov_, measurement_ref_; 
  Satellite satellite_rov_, satellite_ref_,
            satellite_rov_base_, satellite_ref_base_;
  Observation observation_rov_, observation_ref_,
              observation_rov_base_, observation_ref_base_;

  // weighting related
  GnssErrorParameter error_parameter_;
  information_t information_; ///< The DimxDim information matrix.
  information_t square_root_information_; ///< The DimxDim square root information matrix.
  information_t square_root_information_inverse_;
  covariance_t covariance_; ///< The DimxDim covariance matrix.

  // Parameter dimensions
  ceres::internal::StaticParameterDims<Ns...> dims_;

  // Geodetic coordinate
  GeoCoordinatePtr coordinate_;

  // parameter types
  bool is_estimate_body_;
  bool is_estimate_atmosphere_;
  int parameter_block_group_;
};

// Explicitly instantiate template classes
template class PhaserangeErrorDD<3, 1, 1>;  // Group 1
template class PhaserangeErrorDD<7, 3, 1, 1>;  // Group 2
template class PhaserangeErrorDD<3, 1, 1, 1, 1, 1, 1>;  // Group 3
template class PhaserangeErrorDD<7, 3, 1, 1, 1, 1, 1, 1>;  // Group 4

}  

