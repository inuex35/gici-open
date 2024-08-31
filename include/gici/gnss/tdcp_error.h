#pragma once

#include <ceres/ceres.h>
#include <Eigen/Core>

#include "gici/gnss/gnss_common.h"
#include "gici/estimate/error_interface.h"
#include "gici/gnss/geodetic_coordinate.h"
#include "gici/gnss/gnss_types.h"

namespace gici {

template<int... Ns>
class TDCPError :
    public ceres::SizedCostFunction<
    1 /* number of residuals */,
    Ns ... /* parameter blocks */>,
    public ErrorInterface
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef ceres::SizedCostFunction<1, Ns ...> base_t;

  static const int kNumResiduals = 1;

  typedef Eigen::Matrix<double, 1, 1> information_t;
  typedef Eigen::Matrix<double, 1, 1> covariance_t;

  TDCPError() = default;

  TDCPError(double phase1, double phase2, double delta_time)
      : phase1_(phase1), phase2_(phase2), delta_time_(delta_time) {}

  TDCPError(
            const GnssMeasurement& last_measurement,
            const GnssMeasurement& cur_measurement,
            const State& last_state, const State& cur_state,
            const GnssErrorParameter& error_parameter);


  virtual ~TDCPError() {}

  virtual bool Evaluate(double const* const * parameters, double* residuals,
                        double** jacobians) const override;

  bool EvaluateWithMinimalJacobians(double const* const * parameters,
                                    double* residuals, double** jacobians,
                                    double** jacobians_minimal) const;

  void setInformation(const GnssErrorParameter& error_parameter);
  void setCoordinate(const GeoCoordinatePtr& coordinate) {
    coordinate_ = coordinate;
  }
  size_t residualDim() const { return kNumResiduals; }
  size_t parameterBlocks() const { return dims_.kNumParameterBlocks; }
  size_t parameterBlockDim(size_t parameter_block_idx) const
  {
    return dims_.GetDim(parameter_block_idx);
  }

  virtual ErrorType typeInfo() const
  {
    return ErrorType::kTDCPError;
  }

  virtual void deNormalizeResidual(double *residuals) const override
  {
    Eigen::Map<Eigen::Matrix<double, 1, 1>> Residual(residuals);
    Residual = square_root_information_inverse_ * Residual;
  }

protected:
  double phase1_;
  double phase2_;
  double delta_time_;

  information_t information_;
  information_t square_root_information_;
  information_t square_root_information_inverse_;

  ceres::internal::StaticParameterDims<Ns...> dims_;

  // Missing member variables
  GnssMeasurement last_measurement_;
  GnssMeasurement cur_measurement_;
  State last_state_;
  State cur_state_;

  Satellite satellite_;
  Observation observation_;
  Satellite satellite2_;
  Observation observation2_;

  bool is_estimate_body_;
  int parameter_block_group_;
  GnssErrorParameter error_parameter_;
  Eigen::Vector3d angular_velocity_;
  GeoCoordinatePtr coordinate_;
  covariance_t covariance_;
};

template class TDCPError<7, 7, 1, 1>;

}  // namespace gici
