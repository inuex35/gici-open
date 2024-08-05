#include "gici/gnss/tdcp_error.h"
#include "gici/gnss/gnss_common.h"
#include "gici/utility/transform.h"
#include "gici/estimate/pose_local_parameterization.h"

namespace gici {

template<int... Ns>
TDCPError<Ns...>::TDCPError(
    const Satellite& last_sat,
    const Satellite& cur_sat,
    const State& last_state,
    const State& cur_state,
    const GnssErrorParameter& error_parameter)
  : last_sat_(last_sat), cur_sat_(cur_sat),
    last_state_(last_state), cur_state_(cur_state), error_parameter_(error_parameter),
    is_estimate_body_(false), parameter_block_group_(0)
{
  if (dims_.kNumParameterBlocks == 2 && 
      dims_.GetDim(2) == 7 && dims_.GetDim(3) == 7) {
    is_estimate_body_ = false;
    parameter_block_group_ = 1;
  }   if (dims_.kNumParameterBlocks == 2 && 
      dims_.GetDim(0) == 7 && dims_.GetDim(1) == 7) {
    is_estimate_body_ = false;
    parameter_block_group_ = 2;
  } 
    else {
    LOG(FATAL) << "TDCPError parameter blocks setup invalid!";
  }

  // Initialize other necessary variables and data structures
  setInformation(error_parameter_);
  for (const auto& [obs_key, obs] : cur_sat_.observations) {
    observation_ = obs; // Use the first observation found in current satellite
    break;
  }
  for (const auto& [obs_key, obs] : last_sat_.observations) {
    observation2_ = obs; // Use the first observation found in last satellite
    break;
  }

}

template<int... Ns>
void TDCPError<Ns ...>::setInformation(const GnssErrorParameter& error_parameter)
{
  error_parameter_ = error_parameter;
  double factor = error_parameter_.tdcp_error_factor;
  covariance_ = covariance_t(factor * factor);
  //char system = satellite_.getSystem();
  //covariance_ *= error_parameter_.system_error_ratio.at(system) * error_parameter_.system_error_ratio.at(system);

  information_ = covariance_.inverse();
  Eigen::LLT<information_t> lltOfInformation(information_);
  square_root_information_ = lltOfInformation.matrixL().transpose();
  square_root_information_inverse_ = square_root_information_.inverse();
}

template<int... Ns>
bool TDCPError<Ns ...>::Evaluate(double const* const * parameters,
                                 double* residuals, double** jacobians) const
{
  return EvaluateWithMinimalJacobians(parameters, residuals, jacobians, nullptr);
}

template<int... Ns>
bool TDCPError<Ns ...>::EvaluateWithMinimalJacobians(
    double const* const * parameters, double* residuals, double** jacobians,
    double** jacobians_minimal) const
{
    Eigen::Vector3d t_WR_ECEF_1, t_WR_ECEF_2, t_WS_W, t_SR_S;
    Eigen::Quaterniond q_WS;
    Eigen::Vector3d v_WR_ECEF_1, v_WR_ECEF_2;

        t_WR_ECEF_1 = Eigen::Map<const Eigen::Vector3d>(parameters[0]);
        t_WR_ECEF_2 = Eigen::Map<const Eigen::Vector3d>(parameters[1]);

        double rho = gnss_common::satelliteToReceiverDistance(satellite_.sat_position, t_WR_ECEF_1);
        Eigen::Vector3d e = (satellite_.sat_position - t_WR_ECEF_1) / rho;
        Eigen::Vector3d v_sat = satellite_.sat_velocity;
        Eigen::Vector3d p_sat = satellite_.sat_position;

        double dist = (t_WR_ECEF_1 - t_WR_ECEF_2).norm();
        double tdcp = abs((observation_.phaserange - observation2_.phaserange));
        double err = tdcp - dist;
        if (abs(err) > 0.1) err = 0.01;
        Eigen::Matrix<double, 1, 1> error = Eigen::Matrix<double, 1, 1>(err);

        LOG(INFO) << "------------------------";
        LOG(INFO) << std::setprecision(15) << "err" << err;
        LOG(INFO) << std::setprecision(15) << "observation_.wavelength:" << observation_.wavelength;
        LOG(INFO) << std::setprecision(15) << "observation_.pseudorange:" << observation_.pseudorange;
        LOG(INFO) << std::setprecision(15) << "observation_.phaserange:" << observation_.phaserange;
        LOG(INFO) << std::setprecision(15) << "observation_.pseudorange2:" << observation2_.pseudorange;
        LOG(INFO) << std::setprecision(15) << "observation_.phaserange2:" << observation2_.phaserange;
        LOG(INFO) << std::setprecision(15) << "t_WR_ECEF_1:" << t_WR_ECEF_1;
        LOG(INFO) << std::setprecision(15) << "t_WR_ECEF_2:" << t_WR_ECEF_2;
        LOG(INFO) << std::setprecision(15) << "tdcp:" << tdcp << " dist:" << dist << " error:" << error;
        LOG(INFO) << "------------------------";


        Eigen::Map<Eigen::Matrix<double, 1, 1>> weighted_error(residuals);
        weighted_error = square_root_information_ * error;

        if (jacobians != nullptr) {
            if (jacobians[0] != nullptr) {
                Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> J0(jacobians[0]);
                J0.setZero();
                J0.block<1, 3>(0, 0) = -e.transpose(); // Set the jacobian for t_WR_ECEF
            }
            if (jacobians[1] != nullptr) {
                Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> J1(jacobians[1]);
                J1.setZero();
                J1.block<1, 3>(0, 0) = - e.transpose(); // Set the jacobian for v_WR_ECEF_1
            }
    }
    return true;
}

}  // namespace gici
