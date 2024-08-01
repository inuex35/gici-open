#include "gici/gnss/tdcp_error.h"
#include "gici/gnss/gnss_common.h"
#include "gici/utility/transform.h"
#include "gici/estimate/pose_local_parameterization.h"

namespace gici {

template<int... Ns>
TDCPError<Ns ...>::TDCPError(
                       const GnssMeasurement& last_measurement,
                       const GnssMeasurement& cur_measurement,
                       const State& last_state, const State& cur_state, const GnssErrorParameter& error_parameter)
  : last_measurement_(last_measurement), cur_measurement_(cur_measurement),
    last_state_(last_state), cur_state_(cur_state), error_parameter_(error_parameter)
{
  if (dims_.kNumParameterBlocks == 4 && 
      dims_.GetDim(0) == 3 && dims_.GetDim(1) == 3 && 
      dims_.GetDim(2) == 1 && dims_.GetDim(3) == 1) {
    is_estimate_body_ = false;
    parameter_block_group_ = 1;
  } else if (dims_.kNumParameterBlocks == 5 &&
      dims_.GetDim(0) == 7 && dims_.GetDim(1) == 3 &&
      dims_.GetDim(2) == 3 && dims_.GetDim(3) == 1 &&
      dims_.GetDim(4) == 1) {
    is_estimate_body_ = true;
    parameter_block_group_ = 2;
  } else {
    LOG(FATAL) << "TDCPError parameter blocks setup invalid!";
  }

  setInformation(error_parameter);
}

template<int... Ns>
void TDCPError<Ns ...>::setInformation(const GnssErrorParameter& error_parameter)
{
  error_parameter_ = error_parameter;
  double factor = error_parameter_.tdcp_error_factor;
  covariance_ = covariance_t(factor * factor);
  char system = satellite_.getSystem();
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
  Eigen::Vector3d t_WR_ECEF, t_WS_W, t_SR_S;
  Eigen::Quaterniond q_WS;
  Eigen::Vector3d v_WR_ECEF_1, v_WR_ECEF_2;
  double clock_bias_1, clock_bias_2;
  
  if (!is_estimate_body_) 
  {
    t_WR_ECEF = Eigen::Map<const Eigen::Vector3d>(parameters[0]);
    v_WR_ECEF_1 = Eigen::Map<const Eigen::Vector3d>(parameters[1]);
    clock_bias_1 = parameters[2][0];
    clock_bias_2 = parameters[3][0];
  }
  else 
  {
    t_WS_W = Eigen::Map<const Eigen::Vector3d>(&parameters[0][0]);
    q_WS = Eigen::Map<const Eigen::Quaterniond>(&parameters[0][3]);

    v_WR_ECEF_1 = Eigen::Map<const Eigen::Vector3d>(&parameters[1][0]);
    v_WR_ECEF_2 = Eigen::Map<const Eigen::Vector3d>(&parameters[2][0]);

    t_SR_S = Eigen::Map<const Eigen::Vector3d>(parameters[3]);

    clock_bias_1 = parameters[4][0];
    clock_bias_2 = parameters[5][0];

    Eigen::Vector3d t_WR_W = t_WS_W + q_WS * t_SR_S;

    Eigen::Vector3d v_WR_1 = v_WR_ECEF_1 + skewSymmetric(angular_velocity_) * q_WS * t_SR_S;
    Eigen::Vector3d v_WR_2 = v_WR_ECEF_2 + skewSymmetric(angular_velocity_) * q_WS * t_SR_S;

    if (!coordinate_) {
      LOG(FATAL) << "Coordinate not set!";
    }
    if (!coordinate_->isZeroSetted()) {
      LOG(FATAL) << "Coordinate zero not set!";
    }
    t_WR_ECEF = coordinate_->convert(t_WR_W, GeoType::ENU, GeoType::ECEF);
    v_WR_ECEF_1 = coordinate_->rotate(v_WR_1, GeoType::ENU, GeoType::ECEF);
    v_WR_ECEF_2 = coordinate_->rotate(v_WR_2, GeoType::ENU, GeoType::ECEF);
  }

  double rho = gnss_common::satelliteToReceiverDistance(satellite_.sat_position, t_WR_ECEF);
  Eigen::Vector3d e = (satellite_.sat_position - t_WR_ECEF) / rho;
  Eigen::Vector3d v_sat = satellite_.sat_velocity;
  Eigen::Vector3d p_sat = satellite_.sat_position;
  Eigen::Vector3d vs_1 = v_sat - v_WR_ECEF_1;
  Eigen::Vector3d vs_2 = v_sat - v_WR_ECEF_2;

  double range_rate_1 = vs_1.dot(e) + OMGE / CLIGHT *
      (v_sat(1) * t_WR_ECEF(0) + p_sat(1) * v_WR_ECEF_1(0) -
       v_sat(0) * t_WR_ECEF(1) - p_sat(0) * v_WR_ECEF_1(1));
  double range_rate_2 = vs_2.dot(e) + OMGE / CLIGHT *
      (v_sat(1) * t_WR_ECEF(0) + p_sat(1) * v_WR_ECEF_2(0) -
       v_sat(0) * t_WR_ECEF(1) - p_sat(0) * v_WR_ECEF_2(1));
  double tdcp_estimate = 
    (range_rate_2 - range_rate_1) + (clock_bias_2 - clock_bias_1);

  double tdcp = observation_.phaserange - observation2_.phaserange;
  Eigen::Matrix<double, 1, 1> error = 
    Eigen::Matrix<double, 1, 1>(tdcp - tdcp_estimate);

  Eigen::Map<Eigen::Matrix<double, 1, 1> > weighted_error(residuals);
  weighted_error = square_root_information_ * error;

  if (jacobians != nullptr)
  {
    Eigen::Matrix<double, 1, 3> J_t_ECEF = Eigen::Matrix<double, 1, 3>::Zero();
    Eigen::Matrix<double, 1, 3> J_v_ECEF_1 = -e.transpose();
    Eigen::Matrix<double, 1, 3> J_v_ECEF_2 = e.transpose();

    Eigen::Matrix<double, 1, 6> J_T_WS;
    Eigen::Matrix<double, 1, 9> J_speed_and_bias;
    Eigen::Matrix<double, 1, 3> J_t_SR_S;

    if (is_estimate_body_) {
      Eigen::Matrix<double, 1, 3> J_t_W = Eigen::Matrix<double, 1, 3>::Zero();

      Eigen::Matrix3d R_ECEF_ENU = coordinate_->rotationMatrix(GeoType::ENU, GeoType::ECEF);
      Eigen::Matrix<double, 1, 3> J_v_W_1 = J_v_ECEF_1 * R_ECEF_ENU;
      Eigen::Matrix<double, 1, 3> J_v_W_2 = J_v_ECEF_2 * R_ECEF_ENU;

      Eigen::Matrix<double, 1, 3> J_q_WS_1 = J_v_W_1 * skewSymmetric(angular_velocity_) * -skewSymmetric(q_WS.toRotationMatrix() * t_SR_S);
      Eigen::Matrix<double, 1, 3> J_q_WS_2 = J_v_W_2 * skewSymmetric(angular_velocity_) * -skewSymmetric(q_WS.toRotationMatrix() * t_SR_S);

      J_T_WS.setZero();
      J_T_WS.topLeftCorner(1, 3) = J_t_W;
      J_T_WS.topRightCorner(1, 3) = J_q_WS_1 + J_q_WS_2;

      J_speed_and_bias.setZero();
      J_speed_and_bias.topLeftCorner(1, 3) = J_v_W_1;
      J_speed_and_bias.topRightCorner(1, 3) = J_v_W_2;

      J_t_SR_S = J_v_W_1 * skewSymmetric(angular_velocity_) * q_WS.toRotationMatrix();
    }

    Eigen::Matrix<double, 1, 1> J_bias_1 = -Eigen::MatrixXd::Identity(1, 1);
    Eigen::Matrix<double, 1, 1> J_bias_2 = Eigen::MatrixXd::Identity(1, 1);

    if (parameter_block_group_ == 1) 
    {
      if (jacobians[0] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J0(jacobians[0]);
        J0 = square_root_information_ * J_t_ECEF;

        if (jacobians_minimal != nullptr && jacobians_minimal[0] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J0_minimal_mapped(jacobians_minimal[0]);
          J0_minimal_mapped = J0;
        }
      }
      if (jacobians[1] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J1(jacobians[1]);
        J1 = square_root_information_ * J_v_ECEF_1;

        if (jacobians_minimal != nullptr && jacobians_minimal[1] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J1_minimal_mapped(jacobians_minimal[1]);
          J1_minimal_mapped = J1;
        }
      }
      if (jacobians[2] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> J2(jacobians[2]);
        J2 = square_root_information_ * J_bias_1;

        if (jacobians_minimal != nullptr && jacobians_minimal[2] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> J2_minimal_mapped(jacobians_minimal[2]);
          J2_minimal_mapped = J2;
        }
      }
      if (jacobians[3] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> J3(jacobians[3]);
        J3 = square_root_information_ * J_bias_2;

        if (jacobians_minimal != nullptr && jacobians_minimal[3] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> J3_minimal_mapped(jacobians_minimal[3]);
          J3_minimal_mapped = J3;
        }
      }
    }
    if (parameter_block_group_ == 2)
    {
      if (jacobians[0] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> J0(jacobians[0]);
        Eigen::Matrix<double, 1, 6, Eigen::RowMajor> J0_minimal;
        J0_minimal = square_root_information_ * J_T_WS;

        Eigen::Matrix<double, 6, 7, Eigen::RowMajor> J_lift;
        PoseLocalParameterization::liftJacobian(parameters[0], J_lift.data());

        J0 = J0_minimal * J_lift;

        if (jacobians_minimal != nullptr && jacobians_minimal[0] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 6, Eigen::RowMajor>> J0_minimal_mapped(jacobians_minimal[0]);
          J0_minimal_mapped = J0_minimal;
        }
      }
      if (jacobians[1] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J1(jacobians[1]);
        J1 = square_root_information_ * J_v_ECEF_1;

        if (jacobians_minimal != nullptr && jacobians_minimal[1] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J1_minimal_mapped(jacobians_minimal[1]);
          J1_minimal_mapped = J1;
        }
      }
      if (jacobians[2] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J2(jacobians[2]);
        J2 = square_root_information_ * J_v_ECEF_2;

        if (jacobians_minimal != nullptr && jacobians_minimal[2] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J2_minimal_mapped(jacobians_minimal[2]);
          J2_minimal_mapped = J2;
        }
      }
      if (jacobians[3] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J3(jacobians[3]);
        J3 = square_root_information_ * J_t_SR_S;

        if (jacobians_minimal != nullptr && jacobians_minimal[3] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> J3_minimal_mapped(jacobians_minimal[3]);
          J3_minimal_mapped = J3;
        }
      }
      if (jacobians[4] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> J4(jacobians[4]);
        J4 = square_root_information_ * J_bias_1;

        if (jacobians_minimal != nullptr && jacobians_minimal[4] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> J4_minimal_mapped(jacobians_minimal[4]);
          J4_minimal_mapped = J4;
        }
      }
      if (jacobians[5] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> J5(jacobians[5]);
        J5 = square_root_information_ * J_bias_2;

        if (jacobians_minimal != nullptr && jacobians_minimal[5] != nullptr) {
          Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> J5_minimal_mapped(jacobians_minimal[5]);
          J5_minimal_mapped = J5;
        }
      }
    }
  }

  return true;
}

}  // namespace gici
