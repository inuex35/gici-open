/**
* @Function: Align differential measurements
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#pragma once

#include <deque>

#include "gici/gnss/gnss_types.h"
#include "gici/estimate/estimator_types.h"
#include "gici/utility/common.h"

namespace gici {

// Estimator
class DifferentialMeasurementsAlign {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // The default constructor
  DifferentialMeasurementsAlign() {}
  ~DifferentialMeasurementsAlign() {}

  // Set measurements
  inline void add(const EstimatorDataCluster& measurement) {
    if (measurement.gnss && measurement.gnss_role == GnssRole::Rover) {
      measurement_rov_.push_back(*measurement.gnss);
    }
    if (measurement.gnss && measurement.gnss_role == GnssRole::Heading) {
      measurement_heading_.push_back(*measurement.gnss);
    }
    if (measurement.gnss && measurement.gnss_role == GnssRole::Reference) {
      measurement_ref_.push_back(*measurement.gnss);
    }
  }

  // Get aligned
  inline bool get(const double max_age, 
    GnssMeasurement& rov, GnssMeasurement& ref) {

    const double offset = 0.0;  // shift a time offset for test

    if (measurement_rov_.size() == 0) return false;
    if (measurement_ref_.size() == 0) {
      measurement_rov_.clear(); return false;
    }

    rov = measurement_rov_.front();
    measurement_rov_.pop_front();

    // get the nearest timestamp
    double min_dt = 1.0e6;
    int index = -1;
    for (int i = measurement_ref_.size() - 1; i >= 0; i--) {
      double dt = fabs(rov.timestamp - measurement_ref_.at(i).timestamp - offset);
      if (min_dt > dt) {
        min_dt = dt; index = i;
      }
    }
    CHECK(index != -1);
    ref = measurement_ref_.at(index);
    // pop all in front of this
    double cut_timestamp = ref.timestamp;
    while (!checkEqual(cut_timestamp, measurement_ref_.front().timestamp)) {
      measurement_ref_.pop_front();
    }

    // Check timestamps
    if (!checkEqual(rov.timestamp, ref.timestamp, max_age)) {
      LOG(WARNING) << "Max age between two measurements exceeded! "
        << "age = " << fabs(rov.timestamp - ref.timestamp)
        << ", max_age = " << max_age;
      return false;
    }

    return true;
  }

inline bool multi_ant_get(const double max_age, 
                const double max_age_rov_heading,
                GnssMeasurement& rov, 
                GnssMeasurement& heading, 
                GnssMeasurement& ref) {
    const double offset = 0.0;  // shift a time offset for test

    LOG(INFO) << "measurement_rov_.size():" << measurement_rov_.size();
    LOG(INFO) << "measurement_heading_.size():" << measurement_heading_.size();
    LOG(INFO) << "measurement_ref_.size():" << measurement_ref_.size();

    if (measurement_rov_.size() == 0) return false;
    if (measurement_heading_.size() == 0) return false;
    if (measurement_ref_.size() == 0) {
        measurement_rov_.clear(); 
        measurement_heading_.clear();
        return false;
    }

    rov = measurement_rov_.front();
    measurement_rov_.pop_front();

    // Get the nearest timestamp for ref
    double min_dt_ref = 1.0e6;
    int index_ref = -1;
    for (int i = measurement_ref_.size() - 1; i >= 0; i--) {
        double dt = fabs(rov.timestamp - measurement_ref_.at(i).timestamp - offset);
        if (min_dt_ref > dt) {
            min_dt_ref = dt; index_ref = i;
        }
    }
    CHECK(index_ref != -1);
    ref = measurement_ref_.at(index_ref);

    // Get the nearest timestamp for heading
    double min_dt_heading = 1.0e6;
    int index_heading = -1;
    for (int i = measurement_heading_.size() - 1; i >= 0; i--) {
        double dt = fabs(rov.timestamp - measurement_heading_.at(i).timestamp - offset);
        if (min_dt_heading > dt) {
            min_dt_heading = dt; index_heading = i;
        }
    }
    CHECK(index_heading != -1);
    heading = measurement_heading_.at(index_heading);

    // Pop all in front of the nearest ref timestamp
    double cut_timestamp_ref = ref.timestamp;
    while (!checkEqual(cut_timestamp_ref, measurement_ref_.front().timestamp)) {
        measurement_ref_.pop_front();
    }

    // Pop all in front of the nearest heading timestamp
    double cut_timestamp_heading = heading.timestamp;
    while (!checkEqual(cut_timestamp_heading, measurement_heading_.front().timestamp)) {
        measurement_heading_.pop_front();
    }

    // Check timestamps for rov and ref
    if (!checkEqual(rov.timestamp, ref.timestamp, max_age)) {
        LOG(WARNING) << "Max age between rov and ref exceeded! "
                     << "age = " << fabs(rov.timestamp - ref.timestamp)
                     << ", max_age = " << max_age;
        return false;
    }

    // Check timestamps for rov and heading
    if (!checkEqual(rov.timestamp, heading.timestamp, max_age_rov_heading)) {
        LOG(WARNING) << "Max age between rov and heading exceeded! "
                     << "age = " << fabs(rov.timestamp - heading.timestamp)
                     << ", max_age_rov_heading = " << max_age_rov_heading;
        return false;
    }

    return true;
}

protected:
  // Measurement storage for aligning
  std::deque<GnssMeasurement> measurement_rov_;
  std::deque<GnssMeasurement> measurement_heading_;
  std::deque<GnssMeasurement> measurement_ref_;
};

}