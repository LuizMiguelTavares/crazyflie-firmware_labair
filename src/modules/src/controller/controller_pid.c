
#include "stabilizer_types.h"

#include "attitude_controller.h"
#include "position_controller.h"
#include "controller_pid.h"

#include "log.h"
#include "param.h"
#include "math3d.h"

#include "usec_time.h"

#define ATTITUDE_UPDATE_DT    (float)(1.0f/ATTITUDE_RATE)

static attitude_t attitudeDesired;
static attitude_t rateDesired;
static float actuatorThrust;

static float cmd_thrust;
static float cmd_roll;
static float cmd_pitch;
static float cmd_yaw;
static float r_roll;
static float r_pitch;
static float r_yaw;
static float accelz;

static float control_dt_us = 0.0f;

static uint8_t useSMC = 0;
static float k_phi = 50.0f;
static float k_theta = 50.0f;
static float k_psi = 12.5f;

static float sat_roll = 5.0f;
static float sat_pitch = 5.0f;
static float sat_yaw = 3.0f;

static float er_phi = 0.0f;
static float er_theta = 0.0f;
static float er_psi = 0.0f;

#define MA_N 10
static float rb_phi[MA_N]   = {0};
static float rb_theta[MA_N] = {0};
static float rb_psi[MA_N]   = {0};
static uint8_t rb_k = 0;    
static uint8_t rb_fill = 0;

static float phi_filt;
static float theta_filt;
static float psi_filt;

static float psi_prev = 0.0f;
static float psi_unwrapped = 0.0f;
static uint8_t psi_init = 0;

static float R13_error_ = 0;
static float R23_error_ = 0;

#ifndef PI_
#define PI_ 3.14159265358979323846f
#endif

static float deg2rad = 0.01745329251f;
static float rad2deg = 57.2958f;
static float gv = 9.80665f;

void controllerPidInit(void)
{
  attitudeControllerInit(ATTITUDE_UPDATE_DT);
  positionControllerInit();
}

bool controllerPidTest(void)
{
  bool pass = true;

  pass &= attitudeControllerTest();

  return pass;
}

static float capAngle(float angle) {
  float result = angle;

  while (result > 180.0f) {
    result -= 360.0f;
  }

  while (result < -180.0f) {
    result += 360.0f;
  }

  return result;
}

void controllerPid(control_t *control, const setpoint_t *setpoint,
                                         const sensorData_t *sensors,
                                         const state_t *state,
                                         const stabilizerStep_t stabilizerStep)
{
  control->controlMode = controlModeLegacy;
  rb_phi[rb_k]   = state->attitude.roll * deg2rad;
  rb_theta[rb_k] = -state->attitude.pitch * deg2rad;

  float psi_wrapped = state->attitude.yaw * deg2rad;

  if (!psi_init) {
    psi_init = 1;
    psi_prev = psi_wrapped;
    psi_unwrapped = psi_wrapped;
  } else {
    float d = psi_wrapped - psi_prev;

    if (d >  PI_) d -= 2.0f * PI_;
    if (d < -PI_) d += 2.0f * PI_;

    psi_unwrapped += d;
    psi_prev = psi_wrapped;
  }

  rb_psi[rb_k] = psi_unwrapped;

  rb_k++;
  if (rb_k >= MA_N) rb_k = 0;
  if (rb_fill < MA_N) rb_fill++;

  if (RATE_DO_EXECUTE(ATTITUDE_RATE, stabilizerStep)) {
    // Rate-controled YAW is moving YAW angle setpoint
    if (setpoint->mode.yaw == modeVelocity) {
      attitudeDesired.yaw = capAngle(attitudeDesired.yaw + setpoint->attitudeRate.yaw * ATTITUDE_UPDATE_DT);

      float yawMaxDelta = attitudeControllerGetYawMaxDelta();
      if (yawMaxDelta != 0.0f)
      {
      float delta = capAngle(attitudeDesired.yaw-state->attitude.yaw);
      // keep the yaw setpoint within +/- yawMaxDelta from the current yaw
        if (delta > yawMaxDelta)
        {
          attitudeDesired.yaw = state->attitude.yaw + yawMaxDelta;
        }
        else if (delta < -yawMaxDelta)
        {
          attitudeDesired.yaw = state->attitude.yaw - yawMaxDelta;
        }
      }
    } else if (setpoint->mode.yaw == modeAbs) {
      attitudeDesired.yaw = setpoint->attitude.yaw;
    } else if (setpoint->mode.quat == modeAbs) {
      struct quat setpoint_quat = mkquat(setpoint->attitudeQuaternion.x, setpoint->attitudeQuaternion.y, setpoint->attitudeQuaternion.z, setpoint->attitudeQuaternion.w);
      struct vec rpy = quat2rpy(setpoint_quat);
      attitudeDesired.yaw = degrees(rpy.z);
    }

    attitudeDesired.yaw = capAngle(attitudeDesired.yaw);
  }

  if (RATE_DO_EXECUTE(LABAIR_RATE, stabilizerStep)) {
    if (useSMC==1) {
      float ddxref =  setpoint->acceleration.x;
      float ddyref =  setpoint->acceleration.y;
      float ddzref =  setpoint->acceleration.z;

      float sum_phi = 0.0f, sum_theta = 0.0f, sum_psi = 0.0f;
      int n_ang = rb_fill;
      
      for (int i = 0; i < n_ang; i++) {
        sum_phi   += rb_phi[i];
        sum_theta += rb_theta[i];
        sum_psi   += rb_psi[i];
      }

      rb_fill = 0;
      rb_k = 0;

      phi_filt   = sum_phi / (float)n_ang;
      theta_filt = sum_theta / (float)n_ang;
      psi_filt   = sum_psi / (float)n_ang;

      while (psi_filt >  PI_) psi_filt -= 2.0f * PI_;
      while (psi_filt < -PI_) psi_filt += 2.0f * PI_;

      psi_init = 0;
      psi_prev = 0.0f;
      psi_unwrapped = 0.0f;

      // float phi   = state->attitude.roll * deg2rad;
      // float theta = -state->attitude.pitch * deg2rad;
      // float psi   = state->attitude.yaw * deg2rad;

      float phi   = phi_filt;
      float theta = theta_filt;
      float psi   = psi_filt;

      float thrust_raw = (ddzref + gv)/(cosf(theta)*cosf(phi));
      float R13d = ddxref/thrust_raw;
      float R23d = ddyref/thrust_raw;

      float R11 = cosf(theta)*cosf(psi);
      float R12 = sinf(phi)*sinf(theta)*cosf(psi) - cosf(phi)*sinf(psi);
      float R13 = sinf(phi)*sinf(psi) + cosf(phi)*sinf(theta)*cosf(psi);
      float R21 = cosf(theta)*sinf(psi);
      float R22 = cosf(phi)*cosf(psi) + sinf(phi)*sinf(theta)*sinf(psi);
      float R23 = cosf(phi)*sinf(theta)*sinf(psi) - sinf(phi)*cosf(psi);
      float R33 = cosf(phi)*cosf(theta);

      float R13_error = k_phi*(R13d - R13);
      float R23_error = k_theta*(R23d - R23);
      
      R13_error_ = R13_error;
      R23_error_ = R23_error;

      er_phi = (R13_error*R21 + R23_error*(-R11))/R33;
      er_theta = (R13_error*R22 + R23_error*(-R12))/R33;
      // er_psi = -sinf(phi)*sensors->gyro.y*deg2rad + R33*(k_psi*(0 - psi)); // Psi desejado está em 0!

      if (er_phi > sat_roll) {
        er_phi = sat_roll;
      } else if (er_phi < -sat_roll) {
        er_phi = -sat_roll;
      }

      if (er_theta > sat_pitch) {
        er_theta = sat_pitch;
      } else if (er_theta < -sat_pitch) {
        er_theta = -sat_pitch;
      }
     
      // if (er_psi > sat_yaw) {
      //   er_psi = sat_yaw;
      // } else if (er_psi < -sat_yaw) {
      //   er_psi = -sat_yaw;
      // }

      er_phi = er_phi*rad2deg;
      er_theta = er_theta*rad2deg;
      // er_psi = er_psi*rad2deg;
    }
  }

  if (RATE_DO_EXECUTE(POSITION_RATE, stabilizerStep)) {
    positionController(&actuatorThrust, &attitudeDesired, setpoint, state);
  }

  if (RATE_DO_EXECUTE(ATTITUDE_RATE, stabilizerStep)) {
    uint32_t t0 = usecTimestamp();

    // Switch between manual and automatic position control

    if (setpoint->mode.z == modeDisable) {
      actuatorThrust = setpoint->thrust;
    }
    if (setpoint->mode.x == modeDisable || setpoint->mode.y == modeDisable) {
      attitudeDesired.roll = setpoint->attitude.roll;
      attitudeDesired.pitch = setpoint->attitude.pitch;
    }

    attitudeControllerCorrectAttitudePID(state->attitude.roll, state->attitude.pitch, state->attitude.yaw,
                                attitudeDesired.roll, attitudeDesired.pitch, attitudeDesired.yaw,
                                &rateDesired.roll, &rateDesired.pitch, &rateDesired.yaw);

    // For roll and pitch, if velocity mode, overwrite rateDesired with the setpoint
    // value. Also reset the PID to avoid error buildup, which can lead to unstable
    // behavior if level mode is engaged later
    if (setpoint->mode.roll == modeVelocity) {
      rateDesired.roll = setpoint->attitudeRate.roll;
      attitudeControllerResetRollAttitudePID(state->attitude.roll);
    }
    if (setpoint->mode.pitch == modeVelocity) {
      rateDesired.pitch = setpoint->attitudeRate.pitch;
      attitudeControllerResetPitchAttitudePID(state->attitude.pitch);
    }

    // TODO: Investigate possibility to subtract gyro drift.
    if (useSMC==0) {
      attitudeControllerCorrectRatePID(sensors->gyro.x, -sensors->gyro.y, sensors->gyro.z,
                            rateDesired.roll, rateDesired.pitch, rateDesired.yaw);
    } else {
      attitudeControllerCorrectRatePIDSMC(sensors->gyro.x, -sensors->gyro.y, sensors->gyro.z,
                              er_phi, -er_theta, rateDesired.yaw, setpoint, state);
    }

    attitudeControllerGetActuatorOutput(&control->roll,
                                        &control->pitch,
                                        &control->yaw);

    control->yaw = -control->yaw;

    cmd_thrust = control->thrust;
    cmd_roll = control->roll;
    cmd_pitch = control->pitch;
    cmd_yaw = control->yaw;
    r_roll = radians(sensors->gyro.x);
    r_pitch = -radians(sensors->gyro.y);
    r_yaw = radians(sensors->gyro.z);
    accelz = sensors->acc.z;

    uint32_t t1 = usecTimestamp();
    const uint32_t dt_us = t1 - t0;

    control_dt_us = (float)dt_us;  
  }

  control->thrust = actuatorThrust;

  if (control->thrust == 0)
  {
    control->thrust = 0;
    control->roll = 0;
    control->pitch = 0;
    control->yaw = 0;

    cmd_thrust = control->thrust;
    cmd_roll = control->roll;
    cmd_pitch = control->pitch;
    cmd_yaw = control->yaw;

    attitudeControllerResetAllPID(state->attitude.roll, state->attitude.pitch, state->attitude.yaw);
    positionControllerResetAllPID(state->position.x, state->position.y, state->position.z);

    // Reset the calculated YAW angle for rate control
    attitudeDesired.yaw = state->attitude.yaw;
  }
}

/**
 * Logging variables for the command and reference signals for the
 * altitude PID controller
 */
LOG_GROUP_START(controller)
/**
 * @brief Thrust command
 */
LOG_ADD(LOG_FLOAT, cmd_thrust, &cmd_thrust)
/**
 * @brief Roll command
 */
LOG_ADD(LOG_FLOAT, cmd_roll, &cmd_roll)
/**
 * @brief Pitch command
 */
LOG_ADD(LOG_FLOAT, cmd_pitch, &cmd_pitch)
/**
 * @brief yaw command
 */
LOG_ADD(LOG_FLOAT, cmd_yaw, &cmd_yaw)
/**
 * @brief Gyro roll measurement in radians
 */
LOG_ADD(LOG_FLOAT, r_roll, &r_roll)
/**
 * @brief Gyro pitch measurement in radians
 */
LOG_ADD(LOG_FLOAT, r_pitch, &r_pitch)
/**
 * @brief Yaw  measurement in radians
 */
LOG_ADD(LOG_FLOAT, r_yaw, &r_yaw)
/**
 * @brief Acceleration in the zaxis in G-force
 */
LOG_ADD(LOG_FLOAT, accelz, &accelz)
/**
 * @brief Thrust command without (tilt)compensation
 */
LOG_ADD(LOG_FLOAT, actuatorThrust, &actuatorThrust)
/**
 * @brief Desired roll setpoint
 */
LOG_ADD(LOG_FLOAT, roll,      &attitudeDesired.roll)
/**
 * @brief Desired pitch setpoint
 */
LOG_ADD(LOG_FLOAT, pitch,     &attitudeDesired.pitch)
/**
 * @brief Desired yaw setpoint
 */
LOG_ADD(LOG_FLOAT, yaw,       &attitudeDesired.yaw)
/**
 * @brief Desired roll rate setpoint
 */
LOG_ADD(LOG_FLOAT, rollRate,  &rateDesired.roll)
/**
 * @brief Desired pitch rate setpoint
 */
LOG_ADD(LOG_FLOAT, pitchRate, &rateDesired.pitch)
/**
 * @brief Desired yaw rate setpoint
 */
LOG_ADD(LOG_FLOAT, yawRate,   &rateDesired.yaw)
LOG_GROUP_STOP(controller)

LOG_GROUP_START(timeControl)
LOG_ADD(LOG_FLOAT, control_us, &control_dt_us)
LOG_GROUP_STOP(timeControl)

LOG_GROUP_START(smc) //remove
LOG_ADD(LOG_FLOAT, er_phi, &er_phi)
LOG_ADD(LOG_FLOAT, er_theta, &er_theta)
LOG_ADD(LOG_FLOAT, er_psi, &er_psi)

LOG_ADD(LOG_FLOAT, phi_filt, &phi_filt)
LOG_ADD(LOG_FLOAT, theta_filt, &theta_filt)
LOG_ADD(LOG_FLOAT, psi_filt, &psi_filt)

LOG_ADD(LOG_FLOAT, R13_error, &R13_error_)
LOG_ADD(LOG_FLOAT, R23_error, &R23_error_)
LOG_GROUP_STOP(smc)

/**
 * Tuning settings for the gains of the PID controller for the rate angles of
 * the Crazyflie, which consists of the yaw, pitch and roll rates 
 */
PARAM_GROUP_START(smc)

/**
 * @brief Boudary layer for the SMC PID pitch rate controller
 */
PARAM_ADD(PARAM_UINT8, useSMC, &useSMC)

PARAM_ADD(PARAM_FLOAT, k_phi, &k_phi)
PARAM_ADD(PARAM_FLOAT, k_theta, &k_theta)
PARAM_ADD(PARAM_FLOAT, k_psi, &k_psi)

PARAM_ADD(PARAM_FLOAT, sat_roll, &sat_roll)
PARAM_ADD(PARAM_FLOAT, sat_pitch, &sat_pitch)
PARAM_ADD(PARAM_FLOAT, sat_yaw, &sat_yaw)

PARAM_GROUP_STOP(smc)
