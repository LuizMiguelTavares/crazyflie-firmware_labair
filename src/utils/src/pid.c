/**
 *    ||          ____  _ __
 * +------+      / __ )(_) /_______________ _____  ___
 * | 0xBC |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * +------+    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *  ||  ||    /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2011-2012 Bitcraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * pid.c - implementation of the PID regulator
 */

#include "pid.h"
#include "num.h"
#include <math.h>
#include <float.h>
#include "autoconf.h"
#include "param.h"

static int8_t is_PD_ASMC = 1;
static int8_t use_superTwist = 0;

void pidInit(PidObject* pid, const float desired, const float kp,
             const float ki, const float kd, const float kff, const float dt,
             const float samplingRate, const float cutoffFreq,
             bool enableDFilter)
{
  pid->error         = 0;
  pid->prevMeasured  = 0;
  pid->integ         = 0;
  pid->deriv         = 0;
  pid->desired       = desired;
  pid->kp            = kp;
  pid->ki            = ki;
  pid->kd            = kd;
  pid->kff           = kff;
  pid->iLimit        = DEFAULT_PID_INTEGRATION_LIMIT;
  pid->outputLimit   = DEFAULT_PID_OUTPUT_LIMIT;
  pid->dt            = dt;
  pid->enableDFilter = enableDFilter;
  if (pid->enableDFilter)
  {
    lpf2pInit(&pid->dFilter, samplingRate, cutoffFreq);
  }
}

void smcInit(PidObject* pid, const float k_smc, const float lambda_smc,
             const float sigma_smc, const float ki_smc, const float delta_smc, const float kj_smc, const float sigma_sdelta_smc, const int axis_smc)
{
  pid->smc.k_smc = k_smc;
  pid->smc.lambda_smc = lambda_smc;
  pid->smc.sigma_smc = sigma_smc;
  pid->smc.ki_smc = ki_smc;
  pid->smc.delta_smc = delta_smc;
  pid->smc.kj_smc = kj_smc;
  pid->smc.sigma_sdelta_smc = sigma_sdelta_smc;
  pid->smc.axis = axis_smc;
}

float pidUpdateSMC(PidObject* pid, const float measured, const bool isYawAngle, const setpoint_t *setpoint, const state_t *state)
{
  float output = 0.0f;

  pid->error = pid->desired - measured;
  
  if (isYawAngle){
    if (pid->error > 180.0f){
      pid->error -= 360.0f;
    } else if (pid->error < -180.0f){
      pid->error += 360.0f;
    }
  }
  
  pid->outP = pid->kp * pid->error;
  output += pid->outP;

  /*
  * Note: The derivative term in this PID controller is implemented based on the
  * derivative of the measured process variable instead of the error.
  * This approach avoids derivative kick, which can occur due to sudden changes
  * in the setpoint. By using the process variable for the derivative calculation, we achieve
  * smoother and more stable control during setpoint changes.
  */
  float delta = -(measured - pid->prevMeasured);

  // For yaw measurements, take care of spikes when crossing 180deg <-> -180deg  
  if (isYawAngle){
    if (delta > 180.0f){
      delta -= 360.0f;
    } else if (delta < -180.0f){
      delta += 360.0f;
    }
  }
  
  #if CONFIG_CONTROLLER_PID_FILTER_ALL
    pid->deriv = delta / pid->dt;
  #else
    if (pid->enableDFilter){
      pid->deriv = lpf2pApply(&pid->dFilter, delta / pid->dt);
    } else {
      pid->deriv = delta / pid->dt;
    }
  #endif
  if (isnan(pid->deriv)) {
    pid->deriv = 0;
  }
  pid->outD = pid->kd * pid->deriv;
  output += pid->outD;

  pid->integ += pid->error * pid->dt;

  // Constrain the integral (unless the iLimit is zero)
  if(pid->iLimit != 0)
  {
    pid->integ = constrain(pid->integ, -pid->iLimit, pid->iLimit);
  }

  pid->outI = pid->ki * pid->integ;
  if (!is_PD_ASMC){
    output += pid->outI;
  }

  pid->outFF = pid->kff * pid->desired;
  output += pid->outFF;


  // SMC part
  if (setpoint->thrust > 12000){
    float lambda_smc = pid->smc.lambda_smc;
    float sigma_smc = pid->smc.sigma_smc;

    float ew, er;

    float deg2rad = 0.01745329251;

    ew = pid->error*deg2rad;

    float phi = state->attitude.roll*deg2rad;    // phi
    float theta = -state->attitude.pitch*deg2rad; // theta
    float psi = state->attitude.yaw*deg2rad;     // psi

    float ddxref =  setpoint->acceleration.x;
    float ddyref =  setpoint->acceleration.y;
    float ddzref =  setpoint->acceleration.z;

    float gv = 9.80665f;
    float thrust_raw = (ddzref + gv)/(cosf(theta)*cosf(phi));
    float R13d = ddxref/thrust_raw;
    float R23d = ddyref/thrust_raw;

    float psid =  0.0f; // TODO Ainda não implementado

    float R11 = cosf(theta)*cosf(psi);
    float R12 = sinf(phi)*sinf(theta)*cosf(psi) - cosf(phi)*sinf(psi);
    float R13 = sinf(phi)*sinf(psi) + cosf(phi)*sinf(theta)*cosf(psi);
    float R21 = cosf(theta)*sinf(psi);
    float R22 = cosf(phi)*cosf(psi) + sinf(phi)*sinf(theta)*sinf(psi);
    float R23 = cosf(phi)*sinf(theta)*sinf(psi) - sinf(phi)*cosf(psi);
    float R33 = cosf(phi)*cosf(theta);

    float R13_error = R13d - R13;
    float R23_error = R23d - R23;

    float er_phi = (R13_error*R21 + R23_error*(-R11))/R33;
    float er_theta = (R13_error*R22 + R23_error*(-R12))/R33;
    float er_psi = psid - psi;

    if (pid->smc.axis == 1) {
      er = er_phi;
    } else if (pid->smc.axis == 2) {
      er = er_theta;
    } else if (pid->smc.axis == 3) {
      er = er_psi;
    } else {
      er = 0;
    }

    float s_smc = ew + 0*lambda_smc*er;
    pid->smc.s_smc = s_smc;

    // Saturation
    float s_smc_sat;

    if (s_smc/sigma_smc < -1.0f) {
    s_smc_sat = -1.0f;
    } else if(s_smc/sigma_smc > 1.0f) {
    s_smc_sat = 1.0f;
    } else {
    s_smc_sat = s_smc/sigma_smc;
    }

    // if (s_smc < 0.0f) {
    //   s_smc_sat = -1.0f;
    // } else {
    //   s_smc_sat = 1.0f;
    // }


    if (pid->smc.axis == 3) {
      s_smc_sat = 0.0f;
    }

    float dk_smc, s_minus_delta, s_minus_delta_sat;
    s_minus_delta = fabsf(s_smc) - pid->smc.delta_smc;

    if (s_minus_delta < 0.0f) {
      s_minus_delta_sat = -1.0f;
    } else {
      s_minus_delta_sat = 1.0f;
    }

    // Saturation Sdelta
    float s_delta_smc_sat;

    float sigma_sdelta_smc = pid->smc.sigma_sdelta_smc;

    if (s_smc/(sigma_sdelta_smc) < -1.0f) {
    s_delta_smc_sat = -1.0f;
    } else if(s_smc/(sigma_sdelta_smc) > 1.0f) {
    s_delta_smc_sat = 1.0f;
    } else {
    s_delta_smc_sat = s_smc/(sigma_sdelta_smc);
    }

    float s_delta, s_delta_sat;
    s_delta = s_smc - sigma_sdelta_smc*s_delta_smc_sat;
    // s_delta = s_smc;

    if (s_delta > 0.001f) {
      s_delta_sat = 1.0f;
    } else if (s_delta < -0.001f) {
      s_delta_sat = -1.0f;
    } else {
      s_delta_sat = 0.0f;
    }

    // dk_smc = pid->smc.kj_smc * s_delta_sat + pid->smc.ki_smc * s_minus_delta_sat; // fabsf(s_smc)
    dk_smc = pid->smc.ki_smc * s_minus_delta_sat; // fabsf(s_smc)
    // dk_smc = pid->smc.kj_smc * s_delta_sat; // fabsf(s_smc)
    pid->smc.k_smc += dk_smc * pid->dt;

    if (pid->smc.k_smc < pid->smc.k_min_smc){
      pid->smc.k_smc = pid->smc.k_min_smc;
    } else if (pid->smc.k_smc > pid->smc.k_max_smc){
      pid->smc.k_smc = pid->smc.k_max_smc;
    }

    pid->smc.outPID = output;

    if (~use_superTwist){
      pid->smc.outSMC = pid->smc.k_smc*s_smc_sat;
    }


    // // Super Twist Part
    float s_smc_sgn;

    if (s_smc < 0.0f) {
      s_smc_sgn = -1.0f;
    } else {
      s_smc_sgn = 1.0f;
    }

    if (use_superTwist){
      pid->smc.integ_st_smc += s_smc_sgn * pid->dt;
      // Constrain the integral (unless the iLimit is zero)
      float integral_term_st = pid->smc.kj_smc * pid->smc.integ_st_smc;
      if(pid->smc.i_st_limit_smc != 0)
      {
        integral_term_st = constrain(integral_term_st, -pid->smc.i_st_limit_smc, pid->smc.i_st_limit_smc);
      }
      pid->smc.outSMC = pid->smc.ki_smc * sqrtf(fabsf(s_smc)) * s_smc_sgn + integral_term_st;
    }

    if (is_PD_ASMC){
      output += pid->smc.outSMC;
    }
  }

  #if CONFIG_CONTROLLER_PID_FILTER_ALL
    //filter complete output instead of only D component to compensate for increased noise from increased barometer influence
    if (pid->enableDFilter)
    {
      output = lpf2pApply(&pid->dFilter, output);
    }
    else {
      output = output;
    }
    if (isnan(output)) {
      output = 0;
    }
  #endif

  // Constrain the total PID output (unless the outputLimit is zero)
  if(pid->outputLimit != 0)
  {
    output = constrain(output, -pid->outputLimit, pid->outputLimit);
  }

  pid->prevMeasured = measured;

  return output;
}

float pidUpdate(PidObject* pid, const float measured, const bool isYawAngle)
{
  float output = 0.0f;

  pid->error = pid->desired - measured;
  
  if (isYawAngle){
    if (pid->error > 180.0f){
      pid->error -= 360.0f;
    } else if (pid->error < -180.0f){
      pid->error += 360.0f;
    }
  }
  
  pid->outP = pid->kp * pid->error;
  output += pid->outP;

  pid->smc.outSMC = 0.0f;

  /*
  * Note: The derivative term in this PID controller is implemented based on the
  * derivative of the measured process variable instead of the error.
  * This approach avoids derivative kick, which can occur due to sudden changes
  * in the setpoint. By using the process variable for the derivative calculation, we achieve
  * smoother and more stable control during setpoint changes.
  */
  float delta = -(measured - pid->prevMeasured);

  // For yaw measurements, take care of spikes when crossing 180deg <-> -180deg  
  if (isYawAngle){
    if (delta > 180.0f){
      delta -= 360.0f;
    } else if (delta < -180.0f){
      delta += 360.0f;
    }
  }
  
  #if CONFIG_CONTROLLER_PID_FILTER_ALL
    pid->deriv = delta / pid->dt;
  #else
    if (pid->enableDFilter){
      pid->deriv = lpf2pApply(&pid->dFilter, delta / pid->dt);
    } else {
      pid->deriv = delta / pid->dt;
    }
  #endif
  if (isnan(pid->deriv)) {
    pid->deriv = 0;
  }
  pid->outD = pid->kd * pid->deriv;
  output += pid->outD;

  pid->integ += pid->error * pid->dt;

  // Constrain the integral (unless the iLimit is zero)
  if(pid->iLimit != 0)
  {
    pid->integ = constrain(pid->integ, -pid->iLimit, pid->iLimit);
  }

  pid->outI = pid->ki * pid->integ;
  output += pid->outI;

  pid->outFF = pid->kff * pid->desired;
  output += pid->outFF;
  
  #if CONFIG_CONTROLLER_PID_FILTER_ALL
    //filter complete output instead of only D component to compensate for increased noise from increased barometer influence
    if (pid->enableDFilter)
    {
      output = lpf2pApply(&pid->dFilter, output);
    }
    else {
      output = output;
    }
    if (isnan(output)) {
      output = 0;
    }
  #endif

  // Constrain the total PID output (unless the outputLimit is zero)
  if(pid->outputLimit != 0)
  {
    output = constrain(output, -pid->outputLimit, pid->outputLimit);
  }

  pid->prevMeasured = measured;

  return output;
}

void pidSetIntegralLimit(PidObject* pid, const float limit) {
    pid->iLimit = limit;
}

void pidReset(PidObject* pid, const float actual)
{
  pid->error        = 0;
  pid->prevMeasured = actual;
  pid->integ        = 0;
  pid->deriv        = 0;
}

void pidSetDesired(PidObject* pid, const float desired)
{
  pid->desired = desired;
}

float pidGetDesired(PidObject* pid)
{
  return pid->desired;
}

bool pidIsActive(PidObject* pid)
{
  bool isActive = true;

  if (pid->kp < 0.0001f && pid->ki < 0.0001f && pid->kd < 0.0001f)
  {
    isActive = false;
  }

  return isActive;
}

void pidSetKp(PidObject* pid, const float kp)
{
  pid->kp = kp;
}

void pidSetKi(PidObject* pid, const float ki)
{
  pid->ki = ki;
}

void pidSetKd(PidObject* pid, const float kd)
{
  pid->kd = kd;
}

void pidSetKff(PidObject* pid, const float kff)
{
  pid->kff = kff;
}

void pidSetDt(PidObject* pid, const float dt) {
    pid->dt = dt;
}

void filterReset(PidObject* pid, const float samplingRate, const float cutoffFreq, bool enableDFilter) {
  pid->enableDFilter = enableDFilter;
  if (pid->enableDFilter)
  {
    lpf2pInit(&pid->dFilter, samplingRate, cutoffFreq);
  }
}


PARAM_GROUP_START(smc)

/**
 * @brief Boudary layer for the SMC PID pitch rate controller
 */
PARAM_ADD(PARAM_UINT8, is_PD_ASMC, &is_PD_ASMC)

/**
 * @brief Boudary layer for the SMC PID pitch rate controller
 */
PARAM_ADD(PARAM_UINT8, use_superTwist, &use_superTwist)
PARAM_GROUP_STOP(smc)

