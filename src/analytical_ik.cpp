#include <math.h>
#include <stdio.h>
#include "analytical_ik.h"

namespace ur_kinematics {

  namespace {
    const double ZERO_THRESH = 0.00000001;

    int SIGN(double x) {
      return (x > 0) - (x < 0);
    }

    const double PI = M_PI;

    const double d1 =  0.089159;
    const double a2 = -0.42500;
    const double a3 = -0.39225;
    const double d4 =  0.10915;
    const double d5 =  0.09465;
    const double d6 =  0.0823;
  }

  void forward(const double* q, double* T) {
    double s1 = sin(*q), c1 = cos(*q); q++;
    double q23 = *q, q234 = *q, s2 = sin(*q), c2 = cos(*q); q++;
    double s3 = sin(*q), c3 = cos(*q); q23 += *q; q234 += *q; q++;
    double s4 = sin(*q), c4 = cos(*q); q234 += *q; q++;
    double s5 = sin(*q), c5 = cos(*q); q++;
    double s6 = sin(*q), c6 = cos(*q); 
    double s23 = sin(q23), c23 = cos(q23);
    double s234 = sin(q234), c234 = cos(q234);
    *T = c234*c1*s5 - c5*s1; T++;
    *T = c6*(s1*s5 + c234*c1*c5) - s234*c1*s6; T++;
    *T = -s6*(s1*s5 + c234*c1*c5) - s234*c1*c6; T++;
    *T = d6*c234*c1*s5 - a3*c23*c1 - a2*c1*c2 - d6*c5*s1 - d5*s234*c1 - d4*s1; T++;
    *T = c1*c5 + c234*s1*s5; T++;
    *T = -c6*(c1*s5 - c234*c5*s1) - s234*s1*s6; T++;
    *T = s6*(c1*s5 - c234*c5*s1) - s234*c6*s1; T++;
    *T = d6*(c1*c5 + c234*s1*s5) + d4*c1 - a3*c23*s1 - a2*c2*s1 - d5*s234*s1; T++;
    *T = -s234*s5; T++;
    *T = -c234*s6 - s234*c5*c6; T++;
    *T = s234*c5*s6 - c234*c6; T++;
    *T = d1 + a3*s23 + a2*s2 - d5*(c23*c4 - s23*s4) - d6*s5*(c23*s4 + s23*c4); T++;
    *T = 0.0; T++; *T = 0.0; T++; *T = 0.0; T++; *T = 1.0;
  }

  int inverse(const double* T, double* q_sols,  double q6_des=0.0) {
    int num_sols = 0;
    double T02 = -*T; T++; double T00 =  *T; T++; double T01 =  *T; T++; double T03 = -*T; T++; 
    double T12 = -*T; T++; double T10 =  *T; T++; double T11 =  *T; T++; double T13 = -*T; T++; 
    double T22 =  *T; T++; double T20 = -*T; T++; double T21 = -*T; T++; double T23 =  *T;

    // shoulder rotate joint (q1)
    double q1[2];
    {
      double A = d6*T12 - T13;
      double B = d6*T02 - T03;
      double R = A*A + B*B;

      if (fabs(A) < ZERO_THRESH) {
        double div;
        if (fabs(fabs(d4) - fabs(B)) < ZERO_THRESH)
          div = -SIGN(d4) * SIGN(B);
        else
          div = -d4 / B;

        double arcsin = asin(div);
        if (fabs(arcsin) < ZERO_THRESH)
          arcsin = 0.0;

        if (arcsin < 0.0)
          q1[0] = arcsin + 2.0 * PI;
        else
          q1[0] = arcsin;

        q1[1] = PI - arcsin;
      }
      else if (fabs(B) < ZERO_THRESH) {
        double div;
        if (fabs(fabs(d4) - fabs(A)) < ZERO_THRESH)
          div = SIGN(d4) * SIGN(A);
        else
          div = d4 / A;

        double arccos = acos(div);
        q1[0] = arccos;
        q1[1] = 2.0 * PI - arccos;
      }
      else if (d4 * d4 > R) {
        return num_sols;
      }
      else {
        double arccos = acos(d4 / sqrt(R)) ;
        double arctan = atan2(-B, A);
        double pos = arccos + arctan; 
        double neg = -arccos + arctan; 

        if (fabs(pos) < ZERO_THRESH)
          pos = 0.0;
        if (fabs(neg) < ZERO_THRESH)
          neg = 0.0;

        if (pos >= 2.0 * PI)
          q1[0] = pos - 2.0 * PI;
        else if (pos >= 0.0)
          q1[0] = pos;
        else
          q1[0] = 2.0 * PI + pos;

        if (neg >= 2.0 * PI)
          q1[1] = neg - 2.0 * PI;
        else if (neg >= 0.0)
          q1[1] = neg; 
        else
          q1[1] = 2.0 * PI + neg;
      }
    }

    // wrist 2 joint (q5) 
    double q5[2][2];
    {
      for (int i = 0; i < 2; i++) {
        double numer = (T03*sin(q1[i]) - T13*cos(q1[i]) - d4);
        double div;
        if (fabs(fabs(numer) - fabs(d6)) < ZERO_THRESH)
          div = SIGN(numer) * SIGN(d6);
        else 
          div = numer / d6;

        double arccos = acos(div);
        q5[i][0] = arccos;
        q5[i][1] = 2.0 * PI - arccos;
      }
    }

    {
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          double c1 = cos(q1[i]), s1 = sin(q1[i]);
          double c5 = cos(q5[i][j]), s5 = sin(q5[i][j]);
          double q6;

          // wrist 3 joint (q6) 
          if (fabs(s5) < ZERO_THRESH){
            q6 = q6_des;
          }
          else {
            q6 = atan2(SIGN(s5) * -(T01*s1 - T11*c1), 
                       SIGN(s5) * (T00*s1 - T10*c1));

            if (fabs(q6) < ZERO_THRESH)
              q6 = 0.0;

            if (q6 < 0.0)
              q6 += 2.0 * PI;
          }

          double q2[2], q3[2], q4[2];

          // RRR joints (q2,q3,q4) 
          double c6 = cos(q6), s6 = sin(q6);
          double x04x = -s5*(T02*c1 + T12*s1) - c5*(s6*(T01*c1 + T11*s1) - c6*(T00*c1 + T10*s1));
          double x04y = c5*(T20*c6 - T21*s6) - T22*s5;
          double p13x = d5*(s6*(T00*c1 + T10*s1) + c6*(T01*c1 + T11*s1)) - d6*(T02*c1 + T12*s1) + T03*c1 + T13*s1;
          double p13y = T23 - d1 - d6*T22 + d5*(T21*c6 + T20*s6);
          double c3 = (p13x*p13x + p13y*p13y - a2*a2 - a3*a3) / (2.0*a2*a3);

          if (fabs(fabs(c3) - 1.0) < ZERO_THRESH)
            c3 = SIGN(c3);
          else if (fabs(c3) > 1.0)
            continue;
          
          double arccos = acos(c3);
          q3[0] = arccos;
          q3[1] = 2.0 * PI - arccos;

          double denom = a2*a2 + a3*a3 + 2*a2*a3*c3;
          double s3 = sin(arccos);
          double A = (a2 + a3*c3), B = a3*s3;
          q2[0] = atan2((A*p13y - B*p13x) / denom, (A*p13x + B*p13y) / denom);
          q2[1] = atan2((A*p13y + B*p13x) / denom, (A*p13x - B*p13y) / denom);

          double c23_0 = cos(q2[0] + q3[0]);
          double s23_0 = sin(q2[0] + q3[0]);
          double c23_1 = cos(q2[1] + q3[1]);
          double s23_1 = sin(q2[1] + q3[1]);
          q4[0] = atan2(c23_0*x04y - s23_0*x04x, x04x*c23_0 + x04y*s23_0);
          q4[1] = atan2(c23_1*x04y - s23_1*x04x, x04x*c23_1 + x04y*s23_1);

          // solutions
          for (int k = 0; k < 2; k++) {
            if (fabs(q2[k]) < ZERO_THRESH)
              q2[k] = 0.0;
            else if (q2[k] < 0.0) 
              q2[k] += 2.0 * PI;

            if (fabs(q4[k]) < ZERO_THRESH)
              q4[k] = 0.0;
            else if (q4[k] < 0.0) 
              q4[k] += 2.0 * PI;

            q_sols[num_sols * 6 + 0] = q1[i];    q_sols[num_sols * 6 + 1] = q2[k]; 
            q_sols[num_sols * 6 + 2] = q3[k];    q_sols[num_sols * 6 + 3] = q4[k]; 
            q_sols[num_sols * 6 + 4] = q5[i][j]; q_sols[num_sols * 6 + 5] = q6; 
            num_sols++;
          }
        }
      }
    }
    return num_sols;
  }

  void jacobian(double *jacobian, double *q) {
    const double cos1 = cos(q[0]);
    const double cos2 = cos(q[1]);
    const double cos3 = cos(q[2]);
    const double cos4 = cos(q[3]);
    const double cos5 = cos(q[4]);

    const double sin1 = sin(q[0]);
    const double sin2 = sin(q[1]);
    const double sin3 = sin(q[2]);
    const double sin4 = sin(q[3]);
    const double sin5 = sin(q[4]);

    const double sin234 = sin(q[1] + q[2] + q[3]);
    const double cos234 = cos(q[1] + q[2] + q[3]);

    // jacobian[0 * 6 + 0] = -(d5 * (cos1 * cos234 + sin1 * sin234)) / 2 + (d5 * (cos1 * cos234 - sin1 * sin234)) / 2 + d4 * cos1 - (d6 * (-sin1 * cos234 - cos1 * sin234) * sin5) / 2 - (d6 * (-sin1 * cos234 + cos1 * cos234) * sin5) / 2 - (a2 * sin1 * cos2) + (d6 * cos5 * cos1) - (a3 * sin1 * cos2 * cos3) + (a3 * sin1 * sin2 * sin3);
    // jacobian[0 * 6 + 1] = (d5 * (sin1 * sin234 + cos1 * cos234)) / 2 + (d5 * (-sin1 * sin234 + cos1 * cos234)) / 2 + d4 * sin1 - (d6 * (-cos1 * sin234 - sin1 * cos234) * sin5) / 2 - (d6 * (-cos1 * sin234 + sin1 * cos234) * sin5) / 2 - (a2 * cos1 * sin2) + (d6 * cos5 * sin1) - (a3 * cos1 * sin2 * cos3) - (a3 * cos1 * cos2 * sin3);
    // jacobian[0 * 6 + 2] = (d5 * (sin1 * sin234 + cos1 * cos234)) / 2 + (d5 * (-sin1 * sin234 + cos1 * cos234)) / 2 + d4 * sin1 - (d6 * (-cos1 * sin234 - sin1 * cos234) * sin5) / 2 - (d6 * (-cos1 * sin234 + sin1 * cos234) * sin5) / 2 + (a2 * cos1 * cos2) + (d6 * cos5 * sin1) - (a3 * cos1 * cos2 * sin3) - (a3 * cos1 * sin2 * cos3);
    // jacobian[0 * 6 + 3] = (d5 * (sin1 * sin234 + cos1 * cos234)) / 2 + (d5 * (-sin1 * sin234 + cos1 * cos234)) / 2 + d4 * sin1 - (d6 * (-cos1 * sin234 - sin1 * cos234) * sin5) / 2 - (d6 * (-cos1 * sin234 + sin1 * cos234) * sin5) / 2 + (a2 * cos1 * cos2) + (d6 * cos5 * sin1) + (a3 * cos1 * cos2 * cos3) - (a3 * cos1 * sin2 * sin3);
    // jacobian[0 * 6 + 4] = -(d5 * (sin1 * cos234 - cos1 * sin234)) / 2 + (d5 * (sin1 * cos234 + cos1 * sin234)) / 2 + d4 * sin1 - (d6 * (cos1 * cos234 - sin1 * sin234) * cos5) / 2 - (d6 * (cos1 * cos234 + sin1 * sin234) * cos5) / 2 + (a2 * cos1 * cos2) - (d6 * sin5 * sin1) + (a3 * cos1 * cos2 * cos3) - (a3 * cos1 * sin2 * sin3);
    // jacobian[0 * 6 + 5] = 0;


    // jacobian[1 * 6 + 0] = (d5 * (sin1 * cos234 + cos1 * sin234)) / 2 + (d5 * (-sin1 * cos234 + cos1 * sin234)) / 2 - d4 * sin1 - (d6 * ( cos1 * cos234 + sin1 * sin234) * sin5) / 2 - (d6 * ( cos1 * cos234 - sin1 * sin234) * sin5) / 2 + (d6 * sin1 * cos5) + (a2 * cos2 * cos1) + (a3 * cos2 * cos3 * cos1) - (a3 * cos1 * sin2 * sin3);
    // jacobian[1 * 6 + 1] = (d5 * (cos1 * sin234 + sin1 * cos234)) / 2 + (d5 * (-cos1 * sin234 + sin1 * cos234)) / 2 + d4 * cos1 - (d6 * (-sin1 * sin234 - cos1 * cos234) * sin5) / 2 - (d6 * (-sin1 * sin234 + cos1 * cos234) * sin5) / 2 - (d6 * cos1 * cos5) - (a2 * sin2 * sin1) - (a3 * sin2 * cos3 * sin1) - (a3 * sin1 * cos2 * sin3);
    // jacobian[1 * 6 + 2] = (d5 * (cos1 * sin234 + sin1 * cos234)) / 2 + (d5 * (-cos1 * sin234 + sin1 * cos234)) / 2 + d4 * cos1 - (d6 * (-sin1 * sin234 - cos1 * cos234) * sin5) / 2 - (d6 * (-sin1 * sin234 + cos1 * cos234) * sin5) / 2 - (d6 * cos1 * cos5) + (a2 * cos2 * sin1) - (a3 * cos2 * sin3 * sin1) - (a3 * sin1 * sin2 * cos3); 
    // jacobian[1 * 6 + 3] = (d5 * (cos1 * sin234 + sin1 * cos234)) / 2 + (d5 * (-cos1 * sin234 + sin1 * cos234)) / 2 + d4 * cos1 - (d6 * (-sin1 * sin234 - cos1 * cos234) * sin5) / 2 - (d6 * (-sin1 * sin234 + cos1 * cos234) * sin5) / 2 - (d6 * cos1 * cos5) + (a2 * cos2 * sin1) + (a3 * cos2 * cos3 * sin1) - (a3 * sin1 * sin2 * sin3);
    // jacobian[1 * 6 + 4] = -(d5 * (cos1 * cos234 - sin1 * sin234)) / 2 + (d5 * (cos1 * cos234 + sin1 * sin234)) / 2 + d4 * cos1 - (d6 * ( sin1 * cos234 - cos1 * sin234) * cos5) / 2 - (d6 * ( sin1 * cos234 + cos1 * sin234) * cos5) / 2 + (d6 * cos1 * sin5) + (a2 * cos2 * sin1) + (a3 * cos2 * cos3 * sin1) - (a3 * sin1 * sin2 * sin3);
    // jacobian[1 * 6 + 5] = 0;

    jacobian[0 * 6 + 0] = -(d5 * (cos1 * cos234 + sin1 * sin234)) / 2 + (d5 * ( cos1 * cos234 - sin1 * sin234)) / 2 + d4 * cos1 - (d6 * (-sin1 * cos234 - cos1 * sin234) * sin5) / 2 - (d6 * ( -sin1 * cos234 + cos1 * sin234) * sin5) / 2 - (a2 * sin1 * cos2) + (d6 * cos5 * cos1) - (a3 * sin1 * cos2 * cos3) + (a3 * sin1 * sin2 * sin3);
    jacobian[0 * 6 + 1] = (d5 * (sin1 * sin234 + cos1 * cos234)) / 2 + (d5 * (-sin1 * sin234 + cos1 * cos234)) / 2 - (d6 * (-cos1 * sin234 - sin1 * cos234) * sin5) / 2 - (d6 * ( -cos1 * sin234 + sin1 * cos234) * sin5) / 2 - (a2 * cos1 * sin2) - (a3 * cos1 * sin2 * cos3) - (a3 * cos1 * cos2 * sin3);
    jacobian[0 * 6 + 2] = (d5 * (sin1 * sin234 + cos1 * cos234)) / 2 + (d5 * (-sin1 * sin234 + cos1 * cos234)) / 2 - (d6 * (-cos1 * sin234 - sin1 * cos234) * sin5) / 2 - (d6 * ( -cos1 * sin234 + sin1 * cos234) * sin5) / 2 - (a3 * cos1 * cos2 * sin3) - (a3 * cos1 * sin2 * cos3);
    jacobian[0 * 6 + 3] = (d5 * (sin1 * sin234 + cos1 * cos234)) / 2 + (d5 * (-sin1 * sin234 + cos1 * cos234)) / 2 - (d6 * (-cos1 * sin234 - sin1 * cos234) * sin5) / 2 - (d6 * ( -cos1 * sin234 + sin1 * cos234) * sin5) / 2;
    jacobian[0 * 6 + 4] = -(d6 * ( cos1 * cos234 - sin1 * sin234) * cos5) / 2 - (d6 * ( cos1 * cos234 + sin1 * sin234) * cos5) / 2 - (d6 * sin5 * sin1);
    jacobian[0 * 6 + 5] = 0;


    jacobian[1 * 6 + 0] = (d5 * (sin1 * cos234 + cos1 * sin234)) / 2 + (d5 * (-sin1 * cos234 + cos1 * sin234)) / 2 - d4 * sin1 - (d6 * ( cos1 * cos234 + sin1 * sin234) * sin5) / 2 - (d6 * ( cos1 * cos234 - sin1 * sin234) * sin5) / 2 + (d6 * sin1 * cos5) + (a2 * cos2 * cos1) + (a3 * cos2 * cos3 * cos1) - (a3 * cos1 * sin2 * sin3);
    jacobian[1 * 6 + 1] = (d5 * (cos1 * sin234 + sin1 * cos234)) / 2 + (d5 * (-cos1 * sin234 + sin1 * cos234)) / 2 - (d6 * (-sin1 * sin234 - cos1 * cos234) * sin5) / 2 - (d6 * (-sin1 * sin234 + cos1 * cos234) * sin5) / 2 - (a2 * sin2 * sin1) - (a3 * sin2 * cos3 * sin1) - (a3 * sin1 * cos2 * sin3);
    jacobian[1 * 6 + 2] = (d5 * (cos1 * sin234 + sin1 * cos234)) / 2 + (d5 * (-cos1 * sin234 + sin1 * cos234)) / 2 - (d6 * (-sin1 * sin234 - cos1 * cos234) * sin5) / 2 - (d6 * (-sin1 * sin234 + cos1 * cos234) * sin5) / 2 - (a3 * cos2 * sin3 * sin1) - (a3 * sin1 * sin2 * cos3);
    jacobian[1 * 6 + 3] = (d5 * (cos1 * sin234 + sin1 * cos234)) / 2 + (d5 * (-cos1 * sin234 + sin1 * cos234)) / 2 - (d6 * (-sin1 * sin234 - cos1 * cos234) * sin5) / 2 - (d6 * (-sin1 * sin234 + cos1 * cos234) * sin5) / 2;
    jacobian[1 * 6 + 4] = -(d6 * ( sin1 * cos234 - cos1 * sin234) * cos5) / 2 - (d6 * ( sin1 * cos234 + cos1 * sin234) * cos5) / 2 + (d6 * cos1 * sin5);
    jacobian[1 * 6 + 5] = 0;

    jacobian[2 * 6 + 0] = 0;
    jacobian[2 * 6 + 1] = (d6 * (sin234 * (-cos5) - sin5 * cos234)) / 2 + a2 * cos2 + a3 * (cos2 * cos3 - sin2 * sin3) - (d6 * (sin5 * cos234 - sin234 * cos5)) / 2 + d6 * sin234;
    jacobian[2 * 6 + 2] = (d6 * (sin234 * (-cos5) - sin5 * cos234)) / 2 + a3 * (cos2 * cos3 - sin2 * sin3) - (d6 * (sin5 * cos234 - sin234 * cos5)) / 2 + d6 * sin234;
    jacobian[2 * 6 + 3] = (d6 * (sin234 * (-cos5) - sin5 * cos234)) / 2 - (d6 * (sin5 * cos234 - sin234 * cos5)) / 2 + d6 * sin234;
    jacobian[2 * 6 + 4] = (d6 * (sin234 * (-cos5) - sin5 * cos234)) / 2 - (d6 * (sin234 * cos5 - sin5 * cos234)) / 2 ;
    jacobian[2 * 6 + 5] = 6;
  }
};

std::tuple<double, double, double> forward_kinematics(double *q) {
  double *T = new double[16];
  ur_kinematics::forward(q, T);

  return {T[0*4 + 3], T[1*4 + 3], T[2*4 + 3]};
}

int inverse_kinematics_2PI(double *q_sols, double x, double y, double z) {
  double q[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double *T = new double[16];
  ur_kinematics::forward(q, T);

  std::tuple<double, double, double> target = {x, y, z};
  T[0*4 + 3] = std::get<0>(target);
  T[1*4 + 3] = std::get<1>(target);
  T[2*4 + 3] = std::get<2>(target);

  int num_sols;
  num_sols = ur_kinematics::inverse(T, q_sols);

  for (int i = 0; i < num_sols; i++)
    q_sols[i * 6] -= ur_kinematics::PI / 2; 

  return num_sols;
}

int inverse_kinematics(double *q_sols, double x, double y, double z) {
  double q[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double *T = new double[16];
  ur_kinematics::forward(q, T);

  std::tuple<double, double, double> target = {x, y, z};
  T[0*4 + 3] = std::get<0>(target);
  T[1*4 + 3] = std::get<1>(target);
  T[2*4 + 3] = std::get<2>(target);

  int num_sols;
  num_sols = ur_kinematics::inverse(T, q_sols);

  for (int i = 0; i < num_sols; i++)
    q_sols[i * 6] -= ur_kinematics::PI / 2; 

  for (int i = 0; i < 6 * num_sols; ++i) {
    if (q_sols[i] >= ur_kinematics::PI)
        q_sols[i] -= 2.0 * ur_kinematics::PI;
  }

  return num_sols;
}

#include <iostream>

int joint_jacobian(double *jacobian, double *q) {
  // double q_delta[6] = {q[0], q[1], q[2], q[3], q[4], q[5]};
  // double delta = 0.0001; 

  // // double jacobian[6 * 3];
  // for (int i = 0; i < 6; i++) {
  //   double q_delta[6] = {q[0], q[1], q[2], q[3], q[4], q[5]};
  //   q_delta[i] = q[i] + delta;
  //   auto [q1, q2, q3] = forward_kinematics(q);
  //   auto [q1_delta, q2_delta, q3_delta] = forward_kinematics(q_delta);

  //   jacobian[i * 3] = (q1_delta - q1) / delta;
  //   jacobian[i * 3 + 1] = (q2_delta - q2) / delta;
  //   jacobian[i * 3 + 2] = (q3_delta - q3) / delta;
  // }

  // for (int i = 0; i < 6; i++) {
  //   std::cout << jacobian[i * 3] << " ";
  //   std::cout << jacobian[i * 3 +1] << " ";
  //   std::cout << jacobian[i * 3 + 2] << " ";

  //   std::cout << std::endl;

  // }

  // px = -(d5 * (sin1 * cos234 - cos1 * sin234)) / 2 + (d5 * (sin1 * cos234 + cos1 * sin234)) / 2 + d4 * sin1 - (d6 * ( cos1 * cos234 - sin1 * sin234) * sin5) / 2 - (d6 * ( cos1 * cos234 + sin1 * sin234) * sin5) / 2 + (a2 * cos1 * cos2) + (d6 * cos5 * sin1) + (a3 * cos1 * cos2 * cos3) - (a3 * cos1 * sin2 * sin3)

  // py = -(d5 * (cos1 * cos234 - sin1 * sin234)) / 2 + (d5 * (cos1 * cos234 + sin1 * sin234)) / 2 + d4 * cos1 - (d6 * ( sin1 * cos234 - cos1 * sin234) * sin5) / 2 - (d6 * ( sin1 * cos234 + cos1 * sin234) * sin5) / 2 - (d6 * cos1 * cos5) + (a2 * cos2 * sin1) + (a3 * cos2 * cos3 * sin1) - (a3 * sin1 * sin2 * sin3)

  // pz = d_1 + (d6 * (cos234 * cos5 - sin234 * sin5)) / 2 + (a3 * (sin2 * cos3 + cos2 * sin3)) + (a2 * sin2) - (d6 * (cos234 * cos5 + sin234 * sin5)) / 2 - (d5 * cos234)
  q[0] -= ur_kinematics::PI / 2;

  ur_kinematics::jacobian(jacobian, q);

  return 0;
}