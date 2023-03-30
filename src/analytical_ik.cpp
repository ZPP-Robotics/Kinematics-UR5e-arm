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

    const double d1 = 0.16250;
    const double a2 = -0.42500;
    const double a3 = -0.3922;
    const double d4 = 0.1333;
    const double d5 = 0.0997;
    const double d6 = 0.0996;
  }

    void forward_old(const double* q, double* T) {
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

  void forward(const double* q, double* T) {
    double s1 = sin(*q), c1 = cos(*q); q++;
    double q23 = *q, q234 = *q, s2 = sin(*q), c2 = cos(*q); q++;
    double s3 = sin(*q), c3 = cos(*q); q23 += *q; q234 += *q; q++;
    double s4 = sin(*q), c4 = cos(*q); q234 += *q; q++;
    double s5 = sin(*q), c5 = cos(*q); q++;
    double s6 = sin(*q), c6 = cos(*q); 
    double s23 = sin(q23), c23 = cos(q23);
    double s234 = sin(q234), c234 = cos(q234);

    *T = (s1*s2 + c1*c5*c234)*c6 - s6*s234*c1; T++;
    *T = -(s1*s5 + c1*c5*c234)*s6 - s234*c1*c6; T++;
    *T = s1*s5 - s5*c1*c234; T++;
    *T = a2*c1*c2 + a3*c1*c23 + d4*s1 + d5*s234*c1 + d6*s1*c5 - d6*s5*c1*c234; T++;

    *T = (s1*c5*c234 - s5*c1) * c6 - s1*s6*s234; T++;
    *T = (-s1*c5*c234 + s5*c1)*c6 - s1*s234*c6; T++;
    *T = -s1*s5*c234 - c1*c5; T++;
    *T = a2*s1*c2 + a3*s1*c23 - d4*c1 + d5*s1*s234 - d6*s1*s5*c234 - d6*c1*c5; T++;

    *T = s6*c234 + s234*c5*c6; T++;
    *T = -s6*s234*c5 + c6*c234; T++;
    *T = -s5*s234; T++;
    *T = a2*s2 + a3*s23 + d1 - d5*c234 - d6*s5*s234; T++;

    *T = 0.0; T++; *T = 0.0; T++; *T = 0.0; T++; *T = 1.0;
  }

  void forward_6_back(const double* q, double* T) {
    const double d6_neg_div_2 = d6 / -2.0;

    double s1 = sin(*q), c1 = cos(*q); q++;
    double q23 = *q, q234 = *q, s2 = sin(*q), c2 = cos(*q); q++;
    double s3 = sin(*q), c3 = cos(*q); q23 += *q; q234 += *q; q++;
    double s4 = sin(*q), c4 = cos(*q); q234 += *q; q++;
    double s5 = sin(*q), c5 = cos(*q); q++;
    double s6 = sin(*q), c6 = cos(*q); 
    double s23 = sin(q23), c23 = cos(q23);
    double s234 = sin(q234), c234 = cos(q234);

    *T = (s1*s2 + c1*c5*c234)*c6 - s6*s234*c1; T++;
    *T = -(s1*s5 + c1*c5*c234)*s6 - s234*c1*c6; T++;
    *T = s1*s5 - s5*c1*c234; T++;
    *T = a2*c1*c2 + a3*c1*c23 + d4*s1 + d5*s234*c1 + d6_neg_div_2*s1*c5 - d6_neg_div_2*s5*c1*c234; T++;

    *T = (s1*c5*c234 - s5*c1) * c6 - s1*s6*s234; T++;
    *T = (-s1*c5*c234 + s5*c1)*c6 - s1*s234*c6; T++;
    *T = -s1*s5*c234 - c1*c5; T++;
    *T = a2*s1*c2 + a3*s1*c23 - d4*c1 + d5*s1*s234 - d6_neg_div_2*s1*s5*c234 - d6_neg_div_2*c1*c5; T++;

    *T = s6*c234 + s234*c5*c6; T++;
    *T = -s6*s234*c5 + c6*c234; T++;
    *T = -s5*s234; T++;
    *T = a2*s2 + a3*s23 + d1 - d5*c234 - d6_neg_div_2*s5*s234; T++;

    *T = 0.0; T++; *T = 0.0; T++; *T = 0.0; T++; *T = 1.0;
  }

  void forward_4(const double* q, double* T) {
    double s1 = sin(*q), c1 = cos(*q); q++;
    double q23 = *q, q234 = *q, s2 = sin(*q), c2 = cos(*q); q++;
    double s3 = sin(*q), c3 = cos(*q); q23 += *q; q234 += *q; q++;
    double s4 = sin(*q), c4 = cos(*q); q234 += *q; q++;
    double s5 = sin(*q), c5 = cos(*q); q++;
    double s6 = sin(*q), c6 = cos(*q); 
    double s23 = sin(q23), c23 = cos(q23);
    double s234 = sin(q234), c234 = cos(q234);

    *T = c1*c234; T++;
    *T = s1; T++;
    *T = s234*c1; T++;
    *T = a2*c1*c2 + a3*c1*c23 + d4*s1; T++;
    *T = s1*c234; T++;
    *T = -c1; T++;
    *T = s1*s234; T++;
    *T = a2*s1*c2 + a3*s1*c23 - d4*c1; T++;
    *T = s234; T++;
    *T = 0.0; T++;
    *T = -c234; T++;
    *T = a2*s2 + a3*s23 + d1; T++;
    *T = 0.0; T++; *T = 0.0; T++; *T = 0.0; T++; *T = 1.0;
  }

  void forward_3(const double *q, double *T) {
    double s1 = sin(*q), c1 = cos(*q); q++;
    double q23 = *q, q234 = *q, s2 = sin(*q), c2 = cos(*q); q++;
    double s3 = sin(*q), c3 = cos(*q); q23 += *q; q234 += *q; q++;
    double s4 = sin(*q), c4 = cos(*q); q234 += *q; q++;
    double s5 = sin(*q), c5 = cos(*q); q++;
    double s6 = sin(*q), c6 = cos(*q); 
    double s23 = sin(q23), c23 = cos(q23);

    *T = c1*c23; T++;
    *T = -s23*c1; T++;
    *T = s1; T++;
    *T = c1*(a2*c2 + a3*c23); T++;
    *T = c23*s1; T++;
    *T = -s1*s23; T++;
    *T = -c1; T++;
    *T = s1*(a2*c2 + a3*c23); T++;
    *T = s23; T++;
    *T = c23; T++;
    *T = 0.0; T++;
    *T = a2*s2 + a3*s23 + d1; T++;
    *T = 0.0; T++; *T = 0.0; T++; *T = 0.0; T++; *T = 1.0;
  }

  void forward_elbow_joint(const double* q, double* T) {
    const double s1 = sin(*q), c1 = cos(*q); q++; 
    const double s2 = sin(*q), c2 = cos(*q); q++;
    *T = c1*c2; T++;
    *T = -s2*c1; T++;
    *T = s1; T++;
    *T = a2*c1*c2; T++;
    *T = c2*s1; T++;
    *T = -s2*s1; T++;
    *T = -c1; T++;
    *T = a2*c2*s1; T++;
    *T = s2; T++;
    *T = c2; T++;
    *T = 0.0; T++;
    *T = a2*s2 + d1; T++;
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

  void jacobian_elbow_joint(double *jacobian, double *q) {
    const double s1 = sin(*q), c1 = cos(*q); q++; 
    const double s2 = sin(*q), c2 = cos(*q); q++;

    jacobian[0 * 6 + 0] = -a2 * s1 * c2;
    jacobian[0 * 6 + 1] = -a2 * s2 * c1;
    jacobian[0 * 6 + 2] = 0.0;
    jacobian[0 * 6 + 3] = 0.0;
    jacobian[0 * 6 + 4] = 0.0;
    jacobian[0 * 6 + 5] = 0.0;

    jacobian[1 * 6 + 0] = a2 * c1 * c2;
    jacobian[1 * 6 + 1] = -a2 * s1 * s2;
    jacobian[1 * 6 + 2] = 0.0;
    jacobian[1 * 6 + 3] = 0.0;
    jacobian[1 * 6 + 4] = 0.0;
    jacobian[1 * 6 + 5] = 0.0;

    jacobian[2 * 6 + 0] = 0.0;
    jacobian[2 * 6 + 1] = a2 * c2;
    jacobian[2 * 6 + 2] = 0.0;
    jacobian[2 * 6 + 3] = 0.0;
    jacobian[2 * 6 + 4] = 0.0;
    jacobian[2 * 6 + 5] = 0.0;
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

    const double sin23 = sin(q[1] + q[2]);
    const double cos23 = cos(q[1] + q[2]);

    const double sin234 = sin(q[1] + q[2] + q[3]);
    const double cos234 = cos(q[1] + q[2] + q[3]);

    jacobian[0 * 6 + 0] = -a2*sin1*cos2 - a3*sin1*cos23 + d4*cos1 - d5*sin1*sin234 + d6*sin1*sin5*cos234 + d6*cos1*cos5;
    jacobian[0 * 6 + 1] = (-a2*sin2 - a3*sin23 + d5*cos234 + d6*sin5*sin234)*cos1;
    jacobian[0 * 6 + 2] = (-a3*sin23 + d5*cos234 + d6*sin5*sin234)*cos1;
    jacobian[0 * 6 + 3] = (d5*cos234 + d6*sin5*sin234)*cos1;
    jacobian[0 * 6 + 4] = -d6*(sin1*sin5 + cos1*cos5*cos234);
    jacobian[0 * 6 + 5] = 0;

    jacobian[1 * 6 + 0] = a2*cos1*cos2 + a3*cos1*cos23 + d4*sin1 + d5*sin234*cos1 + d6*sin1*cos5 - d6*sin5*cos1*cos234;
    jacobian[1 * 6 + 1] = (-a2*sin2 - a3*sin23 + d5*cos234 + d6*sin5*sin234)*sin1;
    jacobian[1 * 6 + 2] = (-a3*sin23 + d5*cos234 + d6*sin5*sin234)*sin1;
    jacobian[1 * 6 + 3] = (d5*cos234 + d6*sin5*sin234)*sin1;
    jacobian[1 * 6 + 4] = d6*(-sin1*cos5*cos234 + sin5*cos1);
    jacobian[1 * 6 + 5] = 0;

    jacobian[2 * 6 + 0] = 0;
    jacobian[2 * 6 + 1] = a2*cos2 + a3*cos23 + d5*sin234 - d6*sin5*cos234;
    jacobian[2 * 6 + 2] = a3*cos23 + d5*sin234 - d6*sin5*cos234;
    jacobian[2 * 6 + 3] = d5*sin234 - d6*sin5*cos234;
    jacobian[2 * 6 + 4] = -d6*sin234*cos5;
    jacobian[2 * 6 + 5] = 0;

  }

  void jacobian_6_back(double *jacobian, double *q) {
    const double d6_neg_div_2 = d6 / -2.0;

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

    const double sin23 = sin(q[1] + q[2]);
    const double cos23 = cos(q[1] + q[2]);

    const double sin234 = sin(q[1] + q[2] + q[3]);
    const double cos234 = cos(q[1] + q[2] + q[3]);

    jacobian[0 * 6 + 0] = -a2*sin1*cos2 - a3*sin1*cos23 + d4*cos1 - d5*sin1*sin234 + d6_neg_div_2*sin1*sin5*cos234 + d6_neg_div_2*cos1*cos5;
    jacobian[0 * 6 + 1] = (-a2*sin2 - a3*sin23 + d5*cos234 + d6_neg_div_2*sin5*sin234)*cos1;
    jacobian[0 * 6 + 2] = (-a3*sin23 + d5*cos234 + d6_neg_div_2*sin5*sin234)*cos1;
    jacobian[0 * 6 + 3] = (d5*cos234 + d6_neg_div_2*sin5*sin234)*cos1;
    jacobian[0 * 6 + 4] = -d6_neg_div_2*(sin1*sin5 + cos1*cos5*cos234);
    jacobian[0 * 6 + 5] = 0;

    jacobian[1 * 6 + 0] = a2*cos1*cos2 + a3*cos1*cos23 + d4*sin1 + d5*sin234*cos1 + d6_neg_div_2*sin1*cos5 - d6_neg_div_2*sin5*cos1*cos234;
    jacobian[1 * 6 + 1] = (-a2*sin2 - a3*sin23 + d5*cos234 + d6_neg_div_2*sin5*sin234)*sin1;
    jacobian[1 * 6 + 2] = (-a3*sin23 + d5*cos234 + d6_neg_div_2*sin5*sin234)*sin1;
    jacobian[1 * 6 + 3] = (d5*cos234 + d6_neg_div_2*sin5*sin234)*sin1;
    jacobian[1 * 6 + 4] = d6_neg_div_2*(-sin1*cos5*cos234 + sin5*cos1);
    jacobian[1 * 6 + 5] = 0;

    jacobian[2 * 6 + 0] = 0;
    jacobian[2 * 6 + 1] = a2*cos2 + a3*cos23 + d5*sin234 - d6_neg_div_2*sin5*cos234;
    jacobian[2 * 6 + 2] = a3*cos23 + d5*sin234 - d6_neg_div_2*sin5*cos234;
    jacobian[2 * 6 + 3] = d5*sin234 - d6_neg_div_2*sin5*cos234;
    jacobian[2 * 6 + 4] = -d6_neg_div_2*sin234*cos5;
    jacobian[2 * 6 + 5] = 0;
  }

  void jacobian_4(double *jacobian, double *q) {
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

    const double sin23 = sin(q[1] + q[2]);
    const double cos23 = cos(q[1] + q[2]);

    jacobian[0 * 6 + 0] = -a2*sin1*cos2 - a3*sin1*cos23 + d4*cos1;
    jacobian[0 * 6 + 1] = -a2*sin2*cos1 - a3*sin23*cos1;
    jacobian[0 * 6 + 2] = -a3*sin23*cos1;
    jacobian[0 * 6 + 3] = 0;
    jacobian[0 * 6 + 4] = 0;
    jacobian[0 * 6 + 5] = 0;

    jacobian[1 * 6 + 0] = a2*cos1*cos2 + a3*cos1*cos23 + d4*sin1;
    jacobian[1 * 6 + 1] = -a2*sin2*sin1 - a3*sin23*sin1;
    jacobian[1 * 6 + 2] = -a3*sin23*sin1;
    jacobian[1 * 6 + 3] = 0;
    jacobian[1 * 6 + 4] = 0;  
    jacobian[1 * 6 + 5] = 0;

    jacobian[2 * 6 + 0] = 0;
    jacobian[2 * 6 + 1] = a2*cos2 + a3*cos23;
    jacobian[2 * 6 + 2] = a3*cos23;
    jacobian[2 * 6 + 3] = 0;
    jacobian[2 * 6 + 4] = 0;
    jacobian[2 * 6 + 5] = 0;
  }

  void jacobian_3(double *jacobian, double *q) {
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

    const double sin23 = sin(q[1] + q[2]);
    const double cos23 = cos(q[1] + q[2]);

    jacobian[0 * 6 + 0] = -(a2*cos2 + a3*cos23) * sin1;
    jacobian[0 * 6 + 1] = (-a2 * sin2 - a3 * sin23) * cos1;
    jacobian[0 * 6 + 2] = -a3*sin23*cos1;
    jacobian[0 * 6 + 3] = 0;
    jacobian[0 * 6 + 4] = 0;
    jacobian[0 * 6 + 5] = 0;

    jacobian[1 * 6 + 0] = (a2*cos2 + a3*cos23) * cos1;
    jacobian[1 * 6 + 1] = (-a2*sin2 - a3*sin23)*sin1;
    jacobian[1 * 6 + 2] = -a3*sin1*sin23;
    jacobian[1 * 6 + 3] = 0;
    jacobian[1 * 6 + 4] = 0;
    jacobian[1 * 6 + 5] = 0;

    jacobian[2 * 6 + 0] = 0;
    jacobian[2 * 6 + 1] = a2*cos2 + a3*cos23;
    jacobian[2 * 6 + 2] = a3*cos23;
    jacobian[2 * 6 + 3] = 0;
    jacobian[2 * 6 + 4] = 0;
    jacobian[2 * 6 + 5] = 0;
  }

};

std::tuple<double, double, double> forward_kinematics(double *q) {
  double *T = new double[16];
  for(auto i = 0; i < 16; i++)
    T[i] = 0.0;

  ur_kinematics::forward(q, T);
  
  double x = T[0*4 + 3], y = T[1*4 + 3], z = T[2*4 + 3];
  delete[] T;
  return {x, y, z};
}

std::tuple<double, double, double> forward_kinematics_6_back(double *q) {
  double *T = new double[16];
  for(auto i = 0; i < 16; i++)
    T[i] = 0.0;

  ur_kinematics::forward_6_back(q, T);

  double x = T[0*4 + 3], y = T[1*4 + 3], z = T[2*4 + 3];
  delete[] T;
  return {x, y, z};
}

std::tuple<double, double, double> forward_kinematics_4(double *q) {
  double *T = new double[16];
  for(auto i = 0; i < 16; i++)
    T[i] = 0.0;

  ur_kinematics::forward_4(q, T);

  double x = T[0*4 + 3], y = T[1*4 + 3], z = T[2*4 + 3];
  delete[] T;
  return {x, y, z};
}

std::tuple<double, double, double> forward_kinematics_3(double *q) {
  double *T = new double[16];
  for(auto i = 0; i < 16; i++)
    T[i] = 0.0;

  ur_kinematics::forward_3(q, T);

  double x = T[0*4 + 3], y = T[1*4 + 3], z = T[2*4 + 3];
  delete[] T;
  return {x, y, z};
}

std::tuple<double, double, double> forward_kinematics_elbow_joint(double *q) {
  double *T = new double[16];
  for(auto i = 0; i < 16; i++)
    T[i] = 0.0;

  ur_kinematics::forward_elbow_joint(q, T);

  double x = T[0*4 + 3], y = T[1*4 + 3], z = T[2*4 + 3];
  delete[] T;
  return {x, y, z};
}

int inverse_kinematics_2PI(double *q_sols, double x, double y, double z) {
  double q[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double *T = new double[16]{0};
  ur_kinematics::forward(q, T);

  std::tuple<double, double, double> target = {x, y, z};
  T[0*4 + 3] = std::get<0>(target);
  T[1*4 + 3] = std::get<1>(target);
  T[2*4 + 3] = std::get<2>(target);

  int num_sols;
  num_sols = ur_kinematics::inverse(T, q_sols);
  free(T);
  return num_sols;
}

#include <iostream>
int inverse_kinematics(double *q_sols, double x, double y, double z) {
  double *T = new double[16];
  for(auto i = 0; i < 16; i++)
    T[i] = 0.0;

  // ur_kinematics::forward_old(q, T);
  // for(auto row = 0; row < 3; row++) {
  //   for(auto col = 0; col < 3; col++) {
  //     std::cout << T[row*4 + col] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  double eerot[9] = {0, 0,  -1,  0, -1,  0,  -1, 0,  0};
  for(auto i = 0; i < 9; i++)
    T[i/3 * 4 + i % 3] = eerot[i];

  std::tuple<double, double, double> target = {x, y, z};
  T[0*4 + 3] = -std::get<0>(target);
  T[1*4 + 3] = -std::get<1>(target);
  T[2*4 + 3] = std::get<2>(target); 


  // for(auto row = 0; row < 2; row++) {
  //   for(auto col = 0; col < 4; col++) {
  //     T[row*4 + col] *= -1;
  //   }
  // }

  int num_sols;
  num_sols = ur_kinematics::inverse(T, q_sols);
  delete[] T;

  for (int i = 0; i < 6 * num_sols; ++i) {
    if (q_sols[i] >= ur_kinematics::PI)
        q_sols[i] -= 2.0 * ur_kinematics::PI;
  }

  return num_sols;
}

void joint_jacobian(double *jacobian, double *q) {
  ur_kinematics::jacobian(jacobian, q);
}

void joint_jacobian_6_back(double *jacobian, double *q) {
  ur_kinematics::jacobian_6_back(jacobian, q);
}

void joint_jacobian_4(double *jacobian, double *q) {
  ur_kinematics::jacobian_4(jacobian, q);
}

void joint_jacobian_3(double *jacobian, double *q) {
  ur_kinematics::jacobian_3(jacobian, q);
}

void jacobian_elbow_joint(double *jacobian, double *q) {
  ur_kinematics::jacobian_elbow_joint(jacobian, q);
}