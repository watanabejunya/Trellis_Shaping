#ifndef LIB_TS_H
#define LIB_TS_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <stdarg.h>
#include "lib_base.h"

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif

#ifndef AMPLITUDE_GAIN
#define AMPLITUDE_GAIN 1.0
#endif

#ifndef SMOOTHNESS_FACTOR
#define SMOOTHNESS_FACTOR 2.0
#endif


void ideally_linearized_power_amplifier(complex *s_in, complex *s_out, int n, double iob);
void solid_state_power_amplifier(complex *s_in, complex *s_out, int n, double iob);
void travelling_wave_tube_amplifier(complex *s_in, complex *s_out, int n, double iob);
void erf_based_power_amplifier(complex *s_in, complex *s_out, int n, double iob);
void print_in_out_amplitude(FILE *fp, complex *s_in, complex *s_out, int n);
void print_in_out_phase(FILE *fp, complex *s_in, complex *s_out, int n);
void print_gain_compression(FILE *fp, complex *s_in, complex *s_out, int n);
void calc_power_spectral_density(complex *x, double *psd, int n, int j, int m);


void ideally_linearized_power_amplifier(complex *s_in, complex *s_out, int n, double iob)
{
    const double k = 1.1220;                            // 1dB電力が落ちる入力電力dB
    double r_in_max;                                    // 入力振幅の最大値
    double r_out_max;                                   // 出力振幅の最大値
    double r, theta;                                    // 正規化瞬時振幅と位相
    int i;                                              // ループカウンタ

    // 入出力電力を設定する
    r_in_max = sqrt(pow(10.0, iob / 10.0) / pow(k, 2.0) * calc_average_power(s_in, n));
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    for (i = 0; i < n; i++) {
        r = cabs(s_in[i]);
        theta = carg(s_in[i]);

        s_out[i] = cexp(I * theta);

        if (r < r_in_max) {
            s_out[i] *= r_out_max * (r / r_in_max);
        } else {
            s_out[i] *= r_out_max;
        }
    }
}


void solid_state_power_amplifier(complex *s_in, complex *s_out, int n, double iob)
{
    const double p = SMOOTHNESS_FACTOR;                 // スムースネスファクター
    double k = 0.8745;                                  // 1dB電力が落ちる入力電力dB
    double r_in_max;                                    // 入力振幅の最大値
    double r_out_max;                                   // 出力振幅の最大値
    double r, theta;                                    // 正規化瞬時振幅と位相
    int i;                                              // ループカウンタ

    // 入出力電力を設定する
    r_in_max = sqrt(pow(10.0, iob / 10.0) / pow(k, 2.0) * calc_average_power(s_in, n));
    r_out_max = AMPLITUDE_GAIN * r_in_max;


    for (i = 0; i < n; i++) {
        r = cabs(s_in[i]);
        theta = carg(s_in[i]);

        s_out[i] = cexp(I * theta);
        s_out[i] *= r_out_max * (r / r_in_max) / pow(1 + pow(r / r_in_max, 2*p), 1/(2*p));
    }
}


void travelling_wave_tube_amplifier(complex *s_in, complex *s_out, int n, double iob)
{
    const double k = 0.6986;                            // 1dB電力が落ちる入力電力dB
    double r_in_max;                                    // 入力振幅の最大値
    double r_out_max ;                                  // 出力振幅の最大値
    double r, theta;                                    // 正規化瞬時振幅と位相
    int i;                                              // ループカウンタ

    // 入出力電力を設定する
    r_in_max = sqrt(pow(10.0, iob / 10.0) / pow(k, 2.0) * calc_average_power(s_in, n));
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    for (i = 0; i < n; i++) {
        r = cabs(s_in[i]);
        theta = carg(s_in[i]);

        s_out[i] = cexp(I * theta);
        s_out[i] *= r_out_max * (r / r_in_max) / (1.0 + pow(r / r_in_max, 2.0) / 4.0);
        s_out[i] *= cexp(I * (M_PI / 12.0) * pow(r / r_in_max, 2.0) / (1.0 + pow(r / r_in_max, 2.0) / 4.0));
    }
}


void erf_based_power_amplifier(complex *s_in, complex *s_out, int n, double iob)
{
    const double k = 0.6794;                            // 1dB電力が落ちる入力電力dB
    double r_in_max;                                    // 入力振幅の最大値
    double r_out_max ;                                  // 出力振幅の最大値
    double r, theta;                                    // 正規化瞬時振幅と位相
    int i;                                              // ループカウンタ

    // 入出力電力を設定する
    r_in_max = sqrt(pow(10.0, iob / 10.0) / pow(k, 2.0) * calc_average_power(s_in, n));
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    for (i = 0; i < n; i++) {
        r = cabs(s_in[i]);
        theta = carg(s_in[i]);

        s_out[i] = cexp(I * theta);
        s_out[i] *= r_out_max * erf((sqrt(M_PI) / 2.0) * (r / r_in_max));
    }
}


void print_in_out_amplitude(FILE *fp, complex *s_in, complex *s_out, int n)
{
    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%lf %lf\n", cabs(s_in[i]), cabs(s_out[i]));
        }
    }
}


void print_in_out_phase(FILE *fp, complex *s_in, complex *s_out, int n)
{
    const double degree_per_radian = 180.0 / M_PI;          // °/rad

    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%lf %lf\n", 20.0 * log10(cabs(s_in[i])), (carg(s_out[i]) - carg(s_in[i])) * degree_per_radian);
        }
    }
}


void print_gain_compression(FILE *fp, complex *s_in, complex *s_out, int n)
{
    double r_in_max = 1.0;

    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%lf %lf\n", 20.0 * log10(cabs(s_in[i]) / r_in_max), 20.0 * log10(cabs(s_out[i]) / (AMPLITUDE_GAIN * cabs(s_in[i]))));
        }
    }
}

#endif
