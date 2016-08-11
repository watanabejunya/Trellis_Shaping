#ifndef LIB_AMP_H
#define LIB_AMP_H


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
double calc_output_back_off(complex *s_in, complex *s_out, int n, double iob);
double calc_amplifier_efficiency(complex *s_in, complex *s_out, int n, double iob, char *class);
void print_am_am_conversion(FILE *fp, complex *s_in, complex *s_out, int n);
void print_am_pm_conversion(FILE *fp, complex *s_in, complex *s_out, int n);
void print_gain_compression(FILE *fp, complex *s_in, complex *s_out, int n);
void calc_power_spectral_density(complex *x, double *psd, int n, int m);


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
        s_out[i] *= r_out_max * (r / r_in_max) / (1.0 + pow(r / r_in_max, 2.0) / 4.0);
        s_out[i] *= cexp(I * (M_PI / 12.0) * pow(r / r_in_max, 2.0) / (1.0 + pow(r / r_in_max, 2.0) / 4.0));
    }
}


void erf_based_power_amplifier(complex *s_in, complex *s_out, int n, double iob)
{
    const double k = 0.6794;                            // 1dB電力が落ちる入力電力dB
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
        s_out[i] *= r_out_max * erf((sqrt(M_PI) / 2.0) * (r / r_in_max));
    }
}


double calc_output_back_off(complex *s_in, complex *s_out, int n, double iob)
{
    double r_in_max;                                    // 入力振幅の最大値
    double r_out_max;                                   // 出力振幅の最大値
    double obo;                                         // OBO

    // 入出力電力を設定する
    r_in_max = sqrt(pow(10.0, iob / 10.0) * calc_average_power(s_in, n));
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    // 平均電力を求める
    obo = pow(r_out_max, 2.0) / calc_average_power(s_out, n);

    // OBOを返す
    return obo;
}


double calc_amplifier_efficiency(complex *s_in, complex *s_out, int n, double iob, char *class)
{
    double r_in_max;                                    // 入力振幅の最大値
    double r_out_max;                                   // 出力振幅の最大値
    double efficiency;                                  // PA効率

    // 入出力電力を設定する
    r_in_max = sqrt(pow(10.0, iob / 10.0) * calc_average_power(s_in, n));
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    if (strncmp(class, "a", 1) || strncmp(class, "A", 1)) {
        efficiency = (1.0 / 2.0) * calc_average_power(s_out, n) / pow(r_out_max, 2.0);
    } else if (strncmp(class, "b", 1) || strncmp(class, "B", 1)) {
        efficiency = (M_PI / 4.0) * calc_average_power(s_out, n) / (r_out_max * calc_average_amplitude(s_out, n));
    } else {
        fprintf(stderr, "invalid argument for pa class\n");
        exit(-1);
    }

    // PA効率を返す
    return efficiency;
}


void print_am_am_conversion(FILE *fp, complex *s_in, complex *s_out, int n)
{
    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%lf %lf\n", cabs(s_in[i]), cabs(s_out[i]));
        }
    }
}


void print_am_pm_conversion(FILE *fp, complex *s_in, complex *s_out, int n)
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


void calc_power_spectral_density(complex *x, double *psd, int n, int l)
{
    static int memory_flag;             // メモリ管理フラグ
    int i;                              // ループカウンタ

    if (memory_flag == 0) {
        // 初期化
        for (i = 0; i < n; i++) {
            psd[i] = 0;
        }

        memory_flag = 1;
    }

    // average_power = calc_average_power(x, n);
    for (i = 0; i < n/2; i++) {
        psd[i] += pow(cabs(x[i + n/2]), 2.0) / (double)l;
    }

    for (; i < n; i++) {
        psd[i] += pow(cabs(x[i - n/2]), 2.0) / (double)l;
    }
}
#endif
