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

#define AMP_TYPE_IDEAL 1
#define AMP_TYPE_SSPA  2
#define AMP_TYPE_TWTA  3
#define AMP_TYPE_ERF   4

#define AMP_CLASS_A 1
#define AMP_CLASS_B 2

#ifndef AMPLITUDE_GAIN
#define AMPLITUDE_GAIN 1.0
#endif

#ifndef SMOOTHNESS_FACTOR
#define SMOOTHNESS_FACTOR 2.0
#endif

#ifndef IS_EFFECTIVE
#define IS_EFFECTIVE 0
#endif

complex ideally_linearized_power_amplifier(complex s_in, double r_in_max);
complex solid_state_power_amplifier(complex s_in, double r_in_max);
complex travelling_wave_tube_amplifier(complex s_in, double r_in_max);
complex erf_based_power_amplifier(complex s_in, double r_in_max);
double calc_output_back_off(complex *s_in, complex *s_out, int n, double iob);
double calc_pa_efficiency(complex *s_in, complex *s_out, int n, double iob, int class);
double calc_sdr_design(complex *s_in, complex *s_out, int n);
void print_am_am_conversion(FILE *fp, complex *s_in, complex *s_out, int n);
void print_am_pm_conversion(FILE *fp, complex *s_in, complex *s_out, int n);
void print_gain_compression(FILE *fp, complex *s_in, complex *s_out, int n);
void calc_power_spectral_density(complex *x, double *psd, int n, int m);


void amplify_signal(complex *s_in, complex *s_out, int n, double iob, int amp_type)
{
    double r_in_max;                                    // 入力振幅の最大値
    int i;                                              // ループカウンタ

    // 入出力電力を設定する
    r_in_max = sqrt(pow(10.0, iob / 10.0) * calc_average_power(s_in, n));

    for (i = 0; i < n; i++) {
        if (amp_type == AMP_TYPE_IDEAL) {
            s_out[i] = ideally_linearized_power_amplifier(s_in[i], r_in_max);
        } else if (amp_type == AMP_TYPE_SSPA) {
            s_out[i] = solid_state_power_amplifier(s_in[i], r_in_max);
        } else if (amp_type == AMP_TYPE_TWTA) {
            s_out[i] = travelling_wave_tube_amplifier(s_in[i], r_in_max);
        } else if (amp_type == AMP_TYPE_ERF) {
            s_out[i] = erf_based_power_amplifier(s_in[i], r_in_max);
        } else {
            fprintf(stderr, "Invalid Amplitude Type.\n");
            exit(-1);
        }
    }
}


complex ideally_linearized_power_amplifier(complex s_in, double r_in_max)
{
    const double k = 1.1220;                            // 1dB電力が落ちる入力電力dB
    complex s_out;                                      // 増幅後の信号
    double r_out_max;                                   // 出力振幅の最大値
    double r, theta;                                    // 正規化瞬時振幅と位相

    // 入出力電力を設定する
    if (IS_EFFECTIVE) {
         r_in_max /= pow(k, 2.0);
    }
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    r = cabs(s_in);
    theta = carg(s_in);

    s_out = cexp(I * theta);

    if (r < r_in_max) {
        s_out *= r_out_max * (r / r_in_max);
    } else {
        s_out *= r_out_max;
    }

    return s_out;
}


complex solid_state_power_amplifier(complex s_in, double r_in_max)
{
    const double p = SMOOTHNESS_FACTOR;                 // スムースネスファクター
    double k = 0.8745;                                  // 1dB電力が落ちる入力電力dB
    complex s_out;                                      // 増幅後の信号
    double r_out_max;                                   // 出力振幅の最大値
    double r, theta;                                    // 正規化瞬時振幅と位相

    // 入出力電力を設定する
    if (IS_EFFECTIVE) {
         r_in_max /= pow(k, 2.0);
    }
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    r = cabs(s_in);
    theta = carg(s_in);

    s_out = cexp(I * theta);
    s_out *= r_out_max * (r / r_in_max) / pow(1 + pow(r / r_in_max, 2*p), 1/(2*p));

    return s_out;
}


complex travelling_wave_tube_amplifier(complex s_in, double r_in_max)
{
    const double k = 0.6986;                            // 1dB電力が落ちる入力電力dB
    complex s_out;                                      // 増幅後の信号
    double r_out_max;                                   // 出力振幅の最大値
    double r, theta;                                    // 正規化瞬時振幅と位相

    // 入出力電力を設定する
    if (IS_EFFECTIVE) {
         r_in_max /= pow(k, 2.0);
    }
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    r = cabs(s_in);
    theta = carg(s_in);

    s_out = cexp(I * theta);
    s_out *= r_out_max * (r / r_in_max) / (1.0 + pow(r / r_in_max, 2.0) / 4.0);
    s_out *= cexp(I * (M_PI / 12.0) * pow(r / r_in_max, 2.0) / (1.0 + pow(r / r_in_max, 2.0) / 4.0));

    return s_out;
}


complex erf_based_power_amplifier(complex s_in, double r_in_max)
{
    const double k = 0.6794;                            // 1dB電力が落ちる入力電力dB
    complex s_out;                                      // 増幅後の信号
    double r_out_max;                                   // 出力振幅の最大値
    double r, theta;                                    // 正規化瞬時振幅と位相

    // 入出力電力を設定する
    if (IS_EFFECTIVE) {
         r_in_max /= pow(k, 2.0);
    }
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    r = cabs(s_in);
    theta = carg(s_in);

    s_out = cexp(I * theta);
    s_out *= r_out_max * erf((sqrt(M_PI) / 2.0) * (r / r_in_max));

    return s_out;
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


double calc_pa_efficiency(complex *s_in, complex *s_out, int n, double iob, int amp_class)
{
    double r_in_max;                                    // 入力振幅の最大値
    double r_out_max;                                   // 出力振幅の最大値
    double efficiency;                                  // PA効率

    // 入出力電力を設定する
    r_in_max = sqrt(pow(10.0, iob / 10.0) * calc_average_power(s_in, n));
    r_out_max = AMPLITUDE_GAIN * r_in_max;

    if (amp_class == AMP_CLASS_A) {
        efficiency = (1.0 / 2.0) * calc_average_power(s_out, n) / pow(r_out_max, 2.0);
    } else if (amp_class == AMP_CLASS_B) {
        efficiency = (M_PI / 4.0) * calc_average_power(s_out, n) / (r_out_max * calc_average_amplitude(s_out, n));
    } else {
        fprintf(stderr, "invalid argument for power amplifier class\n");
        exit(-1);
    }

    // PA効率を返す
    return efficiency;
}


double calc_sdr_design(complex *s_in, complex *s_out, int n)
{
    double numerator = 0;                               // 相関係数の分子
    double rho;                                         // 相関係数
    double sdr;                                         // SDR
    int i;                                              // ループカウンタ

    // 相関係数の分子を計算する
    for (i = 0; i < n; i++) {
        numerator += conj(s_in[i]) * s_out[i] / (double)n;
    }

    // 相関係数を計算する
    rho = numerator / sqrt(calc_average_power(s_in, n) * calc_average_power(s_out, n));

    sdr = pow(rho, 2.0) / (1 - pow(rho, 2.0));

    // PA効率を返す
    return 10.0 * log10(sdr);
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
    double max_power = 0;               // 最大電力
    int i;                              // ループカウンタ

    if (memory_flag == 0) {
        // 初期化
        for (i = 0; i < n; i++) {
            psd[i] = 0;
        }

        memory_flag = 1;
    }

    for (i = 0; i < n; i++) {
        if (pow(cabs(x[i]), 2.0) > max_power) {
            max_power = pow(cabs(x[i]), 2.0);
        }
    }

    for (i = 0; i < n/2; i++) {
        psd[i] += 10.0 * log10(pow(cabs(x[i + n/2]), 2.0) / max_power) / (double)l;
    }

    for (; i < n; i++) {
        psd[i] += 10.0 * log10(pow(cabs(x[i - n/2]), 2.0) / max_power) / (double)l;
    }
}
#endif
