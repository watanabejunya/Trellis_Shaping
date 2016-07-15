#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "lib_base.h"
#include "lib_amp.h"

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif

#define NUM_ARGUMENT        2                       // 引数の数
#define STRIDE              100                     // グラフの刻み
#define INPUT_MAX           4                       // グラフの右端
#define SMOOTHNESS_FACTOR   2.0                     // スムースネスファクター

void run_ideal_power_amplifier()
{
    const double k = 1.1220;        // 1dB電力が落ちる入力電力dB
    complex *s_in, *s_out;          // 入出力信号
    double ibo;                     // Effective Input Back-off
    FILE *fp1, *fp2;                // 出力ファイル
    int i;                          // ループカウンタ

    // メモリの確保
    s_in = (complex *)malloc(STRIDE * INPUT_MAX * sizeof(complex));
    s_out = (complex *)malloc(STRIDE * INPUT_MAX * sizeof(complex));

    // ファイルオープン
    fp1 = fsopen("w", "./Result/AM_AM_convension(Ideal).dat");
    fp2 = fsopen("w", "./Result/gain_compression(Ideal).dat");

    // 信号の生成
    for (i = 0; i < STRIDE * INPUT_MAX; i++) {
        s_in[i] = (double)i / STRIDE;
    }

    // IBOの設定
    ibo = 10.0 * log10(pow(k, 2.0) / calc_average_power(s_in, STRIDE * INPUT_MAX));

    // アンプに通す
    ideally_linearized_power_amplifier(s_in, s_out, STRIDE * INPUT_MAX, ibo);

    print_in_out_amplitude(fp1, s_in, s_out, STRIDE * INPUT_MAX);
    print_gain_compression(fp2, s_in, s_out, STRIDE * INPUT_MAX);
}


void run_solid_state_power_amplifier()
{
    double k = 0.8745;              // 1dB電力が落ちる入力電力dB
    complex *s_in, *s_out;          // 入出力信号
    double ibo;                     // Effective Input Back-off
    FILE *fp1, *fp2;                // 出力ファイル
    int i;                          // ループカウンタ

    // メモリの確保
    s_in = (complex *)malloc(STRIDE * INPUT_MAX * sizeof(complex));
    s_out = (complex *)malloc(STRIDE * INPUT_MAX * sizeof(complex));

    // ファイルオープン
    fp1 = fsopen("w", "./Result/AM_AM_convension_p%lf(SSPA).dat", SMOOTHNESS_FACTOR);
    fp2 = fsopen("w", "./Result/gain_compression(SSPA).dat");

    for (i = 0; i < STRIDE * INPUT_MAX; i++) {
        s_in[i] = (double)i / STRIDE;
    }

    // IBOの設定
    ibo = 10.0 * log10(pow(k, 2.0) / calc_average_power(s_in, STRIDE * INPUT_MAX));

    // アンプに通す
    solid_state_power_amplifier(s_in, s_out, STRIDE * INPUT_MAX, SMOOTHNESS_FACTOR, ibo);

    print_in_out_amplitude(fp1, s_in, s_out, STRIDE * INPUT_MAX);
    print_gain_compression(fp2, s_in, s_out, STRIDE * INPUT_MAX);
}


void run_travelling_wave_tube_amplifier()
{
    const double k = 0.6986;        // 1dB電力が落ちる入力電力dB
    complex *s_in, *s_out;          // 入出力信号
    double ibo;                     // Effective Input Back-off
    FILE *fp1, *fp2, *fp3;          // 出力ファイル
    int i;                          // ループカウンタ

    // メモリの確保
    s_in = (complex *)malloc(STRIDE * INPUT_MAX * sizeof(complex));
    s_out = (complex *)malloc(STRIDE * INPUT_MAX * sizeof(complex));

    // ファイルオープン
    fp1 = fsopen("w", "./Result/AM_AM_convension(TWTA).dat");
    fp2 = fsopen("w", "./Result/AM_PM_curve(TWTA).dat");
    fp3 = fsopen("w", "./Result/gain_compression(TWTA).dat");

    for (i = 0; i < STRIDE * INPUT_MAX; i++) {
        s_in[i] = (double)i / STRIDE;
    }

    // IBOの設定
    ibo = 10.0 * log10(pow(k, 2.0) / calc_average_power(s_in, STRIDE * INPUT_MAX));


    // アンプに通す
    travelling_wave_tube_amplifier(s_in, s_out, STRIDE * INPUT_MAX, ibo);

    print_in_out_amplitude(fp1, s_in, s_out, STRIDE * INPUT_MAX);
    print_in_out_phase(fp2, s_in, s_out, STRIDE * INPUT_MAX);
    print_gain_compression(fp3, s_in, s_out, STRIDE * INPUT_MAX);
}


void run_erf_based_tube_amplifier()
{
    const double k = 0.6794;        // 1dB電力が落ちる入力電力dB
    complex *s_in, *s_out;          // 入出力信号
    double ibo;                     // Effective Input Back-off
    FILE *fp1, *fp2;                // 出力ファイル
    int i;                          // ループカウンタ

    // メモリの確保
    s_in = (complex *)malloc(STRIDE * INPUT_MAX * sizeof(complex));
    s_out = (complex *)malloc(STRIDE * INPUT_MAX * sizeof(complex));

    // ファイルオープン
    fp1 = fsopen("w", "./Result/AM_AM_convension_p%lf(Erf).dat", SMOOTHNESS_FACTOR);
    fp2 = fsopen("w", "./Result/gain_compression(Erf).dat");

    for (i = 0; i < STRIDE * INPUT_MAX; i++) {
        s_in[i] = (double)i / STRIDE;
    }

    // IBOの設定
    ibo = 20.0 * log10(k / 2.0);

    // アンプに通す
    erf_based_power_amplifier(s_in, s_out, STRIDE * INPUT_MAX, ibo);

    print_in_out_amplitude(fp1, s_in, s_out, STRIDE * INPUT_MAX);
    print_gain_compression(fp2, s_in, s_out, STRIDE * INPUT_MAX);
}


int main(int argc,char *argv[])
{
    // 入力チェック
    if (argc != NUM_ARGUMENT) {
        printf("Wrong number of arguments (%d for %d).\n", argc, NUM_ARGUMENT);
        exit(-1);
    }

    // 入力チェックと実行
    if (strcmp(argv[1], "ideal") == 0) {
        printf("Make ideal PA character graph.\n");
        run_ideal_power_amplifier();
    } else if (strcmp(argv[1], "sspa") == 0) {
        printf("Make SSPA character graph.\n");
        run_solid_state_power_amplifier();
    } else if (strcmp(argv[1], "twta") == 0) {
        printf("Make TWTA character graph.\n");
        run_travelling_wave_tube_amplifier();
    } else if (strcmp(argv[1], "erf") == 0) {
        printf("Make erf based PA character graph.\n");
        run_erf_based_tube_amplifier();
    } else {
        printf("Invalid argument\n");
        exit(-1);
    }
}
