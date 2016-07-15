#ifndef LIB_PEAK_H
#define LIB_PEAK_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "lib_base.h"

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif


double sum_square_autocor(complex *a, int n);
complex calc_autocor(complex *a, int n, int m);
double calc_papr_db(complex *a, int n);
double calc_normalized_ccdf(complex *a, int n, double threshold);
void count_papr_distribution(complex *x, double *pdf, int n, int min, int max, int num_index);
void count_normalized_distribution(complex *x, double *pdf, int n, int min, int max, int num_index);
double integrate_ccdf(double *pdf, double papr, int n, int min, int max, int num_index);


double sum_square_autocor(complex *a, int n)
{
    double autocor = 0;
    int m;

    for (m = 1; m < n; m++) {
        autocor += pow(cabs(calc_autocor(a, n, m)), 2.0);
    }

    return autocor;
}


complex calc_autocor(complex *a, int n, int m)
{
    complex autocor = 0;
    int i;

    for (i = 0; i < n - m; i++) {
        autocor += conj(a[i]) * a[i + m];
    }

    return autocor;
}

// PAPR(dB)を求める
double calc_papr_db (complex *a, int n) {
    int i;                              // ループカウンタ
    double peak = 0, average = 0;       // ピーク電力と平均電力

    for (i = 0; i < n; i++) {
        // 平均の計算
        average += pow(cabs(a[i]), 2.0) / (double)n;

        // ピークの更新
        if (peak < pow(cabs(a[i]), 2.0)) {
            peak = pow(cabs(a[i]), 2.0);
        }
    }

    // 比を返す
    return 10.0 * log10(peak / average);
}


// 正規化瞬時電力のCCDF
double calc_normalized_ccdf (complex *a, int n, double threshold) {
    int i;                              // ループカウンタ
    double average = 0;                 // ピーク電力と平均電力
    int count = 0;                      // カウンタ

    // 平均を求める
    for (i = 0; i < n; i++) {
        average += pow(cabs(a[i]), 2.0) / (double)n;
    }

    // CCDFを求める
    for (i = 0; i < n; i++) {
        if (10.0 * log10(pow(cabs(a[i]), 2.0) / average) > threshold) {
            count++;
        }
    }

    return count / (double)n;
}


// PAPRの分布を求める
void count_papr_distribution (complex *x, double *pdf, int n, int min, int max, int num_index) {
    double stride;                      // 刻み幅
    double papr;                        // PAPR
    int index;                          // 配列のインデックス
    static int memory_flag;             // メモリ管理フラグ
    int i;                              // ループカウンタ

    if (memory_flag == 0) {
        // 初期化
        for (i = 0; i < num_index; i++) {
            pdf[i] = 0;
        }

        memory_flag = 1;
    }

    // 配列１つの幅を設定
    stride = (double)(max - min) / (double)num_index;

    // PAPRを求める
    papr = calc_papr_db(x, n);

    // カウントする配列のインデックスを計算する
    index = (int)round((papr - (double)min) / stride);

    if (0 <= index && index < num_index) {
        pdf[index]++;
    }
}


// 正規化瞬時電力の分布を求める
void count_normalized_distribution (complex *x, double *pdf, int n, int min, int max, int num_index) {
    double stride;                      // 刻み幅
    double average;                     // ピーク電力と平均電力
    double ratio;                       // 平均電力との比
    int index;                          // 配列のインデックス
    static int memory_flag;             // メモリ管理フラグ
    int i;                              // ループカウンタ

    if (memory_flag == 0) {
        // 初期化
        for (i = 0; i < num_index; i++) {
            pdf[i] = 0;
        }

        memory_flag = 1;
    }

    // 配列１つの幅を設定
    stride = (double)(max - min) / (double)num_index;

    // 平均を求める
    average = calc_average_power(x, n);

    // CCDFを求める
    for (i = 0; i < n; i++) {
        // 平均電力との比を求める
        ratio = 10.0 * log10(pow(cabs(x[i]), 2.0) / average);

        // カウントする配列のインデックスを計算する
        index = (int)round((ratio - (double)min) / stride);

        if (0 <= index && index < num_index) {
            pdf[index] += 1.0 / (double)n;
        }
    }

}


// CCDFを計算する
double integrate_ccdf (double *pdf, double papr, int n, int min, int max, int num_index) {
    double stride;                      // 刻み幅
    int index;                          // 配列のインデックス
    double ccdf;                        // CCDF
    int i;                              // ループカウンタ

    // 配列１つの幅を設定
    stride = (double)(max - min) / (double)num_index;

    // 配列のインデックスを計算
    index = (int)round((papr - (double)min) / (double)stride);

    // CCDFを計算
    ccdf = 0;
    for (i = index; i < num_index; i++) {
        ccdf += pdf[i];
    }
    ccdf /= (double)n;

    return ccdf;
}


#endif
