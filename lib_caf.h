#ifndef LIB_CAF_H
#define LIB_CAF_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif


void clipping(complex *x, int n, double r);
void offset_attenuation(complex *x, int n, double r);


// クリッピング
void clipping (complex *x, int n, double r) {
    int i;                      // ループカウンタ
    double e = 0;               // 平均値
    double r_max;               // 振幅の最大値

    // 平均値を求める
    for (i = 0; i < n; i++) {
        e += pow(cabs(x[i]), 2.0) / (double)n;
    }

    // 閾値を求める
    r_max = sqrt(e) * r;

    // クリッピング
    for (i = 0; i < n; i++) {
        if (cabs(x[i]) > r_max) {
            x[i] *= r_max / cabs(x[i]);
        }
    }
}


// 減衰を補償する
void offset_attenuation (complex *x, int n, double r) {
    const double alpha = 1.0 - exp(-pow(r, 2.0)) + sqrt(M_PI)*r*erfc(r)/2.0;        // 減衰係数
    int i;                                                                          // ループカウンタ

    for (i = 0; i < n; i++) {
        x[i] /= alpha;
    }
}

#endif
