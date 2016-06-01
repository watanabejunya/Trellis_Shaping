#ifndef LIB_TS_H
#define LIB_TS_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <stdarg.h>

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif

#define FFTJ_FORWORD 1                      // FFTの向き(FFT)
#define FFTJ_BACKWORD -1                    // FFTの向き(IFFT)

/*
 * 変数
 */

extern int count_add;
extern int count_mul;

/*
 * 構造体
 */

// TSのノード
typedef struct {
    int *a_index;
    complex *autocor;
    double metric;
    int pre_state;
} node;

// FFTJの設定
typedef struct {
    int size;                                   // FFTサイズ
    int *index;                                 // ビットリバースした配列のインデックス
    complex *w;                                 // 回転演算子
    complex **data;                             // 計算用データ配列
} fftj_plan;

/*
 * 関数宣言
 */

FILE *fsopen(const char *mode, const char *format, ...);
double dadd(double a, double b);
double dmul(double a, double b);
double dsqr(double a);
complex cadd(complex a, complex b);
complex cmul(complex a, complex b);
double cmag(complex a);
double cdis(complex a, complex b);
complex cipr(complex a, complex b);
fftj_plan fftj_plan_dft(int, complex *, complex *, int);
void fftj_execute(fftj_plan);
void fftj_destroy(fftj_plan);
void fftj(int, complex *, complex *);
void ifftj(int, complex *, complex *);
void make_signal(int *x, int n);
void gaussian_noise(complex *x, double sigma, int n);
int count_be(int *x, int *y, int n);
void print_map(FILE *fp, fftw_complex *x, int n);
void copy_int (int *x, int *x_copy, int n);
void copy_complex (complex *x, complex *x_copy, int n);
void xor_addition (int *x, int *y, int n);
void demultiplexer (int *d, int *s, int *b, int num_d, int num_s, int num_b, int n);
void multiplexer (int *s, int *b, int *d, int num_s, int num_b, int num_d, int n);
void convolutional_encoding (int *u, int *c, int n);
void convolutional_encoding2 (int *u, int *c, int n);
void convolutional_encoding3 (int *u, int *c, int n);
void parity_check_decoding (int *c, int *u, int n);
void parity_check_decoding2 (int *c, int *u, int n);
void parity_check_decoding3 (int *c, int *u, int n);
void inverse_parity_check_encoding (int *u, int *c, int n);
void inverse_parity_check_encoding2 (int *u, int *c, int n);
void inverse_parity_check_encoding3 (int *u, int *c, int n);
complex map_4_4_qam (int bit1, int bit2);
complex map_4_16_qam (int bit1, int bit2);
complex map_4_64_qam (int bit1, int bit2);
complex map_16qam_type1 (int bit1, int bit2, int bit3, int bit4);
complex map_16qam_type2 (int bit1, int bit2, int bit3, int bit4);
complex map_16_4_qam (int bit1, int bit2, int bit3, int bit4);
complex map_16_16_qam (int bit1, int bit2, int bit3, int bit4);
complex map_64qam_type1 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);
complex map_64qam_type2 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);
complex map_64_4_qam (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);
complex map_256qam_type1 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8);
complex map_256qam_type2 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8);
void unmap_qpsk (complex a, int *bits);
void unmap_16qam_type1 (complex a, int *bits);
void unmap_16qam_type2 (complex a, int *bits);
void unmap_64qam_type1 (complex a, int *bits);
void unmap_64qam_type2 (complex a, int *bits);
void unmap_256qam_type1 (complex a, int *bits);
void unmap_256qam_type2 (complex a, int *bits);
void construct_constellation (complex *constellation, int m, int type);
void qam_modulation_lsb (int *c, complex *a, int n, int m);
void qam_modulation_lsb2 (int *c, complex *a, int n, int m);
void qam_modulation_lsb3 (int *c, complex *a, int n, int m);
void qam_demodulation (complex *a, int *c, int n, int m, int type);
void fft (int n, fftw_complex *in, fftw_complex *out);
void ifft (int n, fftw_complex *in, fftw_complex *out);
void over_sampling (complex *x, fftw_complex *y, int j, int n);
void down_sampling (fftw_complex *x, complex *y, int j, int n);
void clipping (complex *x, int n, double r);
void offset_attenuation (complex *x, int n, double r);
void trellis_shaping (int *c, complex *a, int n, int num_qam, int type);
void trellis_shaping_caf (int *c, complex *a_caf, complex *a, int n, int m);
void trellis_shaping_caf2 (int *c, complex *a_caf, complex *a, int n, int m);
void trellis_shaping_caf3 (int *c, complex *a_caf, complex *a, int n, int m);
double calc_average_power(complex *a, int n);
double calc_papr_db (complex *a, int n);
double calc_normalized_ccdf (complex *a, int n, double threshold);
void count_papr_distribution (complex *x, double *pdf, int n, int min, int max, int num_index);
void count_normalized_distribution (complex *x, double *pdf, int n, int min, int max, int num_index);
double integrate_ccdf (double *pdf, double papr, int n, int min, int max, int num_index);

/**
 * 関数記述
 */

// ファイルオープン
FILE *fsopen(const char *mode, const char *format, ...) {
  va_list list;
  va_start(list, format);

  char fname[256];
  vsprintf(fname, format, list);

  FILE *file;
  file = fopen(fname, mode);

  va_end(list);

  return file;
}


double dadd(double a, double b)
{
    count_add++;
    return a + b;
}


double dmul(double a, double b)
{
    count_mul++;
    return a * b;
}


double dsqr(double a)
{
    return dmul(a, a);
}


complex cadd(complex a, complex b)
{
    return dadd(creal(a), creal(b)) + I * dadd(cimag(a), cimag(b));
}


complex cmul(complex a, complex b)
{
    return dadd(dmul(creal(a), creal(b)), - dmul(cimag(a), cimag(b))) + I * dadd(dmul(creal(a), cimag(b)), dmul(cimag(a), creal(b)));
}


double cmag(complex a)
{
    return dadd(dsqr(creal(a)), dsqr(cimag(a)));
}


double cdis(complex a, complex b)
{
    return cmag(cadd(a, - b));
}

complex cipr(complex a, complex b)
{
    return cmul(conj(a), b);
}


/**
 * FFTの回転因子，ビット逆順を計算し，FFTの設定を行う
 *
 * @param  int                  n            FFTサイズ
 * @param  _Complex double*     in           入力データ
 * @param  _Complex double*     out          出力データ
 * @param  int                  type         FFTかIFFTかの選択
 * @return fftj_plan                         FFTの設定
 */
fftj_plan fftj_plan_dft(int n, complex *in, complex *out, int type)
{
    const int exponent = (int) log2(n);             // nの冪指数
    const double sign = type == 1 ? -1.0 : 1.0;     // 1.0 or -1.0
    fftj_plan plan;                                 // FFTの設定
    int i, j;                                       // ループカウンタ
    int start, end;                                 // ビットリバースループで使う作業用変数

    // 入力チェック
    if (n < 2 || n != (int) pow(2.0, exponent)) {
        fprintf(stderr, "n should be the k-th power of 2.\n");
        exit(-1);
    }
    if (type != 1 && type != -1) {
        fprintf(stderr, "type should be 1 or -1.\n");
        exit(-1);
    }

    // サイズの設定
    plan.size = n;

    // メモリ確保
    plan.index = (int *) malloc(n * sizeof(int));
    plan.w = (complex *) malloc(n * sizeof(complex));
    plan.data = (complex **) malloc((exponent+1) * sizeof(complex *));
    for (i = 0; i < exponent; i++) {
        plan.data[i] = (complex *) malloc(n * sizeof(complex));
    }
    plan.data[exponent] = out;

    // ビットリバースしたときのインデックスを作っておく
    plan.index[0] = 0;
    for (i = 1; i <= exponent; i++) {
        // ループの範囲を定める
        start = (int)pow(2.0, i - 1);
        end = (int)pow(2.0, i);

        for (j = start; j < end; j++) {
            // ビットリバースした値を求める
            plan.index[j] = plan.index[j - start] + n / end;
        }
    }

    // 回転因子の計算
    for (i = 0; i < n; i++) {
        plan.w[i] = cexp(sign * I * 2.0 * M_PI * i / n);
    }

    // 入力をバタフライ演算用に入れ替えてセットしておく
    for (i = 0; i < n; i++) {
        plan.data[0][i] = in[plan.index[i]];
    }

    return plan;
}


/**
 * FFT(IFFT)を実行する
 *
 * @param plan      FFTの設定
 */
void fftj_execute(fftj_plan plan)
{
    const int exponent = (int) log2(plan.size);         // FFTサイズの冪指数
    int j_roop, k_roop, half;                           // バタフライ演算で使う作業用変数
    int tmp_index1, tmp_index2;                         // 一時的な配列のインデックス
    int i, j, k;                                        // ループカウンタ

    for (i = 1; i <= exponent; i++) {
        // ループの範囲を設定
        j_roop = plan.size / (int)pow(2.0, i);
        k_roop = (int)pow(2.0, i);
        half = k_roop / 2;

        for (j = 0; j < j_roop; j++) {
            for (k = 0; k < k_roop; k++) {
                // 奇数行と偶数行でバタフライの演算の向きを変える
                if (k < half) {
                    tmp_index1 = j * k_roop + k;
                    tmp_index2 = j * k_roop + k + half;
                } else {
                    tmp_index1 = j * k_roop + k - half;
                    tmp_index2 = j * k_roop + k;
                }
                plan.data[i][j * k_roop + k] = cadd(plan.data[i-1][tmp_index1], cmul(plan.data[i-1][tmp_index2], plan.w[k * j_roop]));
            }
        }
    }

}


/**
 * FFTの際に使ったメモリを解放する
 *
 * @param   fftj_plan   plan        FFTの設定
 */
void fftj_destroy(fftj_plan plan)
{
    const int exponent = (int) log2(plan.size);             // nの冪指数
    int i;                                          // ループカウンタ

    free(plan.index);
    free(plan.w);
    for (i = 0; i < exponent; i++) {
        free(plan.data[i]);
    }
    free(plan.data);
}


/**
 * FFTを行い正規化する
 * @param  int                  n            FFTサイズ
 * @param  _Complex double*     in           入力データ
 * @param  _Complex double*     out          出力データ
 */
void fftj(int n, complex *in, complex *out)
{
    fftj_plan plan;                     // FFTの設定
    int i;                              // ループカウンタ

    plan = fftj_plan_dft(n, in, out, FFTJ_FORWORD);

    fftj_execute(plan);

    // データを正規化
    for (i = 0; i < n; i++) {
        out[i] /= sqrt(n);
    }

    fftj_destroy(plan);
}

/**
 * IFFTを行い正規化する
 * @param  int                  n            FFTサイズ
 * @param  _Complex double*     in           入力データ
 * @param  _Complex double*     out          出力データ
 */
void ifftj(int n, complex *in, complex *out)
{
    fftj_plan plan;                     // FFTの設定
    int i;                              // ループカウンタ

    plan = fftj_plan_dft(n, in, out, FFTJ_BACKWORD);

    fftj_execute(plan);

    // データを正規化
    for (i = 0; i < n; i++) {
        out[i] /= sqrt(n);
    }

    fftj_destroy(plan);
}


// 0,1信号生成
void make_signal(int *x, int n) {
    for (int i = 0; i < n; i++) {
        x[i] = random() % 2;
    }
}


// ガウス雑音
void gaussian_noise(complex *x, double sigma, int n) {

    double u1, u2;
    complex v;

    // ボックスミュラー法
    for (int i = 0; i < n; i++) {
        u1 = ((double)random() + 1.0) / ((double)RAND_MAX + 1.0);
        u2 = ((double)random() + 1.0) / ((double)RAND_MAX + 1.0);

        v = sqrt(-2.0*log(u1)) * cos(2.0*M_PI*u2) + I * sqrt(-2.0*log(u1)) * sin(2.0*M_PI*u2);

        x[i] += v * sigma;
    }
}


// ビットエラーを数える
int count_be(int *x, int *y, int n) {
    int count= 0;

    for (int i = 0; i < n; i++) {
        if (x[i] != y[i]) {
            count++;
        }
    }

    return count;
}


// マッピング出力
void print_map(FILE *fp, complex *x, int n) {
    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%f %f\n", creal(x[i]), cimag(x[i]));
        }
    }
}


// 整数をコピー
void copy_int (int *x, int *x_copy, int n) {
    int i;          // ループカウンタ

    for (i = 0; i < n; i++) {
        x_copy[i] = x[i];
    }
}

// 複素数をコピー
void copy_complex (complex *x, complex *x_copy, int n) {
    int i;          // ループカウンタ

    for (i = 0; i < n; i++) {
        x_copy[i] = x[i];
    }
}


// インタリーバ
void xor_addition (int *x, int *y, int n) {
    int i;

    for (i = 0; i < n; i++) {
        x[i] = x[i] ^ y[i];
    }
}


//デマルチプレクサ(num_d → num_s + num_b)
void demultiplexer (int *d, int *s, int *b, int num_d, int num_s, int num_b, int n) {
    if (num_d != num_s + num_b) {
        fprintf(stderr, "%d bits cannot be demultiplexed to %d and %d bits.\n", num_d, num_s, num_b);
        exit(-1);
    }

    int i, j, k;                    // ループカウンタ

    for(i = 0; i < n; i++){
        for (j = 0; j < num_s; j++) {
            s[i * num_s + j] = d[i * num_d + j];
        }

        for (k = 0; k < num_b; k++, j++) {
            b[i * num_b + k] = d[i * num_d + j];
        }
    }
}


//マルチプレクサ(num_s + num_b → num_d)
void multiplexer (int *s, int *b, int *d, int num_s, int num_b, int num_d, int n) {
    if (num_d != num_s + num_b) {
        fprintf(stderr, "%d bits cannot be demultiplexed to %d and %d bits.\n", num_d, num_s, num_b);
        exit(-1);
    }

    int i, j, k;                    // ループカウンタ

    for(i = 0; i < n; i++){
        for (j = 0; j < num_s; j++) {
             d[i * num_d + j] = s[i * num_s + j];
        }

        for (k = 0; k < num_b; k++, j++) {
            d[i * num_d + j] = b[i * num_b + k];
        }
    }
}


// 畳み込み符号化(拘束長3)
void convolutional_encoding (int *u, int *c, int n) {
    int i;                              // ループカウンタ
    int register1, register2;           // シフトレジスタ


    // シフトレジスタを初期化
    register1 = 0;
    register2 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[2*i] = u[i] ^ register2;
        c[2*i+1] = u[i] ^ register1 ^ register2;

        // レジスタをシフト
        register2 = register1;
        register1 = u[i];
    }
}


// 畳み込み符号化(拘束長3)
void convolutional_encoding2 (int *u, int *c, int n) {
    int i;                              // ループカウンタ
    int register11, register12;         // シフトレジスタ
    int register21, register22;         // シフトレジスタ


    // シフトレジスタを初期化
    register11 = 0;
    register12 = 0;
    register21 = 0;
    register22 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[4*i] = u[2*i] ^ register12;
        c[4*i+1] = u[2*i] ^ register11 ^ register12;
        c[4*i+2] = u[2*i+1] ^ register22;
        c[4*i+3] = u[2*i+1] ^ register21 ^ register22;

        // レジスタをシフト
        register12 = register11;
        register11 = u[2*i];
        register22 = register21;
        register21 = u[2*i+1];
    }
}


// 畳み込み符号化(拘束長3)
void convolutional_encoding3 (int *u, int *c, int n) {
    int i;                              // ループカウンタ
    int register11, register12;           // シフトレジスタ
    int register21, register22;         // シフトレジスタ
    int register31, register32;         // シフトレジスタ


    // シフトレジスタを初期化
    register11 = 0;
    register12 = 0;
    register21 = 0;
    register22 = 0;
    register31 = 0;
    register32 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[6*i] = u[3*i] ^ register12;
        c[6*i+1] = u[3*i] ^ register11 ^ register12;
        c[6*i+2] = u[3*i+1] ^ register22;
        c[6*i+3] = u[3*i+1] ^ register21 ^ register22;
        c[6*i+4] = u[3*i+2] ^ register32;
        c[6*i+5] = u[3*i+2] ^ register31 ^ register32;

        // レジスタをシフト
        register12 = register11;
        register11 = u[3*i];
        register22 = register21;
        register21 = u[3*i+1];
        register32 = register31;
        register31 = u[3*i+2];
    }
}


// パリティ検査行列による復号(拘束長3)
void parity_check_decoding (int *c, int *u, int n) {
    int i;                                                      // ループカウンタ
    int register1, register2, register3, register4;         // シフトレジスタ

    // シフトレジスタを初期化
    register1 = 0;
    register2 = 0;
    register3 = 0;
    register4 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        u[i] = c[2*i] ^ register1 ^ register2 ^ c[2*i+1] ^ register4;

        // レジスタをシフト
        register4 = register3;
        register3 = c[2*i+1];
        register2 = register1;
        register1 = c[2*i];
    }
}


// パリティ検査行列による復号(拘束長3)
void parity_check_decoding2 (int *c, int *u, int n) {
    int i;                                                      // ループカウンタ
    int register11, register12, register13, register14;         // シフトレジスタ
    int register21, register22, register23, register24;         // シフトレジスタ

    // シフトレジスタを初期化
    register11 = 0;
    register12 = 0;
    register13 = 0;
    register14 = 0;
    register21 = 0;
    register22 = 0;
    register23 = 0;
    register24 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        u[2*i] = c[4*i] ^ register11 ^ register12 ^ c[4*i+1] ^ register14;
        u[2*i+1] = c[4*i+2] ^ register21 ^ register22 ^ c[4*i+3] ^ register24;

        // レジスタをシフト
        register14 = register13;
        register13 = c[4*i+1];
        register12 = register11;
        register11 = c[4*i];
        register24 = register23;
        register23 = c[4*i+3];
        register22 = register21;
        register21 = c[4*i+2];
    }
}


// パリティ検査行列による復号(拘束長3)
void parity_check_decoding3 (int *c, int *u, int n) {
    int i;                                                      // ループカウンタ
    int register11, register12, register13, register14;
    int register21, register22, register23, register24;          // シフトレジスタ
    int register31, register32, register33, register34;          // シフトレジスタ

    // シフトレジスタを初期化
    register11 = 0;
    register12 = 0;
    register13 = 0;
    register14 = 0;
    register21 = 0;
    register22 = 0;
    register23 = 0;
    register24 = 0;
    register31 = 0;
    register32 = 0;
    register33 = 0;
    register34 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        u[3*i] = c[6*i] ^ register11 ^ register12 ^ c[6*i+1] ^ register14;
        u[3*i+1] = c[6*i+2] ^ register21 ^ register22 ^ c[6*i+3] ^ register24;
        u[3*i+2] = c[6*i+4] ^ register31 ^ register32 ^ c[6*i+5] ^ register34;

        // レジスタをシフト
        register14 = register13;
        register13 = c[6*i+1];
        register12 = register11;
        register11 = c[6*i];
        register24 = register23;
        register23 = c[6*i+3];
        register22 = register21;
        register21 = c[6*i+2];
        register34 = register33;
        register33 = c[6*i+5];
        register32 = register31;
        register31 = c[6*i+4];
    }
}


// パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_encoding (int *u, int *c, int n) {
    int i;                              // ループカウンタ
    int register1;                      // シフトレジスタ

    // シフトレジスタを初期化
    register1 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[2*i] = register1;
        c[2*i+1] = u[i] ^ register1;

        // レジスタをシフト
        register1 = u[i];
    }
}


// 並列パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_encoding2 (int *u, int *c, int n) {
    int i;                                  // ループカウンタ
    int register1, register2;               // シフトレジスタ

    // シフトレジスタを初期化
    register1 = 0;
    register2 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[4*i] = register1;
        c[4*i+1] = u[2*i] ^ register1;

        c[4*i+2] = register2;
        c[4*i+3] = u[2*i+1] ^ register2;

        // レジスタをシフト
        register1 = u[2*i];
        register2 = u[2*i+1];
    }
}


// 並列パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_encoding3 (int *u, int *c, int n) {
    int i;                                      // ループカウンタ
    int register1, register2, register3;        // シフトレジスタ

    // シフトレジスタを初期化
    register1 = 0;
    register2 = 0;
    register3 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[6*i] = register1;
        c[6*i+1] = u[3*i] ^ register1;

        c[6*i+2] = register2;
        c[6*i+3] = u[3*i+1] ^ register2;

        c[6*i+4] = register3;
        c[6*i+5] = u[3*i+2] ^ register3;

        // レジスタをシフト
        register1 = u[3*i];
        register2 = u[3*i+1];
        register3 = u[3*i+2];
    }
}


// 下位ビットQPSKマッピング(16QAM対応)
complex map_4_4_qam (int bit1, int bit2) {
    complex a;                      // 変調信号

    a = (2.0 - 4.0 * bit2) + I*(2.0 - 4.0 * bit1);

    return a;
}


// 下位ビットQPSKマッピング(64QAM対応)
complex map_4_16_qam (int bit1, int bit2) {
    complex a;                      // 変調信号

    a = (4.0 - 8.0 * bit2) + I*(4.0 - 8.0 * bit1);

    return a;
}


// 下位ビットQPSKマッピング(256QAM対応)
complex map_4_64_qam (int bit1, int bit2) {
    complex a;                      // 変調信号

    a = (8.0 - 16.0 * bit2) + I*(8.0 - 16.0 * bit1);

    return a;
}


// 16QAMマッピング(type1)
complex map_16qam_type1 (int bit1, int bit2, int bit3, int bit4) {
    complex a;                      // 変調信号

    a = (3.0 - 2.0*bit4) + I*(3.0 - 2.0*bit3);
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }

    return a;
}


// 16QAMマッピング(Type2)
complex map_16qam_type2 (int bit1, int bit2, int bit3, int bit4) {
    complex a;                      // 変調信号

    a = (2.0 - 4.0 * bit2) + I*(2.0 - 4.0 * bit1);
    a += (1.0 - 2.0 * bit4) + I*(1.0 - 2.0 * bit3);

    return a;
}


// 下位ビット16QAMマッピング(64QAM対応)
complex map_16_4_qam (int bit1, int bit2, int bit3, int bit4) {
    complex a;                      // 変調信号

    a = map_4_4_qam(bit3, bit4);
    a += 4.0 + I*4.0;
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }

    return a;
}


// 下位ビット16QAMマッピング(256QAM対応)
complex map_16_16_qam (int bit1, int bit2, int bit3, int bit4) {
    complex a;                      // 変調信号

    a = map_4_16_qam(bit3, bit4);
    a += 8.0 + I*8.0;
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }

    return a;
}


// 64QAMマッピング(Type1)
complex map_64qam_type1 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6) {
    complex a;                      // 変調信号

    a = map_16qam_type1(bit3, bit4, bit5, bit6);
    a += 4.0 + I*4.0;
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }

    return a;
}


// 64QAMマッピング(Type2)
complex map_64qam_type2 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6) {
    complex a;                      // 変調信号

    a = map_16qam_type1(bit3, bit4, bit5, bit6);
    a += (4.0 - 8.0 * bit2) + I*(4.0 - 8.0 * bit1);

    return a;
}


// 下位ビット64QAMマッピング(256QAM対応)
complex map_64_4_qam (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6) {
    complex a;                      // 変調信号

    a = map_16_4_qam(bit3, bit4, bit5, bit6);
    a += 8.0 + I*8.0;
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }


    return a;
}


// 256QAMマッピング
complex map_256qam_type1 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8) {
    complex a;                      // 変調信号

    a = map_64qam_type1(bit3, bit4, bit5, bit6, bit7, bit8);
    a += 8.0 + I*8.0;
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }

    return a;
}


// 256QAMマッピング(Type2)
complex map_256qam_type2 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8) {
    complex a;                      // 変調信号

    a = map_64qam_type1(bit3, bit4, bit5, bit6, bit7, bit8);
    a += (8.0 - 16.0 * bit2) + I*(8.0 - 16.0 * bit1);

    return a;
}


// QPSKアンマッピング
void unmap_qpsk (complex a, int *bits) {

    if (creal(a) > 0) {
        bits[1] = 0;
        if (cimag(a) > 0) {
            bits[0] = 0;
        } else {
            bits[0] = 1;
        }
    } else {
        bits[1] = 1;
        if (cimag(a) > 0) {
            bits[0] = 0;
        } else {
            bits[0] = 1;
        }
    }
}


// 16QAMアンマッピング(Type1)
void unmap_16qam_type1 (complex a, int *bits) {
    if (creal(a) > 0) {
        bits[1] = 0;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += -2.0 - I * 2.0;
        } else {
            bits[0] = 1;
            a += -2.0 + I * 2.0;
            a = conj(a);
        }
    } else {
        bits[1] = 1;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += +2.0 - I * 2.0;
            a = -conj(a);
        } else {
            bits[0] = 1;
            a += +2.0 + I * 2.0;
            a = -a;
        }
    }

    unmap_qpsk(a, bits + 2);
}


// 16QAMアンマッピング(Type2)
void unmap_16qam_type2 (complex a, int *bits) {
    if (creal(a) > 0) {
        bits[1] = 0;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += -2.0 - I * 2.0;
        } else {
            bits[0] = 1;
            a += -2.0 + I * 2.0;
        }
    } else {
        bits[1] = 1;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += +2.0 - I * 2.0;
        } else {
            bits[0] = 1;
            a += +2.0 + I * 2.0;
        }
    }

    unmap_qpsk(a, bits + 2);
}


// 64QAMアンマッピング(Type1)
void unmap_64qam_type1 (complex a, int *bits) {
    if (creal(a) > 0) {
        bits[1] = 0;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += -4.0 - I * 4.0;
        } else {
            bits[0] = 1;
            a += -4.0 + I * 4.0;
            a = conj(a);
        }
    } else {
        bits[1] = 1;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += +4.0 - I * 4.0;
            a = -conj(a);
        } else {
            bits[0] = 1;
            a += +4.0 + I * 4.0;
            a = -a;
        }
    }

    unmap_16qam_type1(a, bits + 2);
}


// 64QAMアンマッピング(Type2)
void unmap_64qam_type2 (complex a, int *bits) {
    if (creal(a) > 0) {
        bits[1] = 0;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += -4.0 - I * 4.0;
        } else {
            bits[0] = 1;
            a += -4.0 + I * 4.0;
        }
    } else {
        bits[1] = 1;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += +4.0 - I * 4.0;
        } else {
            bits[0] = 1;
            a += +4.0 + I * 4.0;
        }
    }

    unmap_16qam_type1(a, bits + 2);
}


// 256QAMアンマッピング(Type1)
void unmap_256qam_type1 (complex a, int *bits) {
    if (creal(a) > 0) {
        bits[1] = 0;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += -8.0 - I * 8.0;
        } else {
            bits[0] = 1;
            a += -8.0 + I * 8.0;
            a = conj(a);
        }
    } else {
        bits[1] = 1;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += +8.0 - I * 8.0;
            a = -conj(a);
        } else {
            bits[0] = 1;
            a += +8.0 + I * 8.0;
            a = -a;
        }
    }

    unmap_64qam_type1(a, bits + 2);
}


// 256QAMアンマッピング(Type2)
void unmap_256qam_type2 (complex a, int *bits) {
    if (creal(a) > 0) {
        bits[1] = 0;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += -8.0 - I * 8.0;
        } else {
            bits[0] = 1;
            a += -8.0 + I * 8.0;
        }
    } else {
        bits[1] = 1;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += +8.0 - I * 8.0;
        } else {
            bits[0] = 1;
            a += +8.0 + I * 8.0;
        }
    }

    unmap_64qam_type1(a, bits + 2);
}


// コンステレーションを初期化
void construct_constellation (complex *constellation, int m, int type) {
    int a, b, c, d, e, f, g, h;                             // ループカウンタ

    // QAM値のチェック
    if (m != 16 && m != 64 && m != 256) {
        printf("Invalid number of constellation.\n");
        exit(-1);
    }

    // タイプのチェック
    if (type != 1 && type != 2) {
        printf("Invalid argument for type (type = 1 or 2).\n");
        exit(-1);
    }

    for (a = 0; a < 2; a++) {
        for (b = 0; b < 2; b++) {
            for (c = 0; c < 2; c++) {
                for (d = 0; d < 2; d++) {
                    if (m == 16) {
                        if (type == 1) {
                            constellation[8*a + 4*b + 2*c + d] = map_16qam_type1(a, b, c, d);
                        } else {
                            constellation[8*a + 4*b + 2*c + d] = map_16qam_type2(a, b, c, d);
                        }
                    } else {
                        for (e = 0; e < 2; e++) {
                            for (f = 0; f < 2; f++) {
                                if (m == 64) {
                                    if (type == 1) {
                                        constellation[32*a + 16*b + 8*c + 4*d + 2*e + f] = map_64qam_type1(a, b, c, d, e, f);
                                    } else {
                                        constellation[32*a + 16*b + 8*c + 4*d + 2*e + f] = map_64qam_type2(a, b, c, d, e, f);
                                    }
                                } else {
                                    for (g = 0; g < 2; g++) {
                                        for (h = 0; h < 2; h++) {
                                            if (m == 256) {
                                                if (type == 1) {
                                                    constellation[128*a + 64*b + 32*c + 16*d + 8*e + 4*f + 2*g + h] = map_256qam_type1(a, b, c, d, e, f, g, h);
                                                } else {
                                                    constellation[128*a + 64*b + 32*c + 16*d + 8*e + 4*f + 2*g + h] = map_256qam_type2(a, b, c, d, e, f, g, h);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}


// トレリスシェイピングの変調
void qam_modulation (int *c, complex *a, int n, int m, int type) {
    // コンステレーション数のチェック
    if (m != 16 && m != 64 && m != 256) {
        printf("Invalid number of constellation.\n");
        exit(-1);
    }

    // コンステレーション数のチェック
    if (type != 1 && type != 2) {
        printf("Invalid argument for type (type = 1 or 2).\n");
        exit(-1);
    }

    for (int i = 0; i < n; i++) {
        if (m == 16) {
            if (type == 1) {
                a[i] = map_16qam_type1(c[4*i], c[4*i+1], c[4*i+2], c[4*i+3]);
            } else {
                a[i] = map_16qam_type2(c[4*i], c[4*i+1], c[4*i+2], c[4*i+3]);
            }
        } else if (m == 64) {
            if (type == 1) {
                a[i] = map_64qam_type1(c[6*i], c[6*i+1], c[6*i+2], c[6*i+3], c[6*i+4], c[6*i+5]);
            } else {
                a[i] = map_64qam_type2(c[6*i], c[6*i+1], c[6*i+2], c[6*i+3], c[6*i+4], c[6*i+5]);
            }
        } else if (m == 256) {
            if (type == 1) {
                a[i] = map_256qam_type1(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3], c[8*i+4], c[8*i+5], c[8*i+6], c[8*i+7]);
            } else {
                a[i] = map_256qam_type2(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3], c[8*i+4], c[8*i+5], c[8*i+6], c[8*i+7]);
            }
        }
    }
}


// トレリスシェイピングの変調
void qam_modulation_lsb (int *c, complex *a, int n, int m) {
    // コンステレーション数のチェック
    if (m != 16 && m != 64 && m != 256) {
        printf("Invalid number of constellation.\n");
        exit(-1);
    }

    for (int i = 0; i < n; i++) {
        if (m == 16) {
            a[i] = map_4_4_qam(c[4*i], c[4*i+1]);
        } else if (m == 64) {
            a[i] = map_16_4_qam(c[6*i], c[6*i+1], c[6*i+2], c[6*i+3]);
        } else if (m == 256) {
            a[i] = map_64_4_qam(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3], c[8*i+4], c[8*i+5]);
        }
    }
}


// トレリスシェイピングの変調
void qam_modulation_lsb2 (int *c, complex *a, int n, int m) {
    // コンステレーション数のチェック
    if (m != 64 && m != 256) {
        printf("Invalid number of constellation.\n");
        exit(-1);
    }

    for (int i = 0; i < n; i++) {
        if (m == 64) {
            a[i] = map_4_16_qam(c[6*i], c[6*i+1]);
        } else if (m == 256) {
            a[i] = map_16_16_qam(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3]);
        }
    }
}


// トレリスシェイピングの変調
void qam_modulation_lsb3 (int *c, complex *a, int n, int m) {
    // コンステレーション数のチェック
    if (m != 256) {
        printf("Invalid number of constellation.\n");
        exit(-1);
    }

    for (int i = 0; i < n; i++) {
        if (m == 256) {
            a[i] = map_4_64_qam(c[8*i], c[8*i+1]);
        }
    }
}


// トレリスシェイピングの復調
void qam_demodulation (complex *a, int *c, int n, int m, int type) {
    const int num_bit = (int)log2(m);             // ビット数
    static int *bits;                                   // ビット列
    static int memory_flag;                             // メモリ管理フラグ
    int i;                                              // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (m != 16 && m != 64 && m != 256) {
            printf("Invalid number of constellation.\n");
            exit(-1);
        }

        // タイプのチェック
        if (type != 1 && type != 2) {
            printf("Invalid argument for type (type = 1 or 2).\n");
            exit(-1);
        }

        // ビット列を生成する
        bits = (int *)malloc(num_bit * sizeof(int));

        // フラグを立てる
        memory_flag = 1;
    }

    // 復調
    for (i = 0; i < n; i++) {
        if (m == 16) {
            if (type == 1) {
                unmap_16qam_type1(a[i], c + i * num_bit);
            } else {
                unmap_16qam_type2(a[i], c + i * num_bit);
            }
        } else if (m == 64) {
            if (type == 1) {
                unmap_64qam_type1(a[i], c + i * num_bit);
            } else {
                unmap_64qam_type2(a[i], c + i * num_bit);
            }
        } else if (m == 256) {
            if (type == 1) {
                unmap_256qam_type1(a[i], c + i * num_bit);
            } else {
                unmap_256qam_type2(a[i], c + i * num_bit);
            }
        }
    }
}


// FFT
void fft (int n, fftw_complex *in, fftw_complex *out) {
    int i;                      // ループカウンタ

    // FFT
    fftw_plan plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // 正規化
    for (i = 0; i < n; i++) {
        out[i] /= sqrt((double)n);
    }

    // planを破棄
    fftw_destroy_plan(plan);
}


// IFFT
void ifft (int n, fftw_complex *in, fftw_complex *out) {
    int i;                      // ループカウンタ

    // IFFT
    fftw_plan plan = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // 正規化
    for (i = 0; i < n; i++) {
        out[i] /= sqrt((double)n);
    }

    // planを破棄
    fftw_destroy_plan(plan);
}


// オーバーサンプリング
void over_sampling (complex *x, complex *y, int j, int n) {
    int i;                      // ループカウンタ
    const int m = j*n;          // yの配列数

    // 0で埋める
    for (i = 0; i < m; i++) {
        y[i] = (complex)0;
    }

    // xの半分に分けて入れる
    for (i = 0; i < n/2; i++) {
        y[i] = x[n/2 + i];
        y[m - n/2 + i] = x[i];
    }
}


// ダウンサンプリング
void down_sampling (complex *y, complex *x, int j, int n) {
    int i;                      // ループカウンタ
    const int m = j*n;          // yの配列数

    // xの半分に分けて入れる
    for (i = 0; i < n/2; i++) {
        x[n/2 + i] = y[i];
        x[i] = y[m - n/2 + i];
    }
}


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
            x[i] = r_max / cabs(x[i]) * x[i];
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


// トレリスシェーピング
void trellis_shaping (int *c, complex *a, int n, int m, int type) {
    const int num_state = 4;                                                            // 状態数
    const int num_trans = 2;                                                            // 遷移数
    const int next_state[num_state][num_trans] = {{0,1}, {2,3}, {0,1}, {2,3}};          // 次状態
    const int output1[num_state][num_trans] = {{0,1}, {0,1}, {1,0}, {1,0}};             // 1ビット目の出力
    const int output2[num_state][num_trans] = {{0,1}, {1,0}, {1,0}, {0,1}};             // 2ビット目の出力
    const double infty = (double)n * 10000000;                                          // 無限
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    int tmp_index;                                                                      // 合成した信号
    static complex *constellation;                                                      // 変調信号のテーブル
    static complex **delta_table1;                                                      // δテーブル
    static double **delta_table2;                                                       // |δ^2|テーブル
    static node **nodes;                                                                // ノード
    static int memory_flag;                                                             // メモリ管理フラグ
    int sub, state, input;                                                              // ループカウンタ
    int i, j, l;                                                                        // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (m != 16 && m != 64 && m != 256) {
            printf("Invalid number of constellation.\n");
            exit(-1);
        }

        // メモリの確保
        constellation = (complex *)malloc(m * sizeof(complex));
        nodes = (node **)malloc((n+1) * sizeof(node *));
        delta_table1 = (complex **)malloc(m * sizeof(complex *));
        delta_table2 = (double **)malloc(m * sizeof(double *));

        for (i = 0; i <= n; i++) {
            nodes[i] = (node *)malloc(num_state * sizeof(node));
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = (int *)malloc(n * sizeof(int));
                nodes[i][j].autocor = (complex *)malloc(n * sizeof(complex));
            }
        }
        for (i = 0; i < m; i++) {
            delta_table1[i] = (complex *)malloc(m * sizeof(complex));
            delta_table2[i] = (double *)malloc(m * sizeof(double));
        }

        // 信号点配置を初期化
        construct_constellation(constellation, m, type);

        // δのテーブルを用意する
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {
                delta_table1[i][j] = constellation[i] * conj(constellation[j]);
                delta_table2[i][j] = pow(cabs(delta_table1[i][j]), 2.0);
            }
        }

        // フラグをチェック
        memory_flag = 1;
    }

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0.0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // サブキャリアの位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 16) {
                    tmp_index = ((c[4*sub] ^ output1[state][input]) << 3) + ((c[4*sub+1] ^ output2[state][input]) << 2) + (c[4*sub+2] << 1) + c[4*sub+3];
                } else if (m == 64) {
                    tmp_index = ((c[6*sub] ^ output1[state][input]) << 5) + ((c[6*sub+1] ^ output2[state][input]) << 4) + (c[6*sub+2] << 3) + (c[6*sub+3] << 2) + (c[6*sub+4] << 1) + c[6*sub+5];
                } else if (m == 256){
                    tmp_index = ((c[8*sub] ^ output1[state][input]) << 7) + ((c[8*sub+1] ^ output2[state][input]) << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + (c[8*sub+6] << 1) + c[8*sub+7];
                }

                // ブランチメトリックを求める
                branch_metric = 0;

                // 第2項
                for (l = 1; l <= i-2; l++) {
                    branch_metric = dadd(branch_metric, dmul(2.0, creal(cipr(nodes[i-1][state].autocor[l], delta_table1[tmp_index][nodes[i-1][state].a_index[i-1-l]]))));
                }

                // 第3項
                if (type == 2) {
                    for (l = 1; l <= i-1; l++) {
                        branch_metric =  dadd(branch_metric, delta_table2[tmp_index][nodes[i-1][state].a_index[i-1-l]]);
                    }
                }

                // パスの選択
                if (nodes[i-1][state].metric + branch_metric < nodes[i][next_state[state][input]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[state][input]].metric = nodes[i-1][state].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[state][input]].pre_state = state;
                    nodes[i][next_state[state][input]].a_index[i-1] = tmp_index;
                }
            }
        }

        // 次状態に情報を渡す
        for (state = 0; state < num_state; state++) {
            if (nodes[i][state].metric > infty - 1.0) continue;

            // 変調信号を更新する
            for (l = 0; l <= i-2; l++) {
                nodes[i][state].a_index[l] = nodes[i-1][nodes[i][state].pre_state].a_index[l];
            }

            // 自己相関を更新する
            for (l = 1; l <= i-2; l++) {
                nodes[i][state].autocor[l] = nodes[i-1][nodes[i][state].pre_state].autocor[l] + delta_table1[nodes[i][state].a_index[i-1]][nodes[i][state].a_index[i-1-l]];
            }
            nodes[i][state].autocor[i-1] = delta_table1[nodes[i][state].a_index[i-1]][nodes[i][state].a_index[0]];
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


// トレリスシェーピング
void trellis_shaping_caf (int *c, complex *a_caf, complex *a, int n, int m) {
    const int num_state = 4;                                                            // 状態数
    const int num_trans = 2;                                                            // 遷移数
    const int next_state[num_state][num_trans] = {{0,1}, {2,3}, {0,1}, {2,3}};          // 次状態
    const int output1[num_state][num_trans] = {{0,1}, {0,1}, {1,0}, {1,0}};             // 1ビット目の出力
    const int output2[num_state][num_trans] = {{0,1}, {1,0}, {1,0}, {0,1}};             // 2ビット目の出力
    const double infty = (double)n * 10000000;                                          // 無限
    int tmp_index;                                                                      // 合成した信号
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    static complex *constellation;                                                      // 変調信号のテーブル
    static node **nodes;                                                                // ノード
    static int memory_flag;                                                             // メモリ管理フラグ
    int sub, state, input;                                                              // ループカウンタ
    int i, j, l;                                                                        // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (m != 16 && m != 64 && m != 256) {
            printf("Invalid number of constellation.\n");
            exit(-1);
        }

        // メモリの確保
        constellation = (complex *)malloc(m * sizeof(complex));
        nodes = (node **)malloc((n+1) * sizeof(node *));

        for (i = 0; i <= n; i++) {
            nodes[i] = (node *)malloc(num_state * sizeof(node));
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = (int *)malloc(n * sizeof(int));
            }
        }

        // 信号点配置を初期化
        construct_constellation(constellation, m, 1);

        // フラグをチェック
        memory_flag = 1;
    }

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 16) {
                    tmp_index = (c[4*sub] << 3) + (c[4*sub+1] << 2) + ((c[4*sub+2] ^ output1[state][input]) << 1) + (c[4*sub+3] ^ output2[state][input]);
                } else if (m == 64) {
                    tmp_index = (c[6*sub] << 5) + (c[6*sub+1] << 4) + (c[6*sub+2] << 3) + (c[6*sub+3] << 2) + ((c[6*sub+4] ^ output1[state][input]) << 1) + (c[6*sub+5] ^ output2[state][input]);
                } else if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + ((c[8*sub+6] ^ output1[state][input]) << 1) + (c[8*sub+7] ^ output2[state][input]);
                }

                // ブランチメトリックを求める
                branch_metric = cdis(a_caf[sub], constellation[tmp_index]);

                // パスの選択
                if (nodes[i-1][state].metric + branch_metric < nodes[i][next_state[state][input]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[state][input]].metric = nodes[i-1][state].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[state][input]].pre_state = state;
                    nodes[i][next_state[state][input]].a_index[sub] = tmp_index;
                }
            }
        }

        // 次状態に情報を渡す
        for (state = 0; state < num_state; state++) {
            if (nodes[i][state].metric > infty - 1.0) continue;

            // 変調信号を更新する
            for (l = 0; l < sub; l++) {
                nodes[i][state].a_index[l] = nodes[i-1][nodes[i][state].pre_state].a_index[l];
            }
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


// トレリスシェーピング
void trellis_shaping_caf2 (int *c, complex *a_caf, complex *a, int n, int m) {
    const int num_state = 4;                                                            // 状態数
    const int num_trans = 2;                                                            // 遷移数
    const int next_state[num_state][num_trans] = {{0,1}, {2,3}, {0,1}, {2,3}};          // 次状態
    const int output1[num_state][num_trans] = {{0,1}, {0,1}, {1,0}, {1,0}};             // 1ビット目の出力
    const int output2[num_state][num_trans] = {{0,1}, {1,0}, {1,0}, {0,1}};             // 2ビット目の出力
    const double infty = (double)n * 10000000;                                          // 無限
    int tmp_index;                                                                      // 合成した信号
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    static complex *constellation;                                                      // 変調信号のテーブル
    static node **nodes;                                                                // ノード
    static int memory_flag;                                                             // メモリ管理フラグ
    int sub, state, input;                                                              // ループカウンタ
    int i, j, l;                                                                        // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (m != 64 && m != 256) {
            printf("Invalid number of constellation.\n");
            exit(-1);
        }

        // メモリの確保
        constellation = (complex *)malloc(m * sizeof(complex));
        nodes = (node **)malloc((n+1) * sizeof(node *));

        for (i = 0; i <= n; i++) {
            nodes[i] = (node *)malloc(num_state * sizeof(node));
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = (int *)malloc(n * sizeof(int));
            }
        }

        // 信号点配置を初期化
        construct_constellation(constellation, m, 1);

        // フラグをチェック
        memory_flag = 1;
    }

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 64) {
                    tmp_index = (c[6*sub] << 5) + (c[6*sub+1] << 4) + ((c[6*sub+2] ^ output1[state][input]) << 3) + ((c[6*sub+3] ^ output2[state][input]) << 2) + (c[6*sub+4] << 1) + c[6*sub+5];
                } else if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + ((c[8*sub+4] ^ output1[state][input]) << 3) + ((c[8*sub+5] ^ output2[state][input]) << 2) + (c[8*sub+6] << 1) + c[8*sub+7];
                }

                // ブランチメトリックを求める
                branch_metric = cdis(a_caf[sub], constellation[tmp_index]);

                // パスの選択
                if (nodes[i-1][state].metric + branch_metric < nodes[i][next_state[state][input]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[state][input]].metric = nodes[i-1][state].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[state][input]].pre_state = state;
                    nodes[i][next_state[state][input]].a_index[sub] = tmp_index;
                }
            }
        }

        // 次状態に情報を渡す
        for (state = 0; state < num_state; state++) {
            if (nodes[i][state].metric > infty - 1.0) continue;

            // 変調信号を更新する
            for (l = 0; l < sub; l++) {
                nodes[i][state].a_index[l] = nodes[i-1][nodes[i][state].pre_state].a_index[l];
            }
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }

    qam_demodulation(a, c, n, m, 1);

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 16) {
                    tmp_index = (c[4*sub] << 3) + (c[4*sub+1] << 2) + ((c[4*sub+2] ^ output1[state][input]) << 1) + (c[4*sub+3] ^ output2[state][input]);
                } else if (m == 64) {
                    tmp_index = (c[6*sub] << 5) + (c[6*sub+1] << 4) + (c[6*sub+2] << 3) + (c[6*sub+3] << 2) + ((c[6*sub+4] ^ output1[state][input]) << 1) + (c[6*sub+5] ^ output2[state][input]);
                } else if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + ((c[8*sub+6] ^ output1[state][input]) << 1) + (c[8*sub+7] ^ output2[state][input]);
                }

                // ブランチメトリックを求める
                branch_metric = cdis(a_caf[sub], constellation[tmp_index]);

                // パスの選択
                if (nodes[i-1][state].metric + branch_metric < nodes[i][next_state[state][input]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[state][input]].metric = nodes[i-1][state].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[state][input]].pre_state = state;
                    nodes[i][next_state[state][input]].a_index[sub] = tmp_index;
                }
            }
        }

        // 次状態に情報を渡す
        for (state = 0; state < num_state; state++) {
            if (nodes[i][state].metric > infty - 1.0) continue;

            // 変調信号を更新する
            for (l = 0; l < sub; l++) {
                nodes[i][state].a_index[l] = nodes[i-1][nodes[i][state].pre_state].a_index[l];
            }
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


// トレリスシェーピング
void trellis_shaping_caf3 (int *c, complex *a_caf, complex *a, int n, int m) {
    const int num_state = 4;                                                            // 状態数
    const int num_trans = 2;                                                            // 遷移数
    const int next_state[num_state][num_trans] = {{0,1}, {2,3}, {0,1}, {2,3}};          // 次状態
    const int output1[num_state][num_trans] = {{0,1}, {0,1}, {1,0}, {1,0}};             // 1ビット目の出力
    const int output2[num_state][num_trans] = {{0,1}, {1,0}, {1,0}, {0,1}};             // 2ビット目の出力
    const double infty = (double)n * 10000000;                                          // 無限
    int tmp_index;                                                                      // 合成した信号
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    static complex *constellation;                                                      // 変調信号のテーブル
    static node **nodes;                                                                // ノード
    static int memory_flag;                                                             // メモリ管理フラグ
    int sub, state, input;                                                              // ループカウンタ
    int i, j, l;                                                                        // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (m != 256) {
            printf("Invalid number of constellation.\n");
            exit(-1);
        }

        // メモリの確保
        constellation = (complex *)malloc(m * sizeof(complex));
        nodes = (node **)malloc((n+1) * sizeof(node *));

        for (i = 0; i <= n; i++) {
            nodes[i] = (node *)malloc(num_state * sizeof(node));
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = (int *)malloc(n * sizeof(int));
            }
        }

        // 信号点配置を初期化
        construct_constellation(constellation, m, 1);

        // フラグをチェック
        memory_flag = 1;
    }

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + ((c[8*sub+2] ^ output1[state][input]) << 5) + ((c[8*sub+3] ^ output2[state][input]) << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + (c[8*sub+6] << 1) + c[8*sub+7];
                }

                // ブランチメトリックを求める
                branch_metric = cdis(a_caf[sub], constellation[tmp_index]);

                // パスの選択
                if (nodes[i-1][state].metric + branch_metric < nodes[i][next_state[state][input]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[state][input]].metric = nodes[i-1][state].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[state][input]].pre_state = state;
                    nodes[i][next_state[state][input]].a_index[sub] = tmp_index;
                }
            }
        }

        // 次状態に情報を渡す
        for (state = 0; state < num_state; state++) {
            if (nodes[i][state].metric > infty - 1.0) continue;

            // 変調信号を更新する
            for (l = 0; l < sub; l++) {
                nodes[i][state].a_index[l] = nodes[i-1][nodes[i][state].pre_state].a_index[l];
            }
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }

    qam_demodulation(a, c, n, m, 1);

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + ((c[8*sub+4] ^ output1[state][input]) << 3) + ((c[8*sub+5] ^ output2[state][input]) << 2) + (c[8*sub+6] << 1) + c[8*sub+7];
                }

                // ブランチメトリックを求める
                branch_metric = cdis(a_caf[sub], constellation[tmp_index]);

                // パスの選択
                if (nodes[i-1][state].metric + branch_metric < nodes[i][next_state[state][input]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[state][input]].metric = nodes[i-1][state].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[state][input]].pre_state = state;
                    nodes[i][next_state[state][input]].a_index[sub] = tmp_index;
                }
            }
        }

        // 次状態に情報を渡す
        for (state = 0; state < num_state; state++) {
            if (nodes[i][state].metric > infty - 1.0) continue;

            // 変調信号を更新する
            for (l = 0; l < sub; l++) {
                nodes[i][state].a_index[l] = nodes[i-1][nodes[i][state].pre_state].a_index[l];
            }
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }

    qam_demodulation(a, c, n, m, 1);

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + ((c[8*sub+6] ^ output1[state][input]) << 1) + (c[8*sub+7] ^ output2[state][input]);
                }

                // ブランチメトリックを求める
                branch_metric = cdis(a_caf[sub], constellation[tmp_index]);

                // パスの選択
                if (nodes[i-1][state].metric + branch_metric < nodes[i][next_state[state][input]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[state][input]].metric = nodes[i-1][state].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[state][input]].pre_state = state;
                    nodes[i][next_state[state][input]].a_index[sub] = tmp_index;
                }
            }
        }

        // 次状態に情報を渡す
        for (state = 0; state < num_state; state++) {
            if (nodes[i][state].metric > infty - 1.0) continue;

            // 変調信号を更新する
            for (l = 0; l < sub; l++) {
                nodes[i][state].a_index[l] = nodes[i-1][nodes[i][state].pre_state].a_index[l];
            }
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


double calc_average_power (complex *a, int n) {
    int i;                             // ループカウンタ
    double average = 0;                // ピーク電力と平均電力

    for (i = 0; i < n; i++) {
        // 平均の計算
        average += pow(cabs(a[i]), 2.0) / (double)n;
    }

    // 比を返す
    return average;
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
