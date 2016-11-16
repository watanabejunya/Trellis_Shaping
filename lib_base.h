#ifndef LIB_BASE_H
#define LIB_BASE_H


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


FILE *fsopen(const char *mode, const char *format, ...);
void make_signal(int *x, int n);
void joint_signal(complex *x, complex *y, int n, int m);
void gaussian_noise(complex *x, double sigma, int n);
int count_bit_error(int *x, int *y, int n);
void print_int(FILE *fp, int *x, int n);
void print_double(FILE *fp, double *x, int n);
void print_two_doubles(FILE *fp, double *x, double *y, int n);
void print_real_imag(FILE *fp, complex *x, int n);
void print_power(FILE *fp, complex *x, int n);
void print_absolute(FILE *fp, complex *x, int n);
void copy_int(int *x, int *x_copy, int n);
void copy_double(double *x, double *x_copy, int n);
void copy_complex(complex *x, complex *x_copy, int n);
void scale_complex(complex *x, int n, double scale);
void xor_addition(int *x, int *y, int n);
void random_interleave(int *x, int *y, int *p, int n);
int convert_binary_into_decimal(int *b, int n);
void convert_decimal_into_binary(int d, int *b, int n);
int check_power_of_2(int n);
int check_power_of_4(int n);
void convolutional_encoding(int *u, int *c, int n);
void parity_check_decoding(int *c, int *u, int n);
void inverse_parity_check_encoding(int *u, int *c, int n);
complex map_qpsk(int bit1, int bit2);
complex map_16qam_type1(int bit1, int bit2, int bit3, int bit4);
complex map_16qam_type2(int bit1, int bit2, int bit3, int bit4);
complex map_64qam_type1(int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);
complex map_64qam_type2(int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);
complex map_256qam_type1(int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8);
complex map_256qam_type2(int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8);
void unmap_qpsk(complex a, int *bits);
void unmap_16qam_type1(complex a, int *bits);
void unmap_16qam_type2(complex a, int *bits);
void unmap_64qam_type1(complex a, int *bits);
void unmap_64qam_type2(complex a, int *bits);
void unmap_256qam_type1(complex a, int *bits);
void unmap_256qam_type2(complex a, int *bits);
void construct_constellation(complex *constellation, int m, int type);
void qam_modulation(int *c, complex *a, int n, int m);
void qam_demodulation(complex *a, int *c, int n, int m);
void fft(int n, fftw_complex *in, fftw_complex *out);
void ifft(int n, fftw_complex *in, fftw_complex *out);
void over_sampling(complex *x, fftw_complex *y, int j, int n);
void down_sampling(fftw_complex *x, complex *y, int j, int n);
double calc_average_power(complex *a, int n);
double calc_average_amplitude(complex *a, int n);
void calc_amp_phase_pmf(complex *x, double *amp_pmf, double *phase_pmf, int n, int m);


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


// 0,1信号生成
void make_signal(int *x, int n) {
    for (int i = 0; i < n; i++) {
        x[i] = random() % 2;
    }
}


void joint_signal(complex *x, complex *y, int n, int m)
{
    for (int i = 0; i < m; i++) {
        x[n + i] = y[i];
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

        v = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2) + I * sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

        x[i] += v * sigma;
    }
}


// ビットエラーを数える
int count_bit_error(int *x, int *y, int n)
{
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (x[i] != y[i]) {
            count++;
        }
    }

    return count;
}


void print_int(FILE *fp, int *x, int n)
{
    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%d %d\n", i, x[i]);
        }
    }
}


void print_double(FILE *fp, double *x, int n)
{
    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%d %lf\n", i, x[i]);
        }
    }
}


void print_two_doubles(FILE *fp, double *x, double *y, int n)
{
    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%lf %lf\n", x[i], y[i]);
        }
    }
}


// マッピング出力
void print_real_imag(FILE *fp, complex *x, int n)
{
    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%lf %lf\n", creal(x[i]), cimag(x[i]));
        }
    }
}


// マッピング出力
void print_power(FILE *fp, complex *x, int n)
{
    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%d %lf\n", i, pow(cabs(x[i]), 2.0));
        }
    }
}


// マッピング出力
void print_absolute(FILE *fp, complex *x, int n)
{
    if (fp != NULL) {
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%d %lf\n", i, cabs(x[i]));
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


// 少数をコピー
void copy_double (double *x, double *x_copy, int n) {
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


void scale_complex(complex *x, int n, double scale)
{
    int i;          // ループカウンタ

    for (i = 0; i < n; i++) {
        x[i] *= scale;
    }
}


// XOR加算
void xor_addition (int *x, int *y, int n) {
    int i;

    for (i = 0; i < n; i++) {
        x[i] = x[i] ^ y[i];
    }
}


// インターリーブ
void random_interleave(int *x, int *y, int *p, int n)
{
    int index1, index2;
    int tmp;
    int i;

    for (i = 0; i < n; i++) {
        p[i] = i;
    }

    for (i = 0; i < n * 10; i++) {
        index1 = random() % n;
        index2 = random() % n;

        tmp = p[index1];
        p[index1] = p[index2];
        p[index2] = tmp;
    }

    for (i = 0; i < n; i++) {
        y[i] = x[p[i]];
    }
}


// インターリーブ
void random_deinterleave(int *x, int *y, int *p, int n)
{
    int i;

    for (i = 0; i < n; i++) {
        y[p[i]] = x[i];
    }
}


int convert_binary_into_decimal(int *b, int n)
{
    int d = 0;              // 10進数
    int i;                  // ループカウンタ

    for (i = 0; i < n; i++) {
        if (b[i] != 0 && b[i] != 1) {
            fprintf(stderr, "cannot convert non-binary number into decimal one.\n");
            exit(-1);
        }

        d += b[i] * (int)pow(2.0, n - i - 1);
    }

    return d;
}


void convert_decimal_into_binary(int d, int *b, int n)
{
    int i;

    for (i = 0; i < n; i++) {
        b[n - 1 - i] = d & 1;
        d >>= 1;
    }
}


int check_power_of_2(int n)
{
    if (n > 0 && (n & (n-1)) == 0) {
        return 1;
    }
    return 0;
}


int check_power_of_4(int n)
{
    if (n > 0) {
        int m = (int)sqrt(n);
        if (n == m * m) {
            return check_power_of_2(m);
        }
    }
    return 0;
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


// 下位ビットqpskマッピング(16qam対応)
complex map_qpsk (int bit1, int bit2) {
    complex a;                      // 変調信号
    a = (1.0 - 2.0 * bit2) + I * (1.0 - 2.0 * bit1);
    return a;
}


// 16qamマッピング(type1)
complex map_16qam_type1 (int bit1, int bit2, int bit3, int bit4) {
    complex a;                      // 変調信号

    a = (3.0 - 2.0*bit4) + I * (3.0 - 2.0*bit3);
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }

    return a;
}


// 16qamマッピング(type2)
complex map_16qam_type2 (int bit1, int bit2, int bit3, int bit4) {
    complex a;                      // 変調信号

    a = (2.0 - 4.0 * bit2) + I * (2.0 - 4.0 * bit1);
    a += (1.0 - 2.0 * bit4) + I * (1.0 - 2.0 * bit3);

    return a;
}


// 64qamマッピング(type1)
complex map_64qam_type1 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6) {
    complex a;                      // 変調信号

    a = map_16qam_type1(bit3, bit4, bit5, bit6);
    a += 4.0 + I * 4.0;
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }

    return a;
}


// 64qamマッピング(type2)
complex map_64qam_type2 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6) {
    complex a;                      // 変調信号

    a = map_16qam_type1(bit3, bit4, bit5, bit6);
    a += (4.0 - 8.0 * bit2) + I * (4.0 - 8.0 * bit1);

    return a;
}


// 256qamマッピング
complex map_256qam_type1 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8) {
    complex a;                      // 変調信号

    a = map_64qam_type1(bit3, bit4, bit5, bit6, bit7, bit8);
    a += 8.0 + I * 8.0;
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }

    return a;
}


// 256qamマッピング(type2)
complex map_256qam_type2 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8) {
    complex a;                      // 変調信号

    a = map_64qam_type1(bit3, bit4, bit5, bit6, bit7, bit8);
    a += (8.0 - 16.0 * bit2) + I * (8.0 - 16.0 * bit1);

    return a;
}


// qpskアンマッピング
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


// 16qamアンマッピング(type1)
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


// 16qamアンマッピング(type2)
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


// 64qamアンマッピング(type1)
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


// 64qamアンマッピング(type2)
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


// 256qamアンマッピング(type1)
void unmap_256qam_type1 (complex a, int *bits) {
    if (creal(a) > 0) {
        bits[1] = 0;
        if (cimag(a) > 0) {
            bits[0] = 0;
            a += -8.0 - I* 8.0;
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


// 256qamアンマッピング(type2)
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


// トレリスシェイピングの変調
void qam_modulation (int *c, complex *a, int n, int m)
{
    // コンステレーション数のチェック
    if (check_power_of_4(m) == 0 || m > 256) {
        printf("invalid number of constellation.\n");
        exit(-1);
    }

    for (int i = 0; i < n; i++) {
        if (m == 4) {
            a[i] = map_qpsk(c[2*i], c[2*i+1]);
        } else if (m == 16) {
            a[i] = map_16qam_type1(c[4*i], c[4*i+1], c[4*i+2], c[4*i+3]);
        } else if (m == 64) {
            a[i] = map_64qam_type1(c[6*i], c[6*i+1], c[6*i+2], c[6*i+3], c[6*i+4], c[6*i+5]);
        } else if (m == 256) {
            a[i] = map_256qam_type1(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3], c[8*i+4], c[8*i+5], c[8*i+6], c[8*i+7]);
        }
    }
}


// QAM復調
void qam_demodulation(complex *a, int *c, int n, int m)
{
    const int num_bit = (int)log2(m);             // ビット数
    int i;                                        // ループカウンタ

    // コンステレーション数のチェック
    if (m != 16 && m != 64 && m != 256) {
        printf("invalid number of constellation.\n");
        exit(-1);
    }

    // 復調
    for (i = 0; i < n; i++) {
        if (m == 16) {
            unmap_16qam_type1(a[i], c + i * num_bit);
        } else if (m == 64) {
            unmap_64qam_type1(a[i], c + i * num_bit);
        } else if (m == 256) {
            unmap_256qam_type1(a[i], c + i * num_bit);
        }
    }
}


// fft
void fft (int n, fftw_complex *in, fftw_complex *out) {
    int i;                      // ループカウンタ

    // fft
    fftw_plan plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // 正規化
    for (i = 0; i < n; i++) {
        out[i] /= sqrt((double)n);
    }

    // planを破棄
    fftw_destroy_plan(plan);
}


// ifft
void ifft (int n, fftw_complex *in, fftw_complex *out) {
    int i;                      // ループカウンタ

    // ifft
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


double calc_average_power (complex *a, int n)
{
    int i;                             // ループカウンタ
    double average = 0;                // ピーク電力と平均電力

    for (i = 0; i < n; i++) {
        // 平均の計算
        average += pow(cabs(a[i]), 2.0) / (double)n;
    }

    // 比を返す
    return average;
}


double calc_average_amplitude (complex *a, int n)
{
    int i;                             // ループカウンタ
    double average = 0;                // ピーク電力と平均電力

    for (i = 0; i < n; i++) {
        // 平均の計算
        average += cabs(a[i]) / (double)n;
    }

    // 比を返す
    return average;
}


void calc_amp_phase_pmf (complex *x, double *amp_pmf, double *phase_pmf, int n, int m)
{
    double *amp, *phase;
    double amp_max = 0;
    double amp_min = 100000;
    int tmp_index;
    int i;

    amp = (double *) malloc(n * sizeof(double));
    phase = (double *) malloc(n * sizeof(double));

    // 振幅と位相を計算
    for (i = 0; i < n; i++) {
        amp[i] = cabs(x[i]);
        phase[i] = carg(x[i]);
    }

    // 振幅の最小と最大を計算
    for (i = 0; i < n; i++) {
        if (amp[i] > amp_max) {
            amp_max = amp[i];
        } else if (amp[i] < amp_min) {
            amp_min = amp[i];
        }
    }

    // 出力用の変数を初期化
    for (i = 0; i < m; i++) {
        amp_pmf[i] = 0;
        phase_pmf[i] = 0;
    }

    // 分布を作成
    for (i = 0; i < n; i++) {
        tmp_index = (int) round((amp[i] - amp_min) / (amp_max - amp_min) * (m - 1));
        amp_pmf[tmp_index] += 1.0 / m;

        tmp_index = (int) round((phase[i] + M_PI) / (2.0 * M_PI) * (m - 1));
        phase_pmf[tmp_index] += 1.0 / m;
    }

    free(amp);
    free(phase);
}


#endif
