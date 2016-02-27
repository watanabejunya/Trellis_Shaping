#ifndef LIB_TS_H
#define LIB_TS_H

#include <complex.h>
#include <fftw3.h>

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif

/*
 * 変数
 */


/*
 * 構造体
 */

// ノードの定義
typedef struct {
    int *a_index;
    complex *autocor;
    double metric;
    int pre_state;
} node;

/*
 * 関数
 */

// ファイルオープン
FILE *fsopen(const char *mode, const char *format, ...);

// 0,1信号生成
void make_signal(int *x, int n);

// ガウス雑音
void gaussian_noise(complex *x, double sigma, int n);

// ビットエラーを数える
int count_be(int *x, int *y, int n);

// マッピング出力
void print_map(FILE *fp, fftw_complex *x, int n);

// 整数をコピー
void copy_int (int *x, int *x_copy, int n);

// 複素数をコピー
void copy_complex (complex *x, complex *x_copy, int n);

// インタリーバ
void xor_addition (int *x, int *y, int n);

//デマルチプレクサ(num_d → num_s + num_b)
void demultiplexer (int *d, int *s, int *b, int num_d, int num_s, int num_b, int n);

//マルチプレクサ(num_s + num_b → num_d)
void multiplexer (int *s, int *b, int *d, int num_s, int num_b, int num_d, int n);

// 畳み込み符号化(拘束長3)
void convolutional_encoding (int *u, int *c, int n);

// 畳み込み符号化(拘束長3)
void convolutional_encoding2 (int *u, int *c, int n);

// 畳み込み符号化(拘束長3)
void convolutional_encoding3 (int *u, int *c, int n);

// パリティ検査行列による復号(拘束長3)
void parity_check_decoding (int *c, int *u, int n);

// パリティ検査行列による復号(拘束長3)
void parity_check_decoding2 (int *c, int *u, int n);

// パリティ検査行列による復号(拘束長3)
void parity_check_decoding3 (int *c, int *u, int n);

// パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_encoding (int *u, int *c, int n);

// 並列パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_encoding2 (int *u, int *c, int n);

// 並列パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_encoding3 (int *u, int *c, int n);

// 下位ビットQPSKマッピング(16QAM対応)
complex map_4_4_qam (int bit1, int bit2);

// 下位ビットQPSKマッピング(64QAM対応)
complex map_4_16_qam (int bit1, int bit2);

// 下位ビットQPSKマッピング(64QAM対応)
complex map_4_64_qam (int bit1, int bit2);

// 16QAMマッピング(type1)
complex map_16qam_type1 (int bit1, int bit2, int bit3, int bit4);

// 16QAMマッピング(Type2)
complex map_16qam_type2 (int bit1, int bit2, int bit3, int bit4);

// 下位ビット16QAMマッピング(64QAM対応)
complex map_16_4_qam (int bit1, int bit2, int bit3, int bit4);

// 下位ビット16QAMマッピング(256QAM対応)
complex map_16_16_qam (int bit1, int bit2, int bit3, int bit4);

// 64QAMマッピング(Type1)
complex map_64qam_type1 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);

// 64QAMマッピング(Type2)
complex map_64qam_type2 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);

// 下位ビット64QAMマッピング(256QAM対応)
complex map_64_4_qam (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);

// 256QAMマッピング
complex map_256qam_type1 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8);

// 256QAMマッピング(Type2)
complex map_256qam_type2 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6, int bit7, int bit8);

// QPSKアンマッピング
void unmap_qpsk (complex a, int *bits);

// 16QAMアンマッピング(Type1)
void unmap_16qam_type1 (complex a, int *bits);

// 16QAMアンマッピング(Type2)
void unmap_16qam_type2 (complex a, int *bits);

// 64QAMアンマッピング(Type1)
void unmap_64qam_type1 (complex a, int *bits);

// 64QAMアンマッピング(Type2)
void unmap_64qam_type2 (complex a, int *bits);

// 256QAMアンマッピング(Type1)
void unmap_256qam_type1 (complex a, int *bits);

// 256QAMアンマッピング(Type2)
void unmap_256qam_type2 (complex a, int *bits);

// コンステレーションを初期化
void construct_constellation (complex *constellation, int m, int type);

void qam_modulation (int *c, complex *a, int n, int m, int type);

// トレリスシェイピングの変調
void qam_modulation_lsb (int *c, complex *a, int n, int m);

// トレリスシェイピングの変調
void qam_modulation_lsb2 (int *c, complex *a, int n, int m);

// トレリスシェイピングの変調
void qam_modulation_lsb3 (int *c, complex *a, int n, int m);

// トレリスシェイピングの復調
void qam_demodulation (complex *a, int *c, int n, int m, int type);

// FFT
void fft (int n, fftw_complex *in, fftw_complex *out);

// IFFT
void ifft (int n, fftw_complex *in, fftw_complex *out);

// オーバーサンプリング
void over_sampling (complex *x, fftw_complex *y, int j, int n);

// ダウンサンプリング
void down_sampling (fftw_complex *x, complex *y, int j, int n);

// クリッピング
void clipping (complex *x, int n, double r);

// 減衰を補償する
void offset_attenuation (complex *x, int n, double r);

// トレリスシェーピング
void trellis_shaping (int *c, complex *a, int n, int num_qam, int type);

// トレリスシェーピング
void trellis_shaping_caf (int *c, complex *a_caf, complex *a, int n, int m);

// トレリスシェーピング
void trellis_shaping_caf2 (int *c, complex *a_caf, complex *a, int n, int m);

// トレリスシェーピング
void trellis_shaping_caf3 (int *c, complex *a_caf, complex *a, int n, int m);

// 平均電力を求める
double calc_average_power(complex *a, int n);

// PAPR(dB)を求める
double calc_papr_db (complex *a, int n);

// 正規化瞬時電力のCCDF
double calc_normalized_ccdf (complex *a, int n, double threshold);

// PAPRの分布を求める
void count_papr_distribution (complex *x, double *pdf, int n, int min, int max, int num_index);

// 正規化瞬時電力の分布を求める
void count_normalized_distribution (complex *x, double *pdf, int n, int min, int max, int num_index);

// CCDFを計算する
double integrate_ccdf (double *pdf, double papr, int n, int min, int max, int num_index);

#endif
