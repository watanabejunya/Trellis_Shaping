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
void convolutional_coding (int *u, int *c, int n);

// パリティ検査行列による復号(拘束長3)
void parity_check_decoding (int *c, int *u, int n);

// パリティ検査行列による復号(拘束長3)
void parallel_parity_check_decoding (int *c, int *u, int n);

// パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_coding (int *u, int *c, int n);

// 並列パリティ検査行列の左逆行列による符号化(拘束長3)
void parallel_inverse_parity_check_coding (int *u, int *c, int n);

// 下位ビットZPSKマッピング(16QAM対応)
complex map_zpsk_lsb_extended ();

// 下位ビットQPSKマッピング(16QAM対応)
complex map_qpsk_lsb (int bit1, int bit2);

// 下位ビットQPSKマッピング(64QAM対応)
complex map_qpsk_lsb_extended (int bit1, int bit2);

// 16QAMマッピング(type1)
complex map_16qam_type1 (int bit1, int bit2, int bit3, int bit4);

// 16QAMマッピング(Type2)
complex map_16qam_type2 (int bit1, int bit2, int bit3, int bit4);

// 下位ビット16QAMマッピング(64QAM対応)
complex map_16qam_lsb (int bit1, int bit2, int bit3, int bit4);

// 下位ビット16QAMマッピング(256QAM対応)
complex map_16qam_lsb_extended (int bit1, int bit2, int bit3, int bit4);

// 64QAMマッピング(Type1)
complex map_64qam_type1 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);

// 64QAMマッピング(Type2)
complex map_64qam_type2 (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);

// 下位ビット64QAMマッピング(256QAM対応)
complex map_64qam_lsb (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6);

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
void construct_constellation (complex *constellation, int num_qam, int type);

void qam_modulation (int *c, complex *a, int n, int num_qam, int type);

// トレリスシェイピングの変調
void qam_modulation_lsb (int *c, complex *a, int n, int num_qam);

// トレリスシェイピングの変調
void qam_modulation_lsb_extended (int *c, complex *a, int n, int num_qam);

// トレリスシェイピングの復調
void qam_demodulation (complex *a, int *c, int n, int num_qam, int type);

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
void trellis_shaping_caf (int *c, complex *a_caf, complex *a, int n, int num_qam);

// トレリスシェーピング
void trellis_shaping_caf_parallel (int *c, complex *a_caf, complex *a, int n, int num_qam);

// 平均電力を求める
double calc_average_power(complex *a, int n);

// PAPR(dB)を求める
double calc_papr_db (complex *a, int n);

// 正規化瞬時電力のCCDF
double calc_normalized_ccdf (complex *a, int n, double threshold);



#endif
