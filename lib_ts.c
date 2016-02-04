#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <stdarg.h>
#include "lib_ts.h"

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif


// 変調信号のテーブル
complex *constellation;


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
void print_map(FILE *fp, fftw_complex *x, int n) {
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
void interleaver (int *x, int *y, int n) {
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
void convolutional_coding (int *info, int *code, int info_length) {
    int i;                              // ループカウンタ
    int register1, register2;           // シフトレジスタ

    // シフトレジスタを初期化
    register1 = 0;
    register2 = 0;

    // 畳み込み
    for (i = 0; i < info_length; i++) {
        code[2*i] = info[i] ^ register2;
        code[2*i + 1] = info[i] ^ register1 ^ register2;

        // レジスタをシフト
        register2 = register1;
        register1 = info[i];
    }
}


// パリティ検査行列による復号(拘束長3)
void parity_check_decoding (int *code, int *info, int info_length) {
    int i;                                                      // ループカウンタ
    int register11, register12, register21, register22;         // シフトレジスタ

    // シフトレジスタを初期化
    register11 = 0;
    register12 = 0;
    register21 = 0;
    register22 = 0;

    // 畳み込み
    for (i = 0; i < info_length; i++) {
        info[i] = code[2*i] ^ register11 ^ register12 ^ code[2*i + 1] ^ register22;

        // レジスタをシフト
        register12 = register11;
        register11 = code[2*i];
        register22 = register21;
        register21 = code[2*i + 1];
    }
}


// パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_coding (int *info, int *code, int info_length) {
    int i;                              // ループカウンタ
    int register1;                      // シフトレジスタ

    // シフトレジスタを初期化
    register1 = 0;

    // 畳み込み
    for (i = 0; i < info_length; i++) {
        code[2*i] = register1;
        code[2*i + 1] = info[i] ^ register1;

        // レジスタをシフト
        register1 = info[i];
    }
}


// QPSKマッピング(16QAM拡張)
complex map_qpsk_extend (int bit1, int bit2) {
    complex a;                      // 変調信号

    a = (2.0 - 4.0 * bit2) + I*(2.0 - 4.0 * bit1);

    return a;
}


// QPSKマッピング(64QAM拡張)
complex map_qpsk_extend2 (int bit1, int bit2) {
    complex a;                      // 変調信号

    a = (4.0 - 8.0 * bit2) + I*(4.0 - 8.0 * bit1);

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


// 16QAMマッピング(64QAM拡張)
complex map_16qam_extend (int bit1, int bit2, int bit3, int bit4) {
    complex a;                      // 変調信号

    a = map_qpsk_extend(bit3, bit4);
    a += 4.0 + I*4.0;
    if (bit1 == 1) {
        a = conj(a);
    }
    if (bit2 == 1) {
        a = -conj(a);
    }

    return a;
}


// 16QAMマッピング(256QAM拡張)
complex map_16qam_extend2 (int bit1, int bit2, int bit3, int bit4) {
    complex a;                      // 変調信号

    a = map_qpsk_extend(bit3, bit4);
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


// 64QAMマッピング(256QAM拡張)
complex map_64qam_extend (int bit1, int bit2, int bit3, int bit4, int bit5, int bit6) {
    complex a;                      // 変調信号

    a = map_16qam_extend(bit3, bit4, bit5, bit6);
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
void construct_constellation (int num_qam, int type) {
    int h, i, j, k, l, m, n, o;                             // ループカウンタ

    // メモリの確保
    constellation = (complex *)malloc(num_qam * sizeof(complex));

    // QAM値のチェック
    if (num_qam != 16 && num_qam != 64 && num_qam != 256) {
        printf("Invalid number of constellation.\n");
        exit(-1);
    }

    // タイプのチェック
    if (type != 1 && type != 2) {
        printf("Invalid argument for type (type = 1 or 2).\n");
        exit(-1);
    }

    for (h = 0; h < 2; h++) {
        for (i = 0; i < 2; i++) {
            for (j = 0; j < 2; j++) {
                for (k = 0; k < 2; k++) {
                    if (num_qam == 16) {
                        if (type == 1) {
                            constellation[8*h + 4*i + 2*j + k] = map_16qam_type1(h, i, j, k);
                        } else {
                            constellation[8*h + 4*i + 2*j + k] = map_16qam_type2(h, i, j, k);
                        }
                    } else {
                        for (l = 0; l < 2; l++) {
                            for (m = 0; m < 2; m++) {
                                if (num_qam == 64) {
                                    if (type == 1) {
                                        constellation[32*h + 16*i + 8*j + 4*k + 2*l + m] = map_64qam_type1(h, i, j, k, l, m);
                                    } else {
                                        constellation[32*h + 16*i + 8*j + 4*k + 2*l + m] = map_64qam_type2(h, i, j, k, l, m);
                                    }
                                } else {
                                    for (n = 0; n < 2; n++) {
                                        for (o = 0; o < 2; o++) {
                                            if (num_qam == 256) {
                                                if (type == 1) {
                                                    constellation[128*h + 64*i + 32*j + 16*k + 8*l + 4*m + 2*n + o] = map_256qam_type1(h, i, j, k, l, m, n, o);
                                                } else {
                                                    constellation[128*h + 64*i + 32*j + 16*k + 8*l + 4*m + 2*n + o] = map_256qam_type2(h, i, j, k, l, m, n, o);
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
void qam_modulation (int *c, complex *a, int n, int num_qam, int type) {
    // コンステレーション数のチェック
    if (num_qam != 16 && num_qam != 64 && num_qam != 256) {
        printf("Invalid number of constellation.\n");
        exit(-1);
    }

    // コンステレーション数のチェック
    if (type != 1 && type != 2) {
        printf("Invalid argument for type (type = 1 or 2).\n");
        exit(-1);
    }

    for (int i = 0; i < n; i++) {
        if (num_qam == 16) {
            if (type == 1) {
                a[i] = map_16qam_type1(c[4*i], c[4*i+1], c[4*i+2], c[4*i+3]);
            } else {
                a[i] = map_16qam_type2(c[4*i], c[4*i+1], c[4*i+2], c[4*i+3]);
            }
        } else if (num_qam == 64) {
            if (type == 1) {
                a[i] = map_64qam_type1(c[6*i], c[6*i+1], c[6*i+2], c[6*i+3], c[6*i+4], c[6*i+5]);
            } else {
                a[i] = map_64qam_type2(c[6*i], c[6*i+1], c[6*i+2], c[6*i+3], c[6*i+4], c[6*i+5]);
            }
        } else if (num_qam == 256) {
            if (type == 1) {
                a[i] = map_256qam_type1(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3], c[8*i+4], c[8*i+5], c[8*i+6], c[8*i+7]);
            } else {
                a[i] = map_256qam_type2(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3], c[8*i+4], c[8*i+5], c[8*i+6], c[8*i+7]);
            }
        }
    }
}


// トレリスシェイピングの変調
void qam_modulation_extend (int *c, complex *a, int n, int num_qam) {
    // コンステレーション数のチェック
    if (num_qam != 16 && num_qam != 64 && num_qam != 256) {
        printf("Invalid number of constellation.\n");
        exit(-1);
    }

    for (int i = 0; i < n; i++) {
        if (num_qam == 16) {
            a[i] = map_qpsk_extend(c[4*i], c[4*i+1]);
        } else if (num_qam == 64) {
            a[i] = map_16qam_extend(c[6*i], c[6*i+1], c[6*i+2], c[6*i+3]);
        } else if (num_qam == 256) {
            a[i] = map_64qam_extend(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3], c[8*i+4], c[8*i+5]);
        }
    }
}


// トレリスシェイピングの変調
void qam_modulation_extend2 (int *c, complex *a, int n, int num_qam) {
    // コンステレーション数のチェック
    if (num_qam != 64 && num_qam != 256) {
        printf("Invalid number of constellation.\n");
        exit(-1);
    }

    for (int i = 0; i < n; i++) {
        if (num_qam == 64) {
            a[i] = map_qpsk_extend2(c[6*i], c[6*i+1]);
        } else if (num_qam == 256) {
            a[i] = map_16qam_extend2(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3]);
        }
    }
}


// トレリスシェイピングの復調
void qam_demodulation (complex *a, int *c, int n, int num_qam, int type) {
    const int num_bit = (int)log2(num_qam);             // ビット数
    static int *bits;                                   // ビット列
    static int memory_flag;                             // メモリ管理フラグ
    int i;                                              // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (num_qam != 16 && num_qam != 64 && num_qam != 256) {
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
        if (num_qam == 16) {
            if (type == 1) {
                unmap_16qam_type1(a[i], c + i * num_bit);
            } else {
                unmap_16qam_type2(a[i], c + i * num_bit);
            }
        } else if (num_qam == 64) {
            if (type == 1) {
                unmap_64qam_type1(a[i], c + i * num_bit);
            } else {
                unmap_64qam_type2(a[i], c + i * num_bit);
            }
        } else if (num_qam == 256) {
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
void over_sampling (complex *x, fftw_complex *y, int j, int n) {
    int i;                      // ループカウンタ
    const int m = j*n;          // yの配列数

    // 0で埋める
    for (i = 0; i < m; i++) {
        y[i] = (complex)0;
    }

    // xの半分に分けて入れる
    for (i = 0; i < n/2; i++) {
        y[i] = x[n/2 + i];
        y[m-1 - i] = x[i];
    }
}


// ダウンサンプリング
void down_sampling (fftw_complex *x, complex *y, int j, int n) {
    int i;                      // ループカウンタ
    const int m = j*n;          // yの配列数

    // xの半分に分けて入れる
    for (i = 0; i < n/2; i++) {
        y[i] = x[m-1 - i];
        y[n/2 + i] = x[i];
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
void trellis_shaping (int *c, complex *a, int n, int num_qam) {
    const int num_state = 4;                                                            // 状態数
    const int next_state[4][2] = {{0,1}, {2,3}, {0,1}, {2,3}};                          // 次状態
    const int output1[4][2] = {{0,1}, {0,1}, {1,0}, {1,0}};                             // 1ビット目の出力
    const int output2[4][2] = {{0,1}, {1,0}, {1,0}, {0,1}};                             // 2ビット目の出力
    const double infty = (double)n * 1000000000;                                        // 無限
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    static complex **delta_table1, *delta_table1_base;                                  // δテーブル
    static complex **delta_table2, *delta_table2_base;                                  // |δ^2|テーブル
    static node **nodes, *nodes_base;                                                   // ノード
    static int *a_index_base;                                                           // 変調信号のメモリ
    static complex *autocor_base;                                                       // 相関関数のメモリ
    static int memory_flag;                                                             // メモリ管理フラグ
    int i, j, k, m;                                                                     // ループカウンタ


    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (num_qam != 16 && num_qam != 64 && num_qam != 256) {
            printf("Invalid number of constellation.\n");
            exit(-1);
        }

        // メモリの確保
        nodes = (node **)malloc((n+1) * sizeof(node *));
        nodes_base = (node *)malloc((n+1) * num_state * sizeof(node));
        a_index_base = (int *)malloc((n+1) * num_state * n * sizeof(int));
        autocor_base = (complex *)malloc((n+1) * num_state * n * sizeof(complex));
        delta_table1 = (complex **)malloc(num_qam * sizeof(complex *));
        delta_table2 = (complex **)malloc(num_qam * sizeof(complex *));
        delta_table1_base = (complex *)malloc(num_qam * num_qam * sizeof(complex));
        delta_table2_base = (complex *)malloc(num_qam * num_qam * sizeof(complex));

        for (i = 0; i <= n; i++) {
            nodes[i] = nodes_base + i * num_state;
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = a_index_base + (i * num_state + j) * n;
                nodes[i][j].autocor = autocor_base + (i * num_state + j) * n;
            }
        }
        for (i = 0; i < num_qam; i++) {
            delta_table1[i] = delta_table1_base + i * num_qam;
            delta_table2[i] = delta_table2_base + i * num_qam;
        }

        // δのテーブルを用意する
        for (i = 0; i < num_qam; i++) {
            for (j = 0; j < num_qam; j++) {
                delta_table1[i][j] = constellation[i] * conj(constellation[j]);
                delta_table2[i][j] = pow(cabs(delta_table1[i][j]), 2.0);
            }
        }

        // フラグをチェック
        memory_flag = 1;
    }

    // 最初のノードメトリックを初期化
    nodes[0][0].metric = 0;
    for (j = 1; j < num_state; j++) {
        nodes[0][j].metric = infty;
    }

    // 適当に埋めておく
    nodes[1][2].a_index[0] = 0;
    nodes[1][3].a_index[0] = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // ノードメトリックを初期化
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }

        // ノードメトリックの計算
        for (j = 0; j < num_state; j++) {
            for (k = 0; k < 2; k++) {
                // 仮の変調信号を求める
                if (num_qam == 16) {
                    nodes[i-1][j].a_index[i-1] = 8*(c[4*(i-1)] ^ output1[j][k]) + 4*(c[4*(i-1)+1] ^ output2[j][k]) + 2*c[4*(i-1)+3] + c[4*(i-1)+4];
                } else if (num_qam == 64) {
                    nodes[i-1][j].a_index[i-1] = 32*(c[6*(i-1)] ^ output1[j][k]) + 16*(c[6*(i-1)+1] ^ output2[j][k]) + 8*c[6*(i-1)+2] + 4*c[6*(i-1)+3] + 2*c[6*(i-1)+4] + c[6*(i-1)+5];
                } else if (num_qam == 256){
                    nodes[i-1][j].a_index[i-1] = 128*(c[8*(i-1)] ^ output1[j][k]) + 64*(c[8*(i-1)+1] ^ output2[j][k]) + 32*c[8*(i-1)+2] + 16*c[8*(i-1)+3] + 8*c[8*(i-1)+4] + 4*c[8*(i-1)+5] + 2*c[8*(i-1)+6] + c[8*(i-1)+7];
                }

                // ブランチメトリックを求める
                branch_metric = 0;

                // 第2項
                for (m = 1; m <= i-2; m++) {
                    branch_metric += 2.0 * creal(conj(nodes[i-1][j].autocor[m]) * delta_table1[nodes[i-1][j].a_index[i-1]][nodes[i-1][j].a_index[i-1-m]]);

                }

                // 第3項
                for (m = 1; m <= i-1; m++) {
                    branch_metric +=  delta_table2[nodes[i-1][j].a_index[i-1]][nodes[i-1][j].a_index[i-1-m]];
                }

                // パスの選択
                if (nodes[i-1][j].metric + branch_metric < nodes[i][next_state[j][k]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[j][k]].metric = nodes[i-1][j].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[j][k]].pre_state = j;
                    nodes[i][next_state[j][k]].a_index[i-1] = nodes[i-1][j].a_index[i-1];
                }
            }
        }

        // 次状態に情報を渡す
        for (j = 0; j < num_state; j++) {
            // 変調信号を更新する
            for (m = 0; m <= i-2; m++) {
                nodes[i][j].a_index[m] = nodes[i-1][nodes[i][j].pre_state].a_index[m];
            }

            // 自己相関関数を更新する
            for (m = 1; m <= i-2; m++) {
                nodes[i][j].autocor[m] = nodes[i-1][nodes[i][j].pre_state].autocor[m] + delta_table1[nodes[i][j].a_index[i-1]][nodes[i][j].a_index[i-1-m]];
            }
            nodes[i][j].autocor[i-1] = delta_table1[nodes[i][j].a_index[i-1]][nodes[i][j].a_index[0]];
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (j = 0; j < num_state; j++) {
        if (nodes[n][j].metric < nodes[n][likely_state].metric) {
            likely_state = j;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


// トレリスシェーピング
void trellis_shaping_caf (int *c, complex *a_caf, complex *a, int n, int num_qam) {
    const int num_state = 4;                                                        // 状態数
    const int next_state[4][2] = {{0,1}, {2,3}, {0,1}, {2,3}};                      // 次状態
    const int output1[4][2] = {{0,1}, {0,1}, {1,0}, {1,0}};                         // 1ビット目の出力
    const int output2[4][2] = {{0,1}, {1,0}, {1,0}, {0,1}};                         // 2ビット目の出力
    const double infty = (double)n * 1000000000;                                    // 無限
    double branch_metric;                                                           // ブランチメトリック
    int likely_state;                                                               // トレースバックの最尤状態
    static node **nodes, *nodes_base;                                               // ノード
    static int *a_index_base;                                                       // 変調信号のメモリ
    static complex *autocor_base;                                                   // 相関関数のメモリ
    static int memory_flag;                                                         // メモリ管理フラグ
    int i, j, k, m;                                                                 // ループカウンタ

    // メモリの確保
    if (memory_flag == 0) {
        nodes = (node **)malloc((n+1) * sizeof(node *));
        nodes_base = (node *)malloc((n+1) * num_state * sizeof(node));
        a_index_base = (int *)malloc((n+1) * num_state * n * sizeof(int));
        autocor_base = (complex *)malloc((n+1) * num_state * n * sizeof(complex));

        for (i = 0; i <= n; i++) {
            nodes[i] = nodes_base + i * num_state;
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = a_index_base + (i * num_state + j) * n;
                nodes[i][j].autocor = autocor_base + (i * num_state + j) * n;
            }
        }

        // フラグをチェック
        memory_flag = 1;
    }

    // 最初のノードメトリックを初期化
    nodes[0][0].metric = 0;
    for (j = 1; j < num_state; j++) {
        nodes[0][j].metric = infty;
    }

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // ノードメトリックを初期化
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }

        // ノードメトリックの計算
        for (j = 0; j < num_state; j++) {
            for (k = 0; k < 2; k++) {
                // 仮の変調信号を求める
                if (num_qam == 16) {
                    nodes[i-1][j].a_index[i-1] = 8*c[4*(i-1)] + 4*c[4*(i-1)+1] + 2*(c[4*(i-1)+2] ^ output1[j][k]) + (c[4*(i-1)+3] ^ output2[j][k]);
                } else if (num_qam == 64) {
                    nodes[i-1][j].a_index[i-1] = 32*c[6*(i-1)] + 16*c[6*(i-1)+1] + 8*c[6*(i-1)+2] + 4*c[6*(i-1)+3] + 2*(c[6*(i-1)+4] ^ output1[j][k]) + (c[6*(i-1)+5] ^ output2[j][k]);
                } else if (num_qam == 256) {
                    nodes[i-1][j].a_index[i-1] = 128*c[8*(i-1)] + 64*c[8*(i-1)+1] + 32*c[8*(i-1)+2] + 16*c[8*(i-1)+3] + 8*c[8*(i-1)+4] + 4*c[8*(i-1)+5] + 2*(c[8*(i-1)+6] ^ output1[j][k]) + (c[8*(i-1)+7] ^ output2[j][k]);
                }

                // ブランチメトリックを求める
                if (i > 1) {
                    branch_metric = pow(cabs(a_caf[i-1] - constellation[nodes[i-1][j].a_index[i-1]]), 2.0);
                }

                // パスの選択
                if (nodes[i-1][j].metric + branch_metric < nodes[i][next_state[j][k]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[j][k]].metric = nodes[i-1][j].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[j][k]].pre_state = j;
                    nodes[i][next_state[j][k]].a_index[i-1] = nodes[i-1][j].a_index[i-1];
                }
            }
        }

        // 次状態に情報を渡す
        for (j = 0; j < num_state; j++) {
            // 変調信号を更新する
            for (m = 0; m <= i-2; m++) {
                nodes[i][j].a_index[m] = nodes[i-1][nodes[i][j].pre_state].a_index[m];
            }
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (j = 0; j < num_state; j++) {
        if (nodes[n][j].metric < nodes[n][likely_state].metric) {
            likely_state = j;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


// トレリスシェーピング
void trellis_shaping_caf2 (int *c, complex *a_caf, complex *a, int n, int num_qam) {
    const int num_state = 4;                                                        // 状態数
    const int num_trans = 2;                                                        // 遷移数
    const int next_state[16][4] = {{0,1,4,5}, {2,3,6,7}, {0,1,4,5}, {2,3,6,7}, {8,9,12,13}, {10,11,14,15}, {8,9,12,13}, {10,11,14,15}, {0,1,4,5}, {2,3,6,7}, {0,1,4,5}, {2,3,6,7}, {8,9,12,13}, {10,11,14,15}, {8,9,12,13}, {10,11,14,15}};                      // 次状態
    const int output1[4][2] = {{0,1}, {0,1}, {1,0}, {1,0}};                         // 1ビット目の出力
    const int output2[4][2] = {{0,1}, {1,0}, {1,0}, {0,1}};                         // 2ビット目の出力
    const double infty = (double)n * 1000000000;                                    // 無限
    double branch_metric;                                                           // ブランチメトリック
    int likely_state;                                                               // トレースバックの最尤状態
    static node **nodes, *nodes_base;                                               // ノード
    static int *a_index_base;                                                       // 変調信号のメモリ
    static complex *autocor_base;                                                   // 相関関数のメモリ
    static int memory_flag;                                                         // メモリ管理フラグ
    int i, j1, j2, k1, k2, m;                                                       // ループカウンタ

    // メモリの確保
    if (memory_flag == 0) {
        nodes = (node **)malloc((n+1) * sizeof(node *));
        nodes_base = (node *)malloc((n+1) * num_state * sizeof(node));
        a_index_base = (int *)malloc((n+1) * num_state * num_state * n * sizeof(int));
        autocor_base = (complex *)malloc((n+1) * num_state * num_state * n * sizeof(complex));
        for (i = 0; i <= n; i++) {
            nodes[i] = nodes_base + i * num_state * num_state;
            for (j1 = 0; j1 < num_state; j1++) {
                for (j2 = 0; j2 < num_state; j2++) {
                    nodes[i][j1 * num_state + j2].a_index = a_index_base + (i * num_state * num_state + j1 * num_state + j2) * n;
                    nodes[i][j1 * num_state + j2].autocor = autocor_base + (i * num_state * num_state + j1 * num_state + j2) * n;
                }
            }
        }

        // フラグをチェック
        memory_flag = 1;
    }

    // 最初のノードメトリックを初期化
    for (j1 = 0; j1 < num_state; j1++) {
        for (j2 = 0; j2 < num_state; j1++) {
            nodes[0][num_state * j1 + j2].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // ノードメトリックを初期化
        for (j1 = 0; j1 < num_state; j1++) {
            for (j2 = 0; j2 < num_state; j2++) {
                nodes[i][num_state * j1 + j2].metric = infty;
            }
        }

        // ノードメトリックの計算
        for (j1 = 0; j1 < num_state; j1++) {
            for (j2 = 0; j2 < num_state; j2++) {
                for (k1 = 0; k1 < num_trans; k1++) {
                    for (k2 = 0; k2 < num_trans; k2++) {
                        // 仮の変調信号を求める
                        if (num_qam == 64) {
                            nodes[i-1][num_state * j1 + j2].a_index[i-1] = 32*c[6*(i-1)] + 16*c[6*(i-1)+1] + 8*(c[6*(i-1)+2] ^ output1[j1][k1]) + 4*(c[6*(i-1)+3] ^ output2[j1][k1]) + 2*(c[6*(i-1)+4] ^ output1[j2][k2]) + (c[4*(i-1)+5] ^ output2[j2][k2]);
                        } else if (num_qam == 256) {
                            nodes[i-1][num_state * j1 + j2].a_index[i-1] = 128*c[8*(i-1)] + 64*c[8*(i-1)+1] + 32*c[8*(i-1)+2] + 16*c[8*(i-1)+3] + 8*(c[8*(i-1)+4] ^ output1[j1][k1]) + 4*(c[8*(i-1)+5] ^ output2[j1][k1]) + 2*(c[8*(i-1)+6] ^ output1[j2][k2]) + (c[8*(i-1)+7] ^ output2[j2][k2]);
                        }

                        // ブランチメトリックを求める
                        if (i > 1) {
                            branch_metric = pow(cabs(a_caf[i-1] - constellation[nodes[i-1][num_state * j1 + j2].a_index[i-1]]), 2.0);
                        }

                        // パスの選択
                        if (nodes[i-1][num_state * j1 + j2].metric + branch_metric < nodes[i][next_state[num_state * j1 + j2][num_trans * k1 + k2]].metric) {
                            // メトリックの更新
                            nodes[i][next_state[num_state * j1 + j2][num_trans * k1 + k2]].metric = nodes[i-1][num_state * j1 + j2].metric + branch_metric;

                            // 前状態を保存
                            nodes[i][next_state[num_state * j1 + j2][num_trans * k1 + k2]].pre_state = num_state * j1 + j2;
                            nodes[i][next_state[num_state * j1 + j2][num_trans * k1 + k2]].a_index[i-1] = nodes[i-1][num_state * j1 + j2].a_index[i-1];
                        }
                    }
                }
            }
        }

        // 次状態に情報を渡す
        for (j1 = 0; j1 < num_state; j1++) {
            for (j2 = 0; j2 < num_state; j2++) {
                // 変調信号を更新する
                for (m = 0; m <= i-2; m++) {
                    nodes[i][num_state * j1 + j2].a_index[m] = nodes[i-1][nodes[i][num_qam * j1 + j2].pre_state].a_index[m];
                }
            }
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (j1 = 0; j1 < num_state; j1++) {
        for (j2 = 0; j2 < num_state; j2++) {
            if (nodes[n][num_qam * j1 + j2].metric < nodes[n][likely_state].metric) {
                likely_state = num_qam * j1 + j2;
            }
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
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
