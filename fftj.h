#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif

extern int count_add;
extern int count_mul;

double dadd(double a, double b);
double dmul(double a, double b);
complex cadd(complex a, complex b);
complex cmul(complex a, complex b);
void fftj(int n, complex *in, complex *out);
void ifftj(int n, complex *in, complex *out);

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

complex cadd(complex a, complex b)
{
    return dadd(creal(a), creal(b)) + I * dadd(cimag(a), cimag(b));
}

complex cmul(complex a, complex b)
{
    return dadd(dmul(creal(a), creal(b)), - dmul(cimag(a), cimag(b))) + I * dadd(dmul(creal(a), cimag(b)), dmul(cimag(a), creal(b)));
}

void fftj(int n, complex *in, complex *out)
{
    int exponent = (int)log2(n);                        // nの冪指数
    static int memory_flag;                             // メモリ管理フラグ
    static int *index;                                  // ビットリバースした配列のインデックス
    static complex *w;                                  // 回転演算子
    static complex **data;                              // 計算用データ配列
    int start, end;                                     // ビットリバースループで使う作業用変数
    int j_roop, k_roop, half;                           // バタフライ演算で使う作業用変数
    int tmp_index1, tmp_index2;                         // 一時的な配列のインデックス
    int i, j, k;                                        // 作業用変数

    if (memory_flag == 0) {
        // 入力チェック
        if (n < 2 || n != (int)pow(2.0, exponent)) {
            fprintf(stderr, "n should be the k-th power of 2.\n");
            exit(-1);
        }

        // メモリ確保
        index = (int *) malloc(n * sizeof(int));
        w = (complex *) malloc(n * sizeof(complex));
        data = (complex **) malloc(n * exponent * sizeof(complex *));
        for (i = 0; i <= exponent; i++) {
            data[i] = (complex *)malloc(n * sizeof(complex));
        }

        // ビットリバースしたときのインデックスを作っておく
        index[0] = 0;
        for (i = 1; i <= exponent; i++) {
            // ループの範囲を定める
            start = (int)pow(2.0, i - 1);
            end = (int)pow(2.0, i);

            for (j = start; j < end; j++) {
                // ビットリバースした値を求める
                index[j] = index[j - start] + n / end;
            }
        }

        // 回転因子の計算
        for (i = 0; i < n; i++) {
            w[i] = cexp(- I * 2.0 * M_PI * i / n);
        }
    }

    // 入力をバタフライ演算用に入れ替える
    for (i = 0; i < n; i++) {
        data[0][i] = in[index[i]];
    }

    // バタフライ演算
    for (i = 1; i <= exponent; i++) {
        // ループの範囲を設定
        j_roop = n / (int)pow(2.0, i);
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
                data[i][j * k_roop + k] = data[i-1][tmp_index1] + data[i-1][tmp_index2] * w[k * j_roop];
            }
        }
    }

    // データを正規化
    for (i = 0; i < n; i++) {
        out[i] = data[exponent][i];// / sqrt(n);
    }
}

void ifftj(int n, complex *in, complex *out)
{
    int exponent = (int)log2(n);                        // nの冪指数
    static int memory_flag;                             // メモリ管理フラグ
    static int *index;                                  // ビットリバースした配列のインデックス
    static complex *w;                                  // 回転演算子
    static complex **data;                              // 計算用データ配列
    int start, end;                                     // ビットリバースループで使う作業用変数
    int j_roop, k_roop, half;                           // バタフライ演算で使う作業用変数
    int tmp_index1, tmp_index2;                         // 一時的な配列のインデックス
    int i, j, k;                                        // 作業用変数

    if (memory_flag == 0) {
        // 入力チェック
        if (n < 2 || n != (int)pow(2.0, exponent)) {
            fprintf(stderr, "n should be the k-th power of 2.\n");
            exit(-1);
        }

        // メモリ確保
        index = (int *) malloc(n * sizeof(int));
        w = (complex *) malloc(n * sizeof(complex));
        data = (complex **) malloc(n * exponent * sizeof(complex *));
        for (i = 0; i <= exponent; i++) {
            data[i] = (complex *)malloc(n * sizeof(complex));
        }

        // ビットリバースしたときのインデックスを作っておく
        index[0] = 0;
        for (i = 1; i <= exponent; i++) {
            // ループの範囲を定める
            start = (int)pow(2.0, i - 1);
            end = (int)pow(2.0, i);

            for (j = start; j < end; j++) {
                // ビットリバースした値を求める
                index[j] = index[j - start] + n / end;
            }
        }

        // 回転因子の計算
        for (i = 0; i < n; i++) {
            w[i] = cexp(I * 2.0 * M_PI * i / n);
        }
    }

    // 入力をバタフライ演算用に入れ替える
    for (i = 0; i < n; i++) {
        data[0][i] = in[index[i]];
    }

    // バタフライ演算
    for (i = 1; i <= exponent; i++) {
        // ループの範囲を設定
        j_roop = n / (int)pow(2.0, i);
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
                data[i][j * k_roop + k] = data[i-1][tmp_index1] + data[i-1][tmp_index2] * w[k * j_roop];
            }
        }
    }

    // データを正規化
    for (i = 0; i < n; i++) {
        out[i] = data[exponent][i];// / sqrt(n);
    }
}
