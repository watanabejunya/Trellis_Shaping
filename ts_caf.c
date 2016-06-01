#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "lib_ts.h"
#include "env.h"

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif

#ifndef ENV_H
#define NUM_ARGUMENT 2                                              // 引数の数
#define NUM_OFDM 100000                                             // OFDMシンボルを送る回数
#define NUM_D 3                                                     // QAMのコンステレーション数
#define NUM_SUBCARRIER 1024                                         // サブキャリア数
#define OVER_SAMPLING_FACTOR 8                                      // オーバーサンプリング係数
#define CLIPPING_RATIO 1.6                                          // クリッピングの閾値
#define MAPPING_TYPE 1                                              // マッピングタイプ
#endif

#define NUM_C (NUM_D + 1)                                           // cのビット数
#define NUM_QAM ((int)pow(2.0, NUM_C))                              // QAMのコンステレーション数
#define NUM_S (NUM_C - 2)                                           // sのビット数
#define NUM_B 1                                                     // bのビット数
#define NUM_Z 2                                                     // zのビット数

int count_add = 0;
int count_mul = 0;

// マッピングを出力する
void run_mapping () {
    int *d;                                     // 情報ビット
    int *s, *z;                                 // 上位情報ビット
    int *b;                                     // 下位情報ビット
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
    complex *a_caf;                             // CAF後のOFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    s = (int *)malloc(NUM_S * NUM_SUBCARRIER * sizeof(int));
    z = (int *)malloc(NUM_Z * NUM_SUBCARRIER * sizeof(int));
    b = (int *)malloc(NUM_B * NUM_SUBCARRIER * sizeof(int));
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    a_caf = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 信号を生成
    make_signal(d, NUM_D * NUM_SUBCARRIER);

    // 信号を分離
    demultiplexer(d, s, b, NUM_D, NUM_S, NUM_B, NUM_SUBCARRIER);

    // 符号化
    inverse_parity_check_encoding(b, z, NUM_SUBCARRIER);

    // 信号を合成
    multiplexer(s, z, c, NUM_S, NUM_Z, NUM_C, NUM_SUBCARRIER);

    // 変調
    qam_modulation(c, a, NUM_SUBCARRIER, NUM_QAM, 1);

    // オーバーサンプリング
    over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

    // IFFT
    ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

    // マッピングを出力
    fp = fsopen("w", "./Result/raw_mapping_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);
    print_map(fp, t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

    // 変調
    qam_modulation_lsb(c, a, NUM_SUBCARRIER, NUM_QAM);

    // オーバーサンプリング
    over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

    // IFFT
    ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

    // クリッピング
    clipping(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

    // FFT
    fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

    // 減衰を補償する
    offset_attenuation(f, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

    // ダウンサンプリング
    down_sampling(f, a_caf, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

    // オーバーサンプリング
    over_sampling(a_caf, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

    // IFFT
    ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

    // マッピングを出力
    fp = fsopen("w", "./Result/caf_mapping_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);
    print_map(fp, t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

    // トレリスシェーピング
    trellis_shaping_caf(c, a_caf, a, NUM_SUBCARRIER, NUM_QAM);

    // オーバーサンプリング
    over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

    // IFFT
    ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

    // マッピングを出力
    fp = fsopen("w", "./Result/ts_mapping_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);
    print_map(fp, t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

    // メモリ解放
    free(d);
    free(s);
    free(z);
    free(b);
    free(a);
    free(a_caf);
    fftw_free(f);
    fftw_free(t);
}


// PAPR-CCDFを計算する
void run_calc_papr_ccdf () {
    int *d;                                     // 情報ビット
    int *s, *z;                                 // 上位情報ビット
    int *b;                                     // 下位情報ビット
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
    complex *a_caf;                             // CAF後のOFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double papr;                                // PAPR
    double ccdf;                                // CCDF
    double *pdf;                                // PAPRのPDF
    int min = 0, max = 10, num_index = 100;     // CCDFグラフの設定
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    s = (int *)malloc(NUM_S * NUM_SUBCARRIER * sizeof(int));
    b = (int *)malloc(NUM_B * NUM_SUBCARRIER * sizeof(int));
    z = (int *)malloc(NUM_Z * NUM_SUBCARRIER * sizeof(int));
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    a_caf = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    pdf = (double *)malloc(num_index * sizeof(double));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/prpr_ccdf_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);

    for (i = 0; i < NUM_OFDM; i++) {
        // 信号を生成
        make_signal(d, NUM_D * NUM_SUBCARRIER);

        // 信号を分離
        demultiplexer(d, s, b, NUM_D, NUM_S, NUM_B, NUM_SUBCARRIER);

        // 符号化
        inverse_parity_check_encoding(b, z, NUM_SUBCARRIER);

        // 信号を合成
        multiplexer(s, z, c, NUM_S, NUM_Z, NUM_C, NUM_SUBCARRIER);

        // 変調
        qam_modulation_lsb(c, a, NUM_SUBCARRIER, NUM_QAM);

        // オーバーサンプリング
        over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // IFFT
        ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

        // クリッピング
        clipping(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

        // FFT
        fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

        // 減衰を補償する
        offset_attenuation(f, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

        // ダウンサンプリング
        down_sampling(f, a_caf, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // トレリスシェーピング
        trellis_shaping_caf(c, a_caf, a, NUM_SUBCARRIER, NUM_QAM);

        // オーバーサンプリング
        over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // IFFT
        ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

        // PAPRの分布を入力
        count_papr_distribution(t, pdf, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, min, max, num_index);

        // 進捗を出力
        fprintf(stderr, "trial = %d   \r", i+1);
    }

    for (i = 0; i < num_index; i++) {
        // thresholdを設定
        papr = (double)min + (double)(i * (max - min)) / (double)num_index;

        // CCDFを求める
        ccdf = integrate_ccdf(pdf, papr, NUM_OFDM, min, max, num_index);

        // ファイル出力
        fprintf(fp, "%lf %e\n", papr, ccdf);
    }

    // 改行
    printf("\n");

    // メモリ解放
    free(d);
    free(s);
    free(b);
    free(c);
    free(z);
    free(a);
    free(a_caf);
    free(pdf);
    fftw_free(f);
    fftw_free(t);
}


// CCDFを計算する
void run_calc_normalized_ccdf () {
    int *d;                                     // 情報ビット
    int *s, *z;                                 // 上位情報ビット
    int *b;                                     // 下位情報ビット
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
    complex *a_caf;                             // CAF後のOFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double power;                               // 正規化瞬時電力
    double ccdf;                                // CCDF
    double *pdf;                                // PAPRのPDF
    int min = 0, max = 10, num_index = 100;     // CCDFグラフの設定
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    s = (int *)malloc(NUM_S * NUM_SUBCARRIER * sizeof(int));
    b = (int *)malloc(NUM_B * NUM_SUBCARRIER * sizeof(int));
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    z = (int *)malloc(NUM_Z * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    a_caf = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    pdf = (double *)malloc(num_index * sizeof(double));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/ccdf_normalized_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);

    for (i = 0; i < NUM_OFDM; i++) {
        // 信号を生成
        make_signal(d, NUM_D * NUM_SUBCARRIER);

        // 信号を分離
        demultiplexer(d, s, b, NUM_D, NUM_S, NUM_B, NUM_SUBCARRIER);

        // 符号化
        inverse_parity_check_encoding(b, z, NUM_SUBCARRIER);

        // 信号を合成
        multiplexer(s, z, c, NUM_S, NUM_Z, NUM_C, NUM_SUBCARRIER);

        // 変調
        qam_modulation_lsb(c, a, NUM_SUBCARRIER, NUM_QAM);

        // オーバーサンプリング
        over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // IFFT
        ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

        // クリッピング
        clipping(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

        // FFT
        fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

        // 減衰を補償する
        offset_attenuation(f, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

        // ダウンサンプリング
        down_sampling(f, a_caf, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // トレリスシェーピング
        trellis_shaping_caf(c, a_caf, a, NUM_SUBCARRIER, NUM_QAM);

        // オーバーサンプリング
        over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // IFFT
        ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

        // 正規化瞬時電力のpdfを求める
        count_normalized_distribution(t, pdf, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, min, max, num_index);

        // 進捗を出力
        fprintf(stderr, "trial = %d   \r", i+1);
    }

    for (i = 0; i < num_index; i++) {
        // thresholdを設定
        power = (double)min + (double)(i * (max - min)) / (double)num_index;

        // CCDFを求める
        ccdf = integrate_ccdf(pdf, power, NUM_OFDM, min, max, num_index);

        // ファイル出力
        fprintf(fp, "%lf %e\n", power, ccdf);
    }

    // 改行
    printf("\n");

    // メモリ解放
    free(d);
    free(s);
    free(b);
    free(c);
    free(z);
    free(a);
    free(a_caf);
    free(pdf);
    fftw_free(f);
    fftw_free(t);
}


// クリッピング比対PAPRのグラフを作成
void run_calc_clipping_ratio_characteristic () {
    int *d;                                     // 情報ビット
    int *s, *z;                                 // 上位情報ビット
    int *b;                                     // 下位情報ビット
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
    complex *a_caf;                             // CAF後のOFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double ratio;                               // クリッピング比
    double papr;                                // PAPR(dB)
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    s = (int *)malloc(NUM_S * NUM_SUBCARRIER * sizeof(int));
    z = (int *)malloc(NUM_Z * NUM_SUBCARRIER * sizeof(int));
    b = (int *)malloc(NUM_B * NUM_SUBCARRIER * sizeof(int));
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    a_caf = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/clipping_ratio_chara_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);

    for (ratio = 0.1; ratio < 2.0; ratio += 0.1) {

        // PAPRを初期化
        papr = 0.0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(d, NUM_D * NUM_SUBCARRIER);

            // 信号を分離
            demultiplexer(d, s, b, NUM_D, NUM_S, NUM_B, NUM_SUBCARRIER);

            // 符号化
            inverse_parity_check_encoding(b, z, NUM_SUBCARRIER);

            // 信号を合成
            multiplexer(s, z, c, NUM_S, NUM_Z, NUM_C, NUM_SUBCARRIER);

            // 変調
            qam_modulation_lsb(c, a, NUM_SUBCARRIER, NUM_QAM);

            // オーバーサンプリング
            over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // IFFT
            ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

            // クリッピング
            clipping(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, ratio);

            // FFT
            fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

            // 減衰を補償する
            offset_attenuation(f, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, ratio);

            // ダウンサンプリング
            down_sampling(f, a_caf, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // トレリスシェーピング
            trellis_shaping_caf(c, a_caf, a, NUM_SUBCARRIER, NUM_QAM);

            // オーバーサンプリング
            over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // IFFT
            ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

            // PAPRを計算
            papr += calc_papr_db(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

            // 進捗を出力
            fprintf(stderr, "ratio = %lf, trial = %d papr = %lf   \r", ratio, i+1, papr / (double)(i+1));
        }

        //PAPRの平均を求める
        papr /= (double)NUM_OFDM;

        // ファイル出力
        fprintf(fp, "%lf %lf\n", ratio, papr);

        // 改行
        printf("\n");
    }

    // メモリ解放
    free(d);
    free(s);
    free(z);
    free(b);
    free(a);
    free(c);
    free(a_caf);
    fftw_free(f);
    fftw_free(t);
}


// SLM候補点数対PAPRのグラフを作成
void run_calc_slm_characteristic () {
    int *d;                                     // 情報ビット
    int *s, *z;                                 // 上位情報ビット
    int *b;                                     // 下位情報ビット
    int *c;                                     // 符号語
    int *x, *y;                                 // インタリーバ
    complex *a;                                 // OFDMシンボル
    complex *a_caf;                             // CAF後のOFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    complex *t_opt;                             // SLMで選ばれたOFDMシンボル
    int num_slm;                                // SLM候補点数
    double papr_mean, papr_tmp, papr_opt;       // PAPR(dB)
    int i, j;                                   // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    s = (int *)malloc(NUM_S * NUM_SUBCARRIER * sizeof(int));
    b = (int *)malloc(NUM_B * NUM_SUBCARRIER * sizeof(int));
    z = (int *)malloc(NUM_Z * NUM_SUBCARRIER * sizeof(int));
    x = (int *)malloc(NUM_B * NUM_SUBCARRIER * sizeof(int));
    y = (int *)malloc(NUM_Z * NUM_SUBCARRIER * sizeof(int));
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    a_caf = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t_opt = (complex *)malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/slm_chara_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);

    for (num_slm = 1; num_slm <= 10; num_slm++) {

        // PAPRを初期化
        papr_mean = 0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(d, NUM_D * NUM_SUBCARRIER);

            // 信号を分離
            demultiplexer(d, s, b, NUM_D, NUM_S, NUM_B, NUM_SUBCARRIER);

            // 符号化
            inverse_parity_check_encoding(b, z, NUM_SUBCARRIER);

            // PAPRを初期化
            papr_opt = 100;

            for (j = 0; j < num_slm; j++) {
                // インタリーバを生成
                make_signal(x, NUM_B * NUM_SUBCARRIER);

                // インタリーバを符号化
                convolutional_encoding(x, y, NUM_SUBCARRIER);

                // 信号をインタリーブする
                xor_addition(z, y, NUM_Z * NUM_SUBCARRIER);

                // 信号を合成
                multiplexer(s, z, c, NUM_S, NUM_Z, NUM_C, NUM_SUBCARRIER);

                // 変調
                qam_modulation_lsb(c, a, NUM_SUBCARRIER, NUM_QAM);

                // オーバーサンプリング
                over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

                // IFFT
                ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

                // クリッピング
                clipping(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

                // FFT
                fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

                // 減衰を補償する
                offset_attenuation(f, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

                // ダウンサンプリング
                down_sampling(f, a_caf, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

                // トレリスシェーピング
                trellis_shaping_caf(c, a_caf, a, NUM_SUBCARRIER, NUM_QAM);

                // オーバーサンプリング
                over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

                // IFFT
                ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

                // PAPRを求める
                papr_tmp = calc_papr_db(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

                // 最もPAPRが低いものを保存
                if (papr_tmp < papr_opt) {
                    papr_opt = papr_tmp;
                    copy_complex(t, t_opt, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);
                }
            }

            papr_mean += papr_opt;

            // 進捗を出力
            fprintf(stderr, "num_slm = %d, trial = %d papr = %lf   \r", num_slm, i+1, papr_mean / (double)(i+1));
        }

        //PAPRの平均を求める
        papr_mean /= (double)NUM_OFDM;

        // ファイル出力
        fprintf(fp, "%d %lf\n", num_slm, papr_mean);

        // 改行
        printf("\n");
    }

    // メモリ解放
    free(d);
    free(s);
    free(z);
    free(b);
    free(x);
    free(y);
    free(c);
    free(a);
    free(a_caf);
    free(t_opt);
    fftw_free(f);
    fftw_free(t);
}


// BER特性のグラフを作成
void run_calc_ber () {
    int *d1, *d2;                               // 情報ビット
    int *s1, *s2;                               // 上位情報ビット
    int *z1, *z2;                               // 畳み込み符号
    int *b1, *b2;                               // 下位情報ビット
    int *c1, *c2;                               // 符号語
    complex *a;                                 // OFDMシンボル
    complex *a_caf;                             // CAF後のOFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double rate;                                // 符号化率
    int ebn0;                                   // E_b/N_0(DB)
    double snr;                                 // SNR(真値)
    double sigma;                               // ガウス雑音の分散
    double ber;                                 // BER
    double average_power;                       // 平均電力
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d1 = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    d2 = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    s1 = (int *)malloc(NUM_S * NUM_SUBCARRIER * sizeof(int));
    s2 = (int *)malloc(NUM_S * NUM_SUBCARRIER * sizeof(int));
    z1 = (int *)malloc(NUM_Z * NUM_SUBCARRIER * sizeof(int));
    z2 = (int *)malloc(NUM_Z * NUM_SUBCARRIER * sizeof(int));
    b1 = (int *)malloc(NUM_B * NUM_SUBCARRIER * sizeof(int));
    b2 = (int *)malloc(NUM_B * NUM_SUBCARRIER * sizeof(int));
    c1 = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    c2 = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    a_caf = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 符号化率を計算
    rate = (double)(NUM_D) / (double)(NUM_C);

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/ber_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);

    for (ebn0 = 6; ebn0 < 25; ebn0++) {
        // SNRを計算
        snr = pow(10, ((double)ebn0 / 10.0)) * rate * log2(NUM_QAM);

        // BERを初期化
        ber = 0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(d1, NUM_D * NUM_SUBCARRIER);

            // 信号を分離
            demultiplexer(d1, s1, b1, NUM_D, NUM_S, NUM_B, NUM_SUBCARRIER);

            // 符号化
            inverse_parity_check_encoding(b1, z1, NUM_SUBCARRIER);

            // 信号を合成
            multiplexer(s1, z1, c1, NUM_S, NUM_Z, NUM_C, NUM_SUBCARRIER);

            // 変調
            qam_modulation_lsb(c1, a, NUM_SUBCARRIER, NUM_QAM);

            // オーバーサンプリング
            over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // IFFT
            ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

            // クリッピング
            clipping(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

            // FFT
            fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

            // 減衰を補償する
            offset_attenuation(f, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

            // ダウンサンプリング
            down_sampling(f, a_caf, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // トレリスシェーピング
            trellis_shaping_caf(c1, a_caf, a, NUM_SUBCARRIER, NUM_QAM);

            // オーバーサンプリング
            over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // IFFT
            ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

            // 平均電力を求める
            average_power = calc_average_power(a, NUM_SUBCARRIER);

            // 雑音の分散を計算
            sigma = sqrt(average_power / (2.0*snr));

            // ガウス雑音を加える
            gaussian_noise(t, sigma, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

            // FFT
            fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

            // ダウンサンプリング
            down_sampling(f, a, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // 復調
            qam_demodulation(a, c2, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

            // 信号を分離
            demultiplexer(c2, s2, z2, NUM_C, NUM_S, NUM_Z, NUM_SUBCARRIER);

            // 復号
            parity_check_decoding(z2, b2, NUM_SUBCARRIER);

            // 信号を合成
            multiplexer(s2, b2, d2, NUM_S, NUM_B, NUM_D, NUM_SUBCARRIER);

            // BERを計算
            ber += count_be(d1, d2, NUM_D * NUM_SUBCARRIER) / (double)(NUM_D * NUM_SUBCARRIER);

            // 進捗を出力
            fprintf(stderr, "Eb/N0 = %.0d, trial = %d, BER = %e   \r", ebn0, i+1, ber / (double)(i+1));
        }

        // BERを計算
        ber = ber / (double)(NUM_OFDM);

        // ファイル出力
        fprintf(fp, "%d %e\n", ebn0, ber);

        // 改行
        printf("\n");
    }

    // メモリ解放
    free(d1);
    free(d2);
    free(s1);
    free(s2);
    free(z1);
    free(z2);
    free(b1);
    free(b2);
    free(c1);
    free(c2);
    free(a);
    free(a_caf);
    fftw_free(f);
    fftw_free(t);
}


// 時間を計算する
void run_calc_time () {
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
    complex *a_caf;                             // CAF後のOFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    clock_t start_time, end_time;               // 開始時間と終了時間
    double average_time;                        // 平均実行時間
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    a_caf = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/time_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);

    // 平均時間を初期化
    average_time = 0.0;

    for (i = 0; i < NUM_OFDM; i++) {
        // 信号を生成
        make_signal(c, NUM_C * NUM_SUBCARRIER);

        // 計測開始
        start_time = clock();

        // 変調
        qam_modulation_lsb(c, a, NUM_SUBCARRIER, NUM_QAM);

        // オーバーサンプリング
        over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // IFFT
        ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

        // クリッピング
        clipping(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

        // FFT
        fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

        // 減衰を補償する
        offset_attenuation(f, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

        // ダウンサンプリング
        down_sampling(f, a_caf, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // トレリスシェーピング
        trellis_shaping_caf(c, a_caf, a, NUM_SUBCARRIER, NUM_QAM);

        // 計測終了
        end_time = clock();

        // 実行時間を足す
        average_time += (double)(end_time - start_time) / CLOCKS_PER_SEC;

        // 進捗を出力
        fprintf(stderr, "trial = %d, average time = %e   \r", i+1, average_time / (double)(i+1));
    }

    // 平均時間を求める
    average_time /= (double)NUM_OFDM;

    // ファイル出力
    fprintf(fp, "%e\n", average_time);

    // 改行
    printf("\n");

    // メモリ解放
    free(c);
    free(a);
    free(a_caf);
    fftw_free(f);
    fftw_free(t);
}


// 計算量を計算する
void run_calc_cost () {
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
    complex *a_caf;                             // CAF後のOFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    a_caf = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/cost_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);

    // 信号を生成
    make_signal(c, NUM_C * NUM_SUBCARRIER);

    // 変調
    qam_modulation_lsb(c, a, NUM_SUBCARRIER, NUM_QAM);

    // オーバーサンプリング
    over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

    // IFFT
    ifftj(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

    // クリッピング
    clipping(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

    // FFT
    fftj(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

    // 減衰を補償する
    offset_attenuation(f, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, CLIPPING_RATIO);

    // ダウンサンプリング
    down_sampling(f, a_caf, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

    // トレリスシェーピング
    trellis_shaping_caf(c, a_caf, a, NUM_SUBCARRIER, NUM_QAM);

    // ファイル出力
    printf("add: %d\n", count_add);
    printf("mul: %d\n", count_mul);
    fprintf(fp, "add: %d\n", count_add);
    fprintf(fp, "mul: %d\n", count_mul);

    // メモリ解放
    free(c);
    free(a);
    free(a_caf);
    fftw_free(f);
    fftw_free(t);
}


int main (int argc,char *argv[]) {
    // 入力チェック
    if (argc != NUM_ARGUMENT) {
        printf("Wrong number of arguments (%d for %d).\n", argc, NUM_ARGUMENT);
        exit(-1);
    }

    // 入力チェックと実行
    if (strcmp(argv[1], "map") == 0) {
        printf("Make mapping graph.\n");
        run_mapping();
    } else if (strcmp(argv[1], "ccdf") == 0) {
        printf("Make PAPR - CCDF graph.\n");
        run_calc_papr_ccdf();
    } else if (strcmp(argv[1], "normal") == 0) {
        printf("Make normalized power - CCDF graph.\n");
        run_calc_normalized_ccdf();
    } else if (strcmp(argv[1], "ratio") == 0) {
        printf("Make clipping ratio - PAPR graph.\n");
        run_calc_clipping_ratio_characteristic();
    } else if (strcmp(argv[1], "slm") == 0) {
        printf("Make # of SLM - PAPR graph.\n");
        run_calc_slm_characteristic();
    } else if (strcmp(argv[1], "ber") == 0) {
        printf("Make Eb/N0 - BER graph.\n");
        run_calc_ber();
    } else if (strcmp(argv[1], "time") == 0) {
        printf("Calculate time\n");
        run_calc_time();
    } else if (strcmp(argv[1], "cost") == 0) {
        printf("Calculate cost\n");
        run_calc_cost();
    } else {
        printf("Invalid argument\n");
        exit(-1);
    }

    return 0;
}
