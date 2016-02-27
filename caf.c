#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "lib_ts.h"

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif

#define PI M_PI                                                     // PI
#define NUM_ARGUMENT 2                                              // 引数の数
#define NUM_QAM 64                                                  // QAMのコンステレーション数
#define NUM_D ((int)log2(NUM_QAM))                                  // dのビット数
#define NUM_SUBCARRIER 64                                           // サブキャリア数
#define NUM_OFDM 100000                                              // OFDMシンボルを送る回数
#define OVER_SAMPLING_FACTOR 8                                      // オーバーサンプリング係数
#define CLIPPING_RATIO 1.2                                          // クリッピング比
#define MAPPING_TYPE 1                                              // トレリスシェーピングのマッピングタイプ


// PAPR-CCDFを計算する
void run_calc_papr_ccdf () {
    int *d;                                     // 情報ビット
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double papr, papr_tmp;                      // PAPR
    double ccdf;                                // CCDF
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/prpr_ccdf_%d-QAM_%d-subs(CAF).dat", NUM_QAM, NUM_SUBCARRIER);

    for (papr = 3.0; papr < 13.0; papr += 0.5) {

        // CCDFを初期化
        ccdf = 0.0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(d, NUM_D * NUM_SUBCARRIER);

            // 変調
            qam_modulation(d, a, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

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
            down_sampling(f, a, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // オーバーサンプリング
            over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // IFFT
            ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

            // PAPRを求める
            papr_tmp = calc_papr_db(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

            // CCDFを計算
            if (papr < papr_tmp) {
                ccdf += 1.0;
            }

            // 進捗を出力
            fprintf(stderr, "papr = %lf, trial = %d ccdf = %e   \r", papr, i+1, ccdf / (double)(i+1));
        }

        // CCDFの平均を求める
        ccdf /= (double)NUM_OFDM;

        // ファイル出力
        fprintf(fp, "%lf %e\n", papr, ccdf);

        // 改行
        printf("\n");
    }

    // メモリ解放
    free(d);
    free(a);
    fftw_free(f);
    fftw_free(t);
}


// 正規化瞬時電力-CCDFを計算する
void run_calc_normalized_ccdf () {
    int *d;                                     // 情報ビット
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double power;                               // 正規化瞬時電力
    double ccdf;                                // CCDF
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/ccdf_normalized_%d-QAM_%d-subs(CAF).dat", NUM_QAM, NUM_SUBCARRIER);

    for (power = 2.0; power < 9.5; power += 0.5) {

        // CCDFを初期化
        ccdf = 0.0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(d, NUM_D * NUM_SUBCARRIER);

            // 変調
            qam_modulation(d, a, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

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
            down_sampling(f, a, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // オーバーサンプリング
            over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // IFFT
            ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

            // CCDFを計算
            ccdf += calc_normalized_ccdf(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, power);

            // 進捗を出力
            fprintf(stderr, "power = %lf, trial = %d ccdf = %e   \r", power, i, ccdf / (double)(i+1));
        }

        // CCDFの平均を求める
        ccdf /= (double)NUM_OFDM;

        // ファイル出力
        fprintf(fp, "%lf %e\n", power, ccdf);

        // 改行
        printf("\n");
    }

    // メモリ解放
    free(d);
    free(a);
    fftw_free(f);
    fftw_free(t);
}


// 平均PAPRを計算
void run_calc_clipping_ratio_characteristic () {
    int *d;                                     // 情報
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double ratio;                               // クリッピング比
    double papr;                                // PAPR
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/clipping_ratio_chara_%d-QAM_%d-subs(CAF).dat", NUM_QAM, NUM_SUBCARRIER);

    for (ratio = 0.1; ratio < 2.0; ratio += 0.1) {

        // PAPRを初期化
        papr = 0.0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(d, NUM_D * NUM_SUBCARRIER);

            // 変調
            qam_modulation(d, a, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

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
            down_sampling(f, a, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // オーバーサンプリング
            over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // IFFT
            ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

            // PAPRを求める
            papr += calc_papr_db(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

            // 進捗を出力
            fprintf(stderr, "ratio = %lf, trial = %d papr = %lf   \r", ratio, i, papr / (double)(i+1));
        }

        // PAPRを計算
        papr /= (double)(NUM_OFDM);

        // ファイル出力
        fprintf(fp, "%lf %lf\n", ratio, papr);

        // 改行
        printf("\n");
    }

    // メモリ解放
    free(d);
    free(a);
    fftw_free(f);
    fftw_free(t);
}


// BER特性のグラフを作成
void run_calc_ber () {
    int *d1, *d2;                               // 情報
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double rate;                                // 符号化率
    int ebn0;                                   // E_b/N_0(DB)
    double snr;                                 // SNR(真値)
    double sigma;                               // ガウス雑音の分散
    double average_power;                       // 平均電力
    double ber;                                 // BER
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d1 = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    d2 = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 符号化率を計算
    rate = 1.0;

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/ber_%d-QAM_%d-subs(CAF).dat", NUM_QAM, NUM_SUBCARRIER);

    for (ebn0 = 6; ebn0 < 25; ebn0++) {
        // SNRを計算
        snr = pow(10, ((double)ebn0 / 10.0)) * rate * log2(NUM_QAM);

        // BERを初期化
        ber = 0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(d1, NUM_D * NUM_SUBCARRIER);

            // 変調
            qam_modulation(d1, a, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

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
            down_sampling(f, a, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

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
            qam_demodulation(a, d2, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

            // BERを計算
            ber += count_be(d1, d2, NUM_D * NUM_SUBCARRIER) / (double)(NUM_D * NUM_SUBCARRIER);

            // 進捗を出力
            fprintf(stderr, "Eb/N0 = %d, trial = %d, BER = %e   \r", ebn0, i+1, ber / (double)(i+1));
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
    free(a);
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
    if (strcmp(argv[1], "ccdf") == 0) {
        printf("Make PAPR - CCDF graph.\n");
        run_calc_papr_ccdf();
    } else if (strcmp(argv[1], "normal") == 0) {
        printf("Make normalized power - CCDF graph.\n");
        run_calc_normalized_ccdf();
    } else if (strcmp(argv[1], "ratio") == 0) {
        printf("Make clipping ratio - PAPR graph.\n");
        run_calc_clipping_ratio_characteristic();
    } else if (strcmp(argv[1], "ber") == 0) {
        printf("Make Eb/N0 - BER graph.\n");
        run_calc_ber();
    } else {
        printf("Invalid argument\n");
        exit(-1);
    }

    return 0;
}
