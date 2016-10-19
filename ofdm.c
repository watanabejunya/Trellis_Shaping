#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "env.h"
#include "lib_base.h"
#include "lib_ts.h"
#include "lib_caf.h"
#include "lib_peak.h"
#include "lib_amp.h"

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
#define MAPPING_TYPE 1                                              // トレリスシェーピングのマッピングタイプ
#define AMP_TYPE AMP_TYPE_IDEAL                                     // アンプの種類
#define AMP_CLASS AMP_CLASS_A                                       // アンプの級数
#define IS_EFFECTIVE 0                                              // effectiveフラグ
#endif

#define NUM_C (NUM_D + 1)                                           // cのビット数
#define NUM_QAM ((int)pow(2.0, NUM_C))                              // QAMのコンステレーション数


// PAPR-CCDFを計算する
void run_calc_papr_ccdf () {
    int *c;                                     // 情報ビット
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double papr;                                // PAPR
    double ccdf;                                // CCDF
    double *pdf;                                // PAPRのPDF
    int min = 0, max = 10, num_index = 100;     // CCDFグラフの設定
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    pdf = (double *)malloc(num_index * sizeof(double));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/prpr_ccdf_%d-QAM_%d-subs(OFDM).dat", NUM_QAM, NUM_SUBCARRIER);

    for (i = 0; i < NUM_OFDM; i++) {
        // 信号を生成
        make_signal(c, NUM_C * NUM_SUBCARRIER);

        // 変調
        qam_modulation(c, a, NUM_SUBCARRIER, NUM_QAM);

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
    free(c);
    free(a);
    free(pdf);
    fftw_free(f);
    fftw_free(t);
}


// 正規化瞬時電力-CCDFを計算する
void run_calc_normalized_ccdf () {
    int *c;                                     // 情報ビット
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double power;                               // 正規化瞬時電力
    double ccdf;                                // CCDF
    double *pdf;                                // PAPRのPDF
    int min = 0, max = 10, num_index = 100;     // CCDFグラフの設定
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    pdf = (double *)malloc(num_index * sizeof(double));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/ccdf_normalized_%d-QAM_%d-subs(OFDM).dat", NUM_QAM, NUM_SUBCARRIER);

    for (i = 0; i < NUM_OFDM; i++) {
        // 信号を生成
        make_signal(c, NUM_C * NUM_SUBCARRIER);

        // 変調
        qam_modulation(c, a, NUM_SUBCARRIER, NUM_QAM);

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
    free(c);
    free(a);
    free(pdf);
    fftw_free(f);
    fftw_free(t);
}


// 平均PAPRを計算
void run_calc_mean_papr () {
    int *c;                                     // 情報
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double papr;                                // PAPR
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/mean_papr_%d-QAM_%d-subs(OFDM).dat", NUM_QAM, NUM_SUBCARRIER);

    // PAPRを初期化
    papr = 0.0;

    for (i = 0; i < NUM_OFDM; i++) {
        // 信号を生成
        make_signal(c, NUM_C * NUM_SUBCARRIER);

        // 変調
        qam_modulation(c, a, NUM_SUBCARRIER, NUM_QAM);

        // オーバーサンプリング
        over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // IFFT
        ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

        // PAPRを求める
        papr += calc_papr_db(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

        // 進捗を出力
        fprintf(stderr, "trial = %d, PAPR = %lf   \r", i+1, papr / (double)(i+1));
    }

    // PAPRを計算
    papr /= (double)(NUM_OFDM);

    // ファイル出力
    fprintf(fp, "%lf\n", papr);

    // 改行
    printf("\n");

    // メモリ解放
    free(c);
    free(a);
    fftw_free(f);
    fftw_free(t);
}


// 平均PAPRを計算
void run_calc_average_power () {
    int *c;                                     // 情報
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double power;                               // 平均電力
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/average_power_%d-QAM_%d-subs(OFDM).dat", NUM_QAM, NUM_SUBCARRIER);

    // PAPRを初期化
    power = 0.0;

    for (i = 0; i < NUM_OFDM; i++) {
        // 信号を生成
        make_signal(c, NUM_C * NUM_SUBCARRIER);

        // 変調
        qam_modulation(c, a, NUM_SUBCARRIER, NUM_QAM);

        // オーバーサンプリング
        over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // IFFT
        ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

        // PAPRを求める
        power += calc_average_power(t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

        // 進捗を出力
        fprintf(stderr, "trial = %d, power = %lf   \r", i+1, power / (double)(i+1));
    }

    // PAPRを計算
    power /= (double)(NUM_OFDM);

    // ファイル出力
    fprintf(fp, "%lf\n", power);

    // 改行
    printf("\n");

    // メモリ解放
    free(c);
    free(a);
    fftw_free(f);
    fftw_free(t);
}


// BER特性のグラフを作成
void run_calc_ber () {
    int *c1, *c2;                               // 情報
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double rate;                                // 符号化率
    int ebn0;                                   // E_b/N_0(DB)
    double snr;                                 // SNR(真値)
    double sigma;                               // ガウス雑音の分散
    double ber;                                 // BER
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    c1 = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    c2 = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 符号化率を計算
    rate = 1.0;

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/ber_%d-QAM_%d-subs(OFDM).dat", NUM_QAM, NUM_SUBCARRIER);

    for (ebn0 = 6; ebn0 < 25; ebn0++) {
        // SNRを計算
        snr = pow(10, ((double)ebn0 / 10.0)) * rate * log2(NUM_QAM);
        sigma = sqrt(((double)(NUM_QAM - 1) * 2.0 / 3.0) / (double)(2.0*snr));

        // BERを初期化
        ber = 0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(c1, NUM_C * NUM_SUBCARRIER);

            // 変調
            qam_modulation(c1, a, NUM_SUBCARRIER, NUM_QAM);

            // オーバーサンプリング
            over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // IFFT
            ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

            // ガウス雑音を加える
            gaussian_noise(t, sigma, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

            // FFT
            fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t, f);

            // ダウンサンプリング
            down_sampling(f, a, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // 復調
            qam_demodulation(a, c2, NUM_SUBCARRIER, NUM_QAM);

            // BERを計算
            ber += count_bit_error(c1, c2, NUM_C * NUM_SUBCARRIER) / (double)(NUM_C * NUM_SUBCARRIER);

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
    free(c1);
    free(c2);
    free(a);
    fftw_free(f);
    fftw_free(t);
}


// 増幅器効率を計算する
void run_calc_pa_efficiency()
{
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    complex *t_amp;                             // 増幅後の信号
    double ibo;                                 // IBO
    double sdr;                                 // SDR
    double efficiency;                          // 増幅器効率
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t_amp = (complex *)malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/pa_efficiency_%d-QAM_%d-subs_amp_type-%d_class-%d(OFDM).dat", NUM_QAM, NUM_SUBCARRIER, AMP_TYPE, AMP_CLASS);

    for (ibo = -15.0; ibo <= 18.0; ibo++) {
        // 初期化
        efficiency = 0;
        sdr = 0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(c, NUM_C * NUM_SUBCARRIER);

            // 変調
            qam_modulation(c, a, NUM_SUBCARRIER, NUM_QAM);

            // オーバーサンプリング
            over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

            // IFFT
            ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

            // 増幅器に通す
            amplify_signal(t, t_amp, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, ibo, AMP_TYPE);

            // SDRを計算
            sdr += calc_sdr_design(t, t_amp, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

            // 増幅器効率を計算
            efficiency += calc_pa_efficiency(t, t_amp, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, ibo, AMP_CLASS);

            // 進捗を出力
            fprintf(stderr, "ibo = %.1lf [dB] trial = %d sdr = %lf [dB] efficiency = %lf [%%]   \r", ibo, i+1, sdr / (double)(i+1), efficiency / (double)(i+1) * 100.0);
        }

        // ファイル出力
        fprintf(fp, "%lf %lf\n", sdr / (double)NUM_OFDM, efficiency / (double)NUM_OFDM * 100.0);

        // 改行
        printf("\n");
    }

    // メモリ解放
    free(c);
    free(a);
    fftw_free(f);
    fftw_free(t);
    free(t_amp);
}


void run_calc_out_of_band_radiation()
{
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    complex *t_amp;                             // 増幅後の信号
    double *psd;                                // 電力密度スペクトル
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t_amp = (complex *)malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(complex));
    psd = (double *)malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(double));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/out_of_band_%d-QAM_%d-subs_amp_type-%d(OFDM).dat", NUM_QAM, NUM_SUBCARRIER, AMP_TYPE);

    for (i = 0; i < NUM_OFDM; i++) {
        // 信号を生成
        make_signal(c, NUM_C * NUM_SUBCARRIER);

        // 変調
        qam_modulation(c, a, NUM_SUBCARRIER, NUM_QAM);

        // オーバーサンプリング
        over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

        // IFFT
        ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

        // 増幅器に通す
        amplify_signal(t, t_amp, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, IBO, AMP_TYPE);

        // FFT
        fft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, t_amp, f);

        // PSDを算出する
        calc_power_spectral_density(f, psd, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, NUM_OFDM);

        // 進捗を出力
        fprintf(stderr, "trial = %d  \r", i+1);
    }

    // 出力
    print_double(fp, psd, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

    // 改行
    printf("\n");

    // メモリ解放
    free(c);
    free(a);
    fftw_free(f);
    fftw_free(t);
    free(t_amp);
    free(psd);
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
    } else if (strcmp(argv[1], "papr") == 0) {
        printf("Calculate mean PAPR.\n");
        run_calc_mean_papr();
    } else if (strcmp(argv[1], "power") == 0) {
        printf("Calculate average power.\n");
        run_calc_average_power();
    } else if (strcmp(argv[1], "ber") == 0) {
        printf("Make Eb/N0 - BER graph.\n");
        run_calc_ber();
    } else if (strcmp(argv[1], "eff") == 0) {
        printf("Make IBO - efficiency graph.\n");
        run_calc_pa_efficiency();
    } else if (strcmp(argv[1], "out") == 0) {
        printf("Make PSD graph.\n");
        run_calc_out_of_band_radiation();
    } else {
        printf("Invalid argument\n");
        exit(-1);
    }

    return 0;
}
