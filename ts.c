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

#define NUM_ARGUMENT 2                                              // 引数の数
#define NUM_OFDM 30000                                              // OFDMシンボルを送る回数
#define NUM_QAM 16                                                 // QAMのコンステレーション数
#define NUM_D 3                                                     // dのビット数
#define NUM_S 1                                                     // sのビット数
#define NUM_Z 2                                                     // zのビット数
#define NUM_B 2                                                     // bのビット数
#define NUM_C (NUM_Z + NUM_B)                                       // cのビット数
#define NUM_SUBCARRIER 64                                           // サブキャリア数
#define OVER_SAMPLING_FACTOR 8                                      // オーバーサンプリング係数
#define MAPPING_TYPE 1                                              // トレリスシェーピングのマッピングタイプ


// マッピングを出力する
void run_mapping () {
    int *d;                                     // 情報ビット
    int *s, *z;                                 // 上位情報ビット
    int *b;                                     // 下位情報ビット
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
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
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // コンステレーションを初期化
    construct_constellation(NUM_QAM, MAPPING_TYPE);

    // 信号を生成
    make_signal(d, NUM_D * NUM_SUBCARRIER);

    // 信号を分離
    demultiplexer(d, s, b, NUM_D, NUM_S, NUM_B, NUM_SUBCARRIER);

    // 符号化
    inverse_parity_check_coding(s, z, NUM_S * NUM_SUBCARRIER);

    // 信号を合成
    multiplexer(z, b, c, NUM_Z, NUM_B, NUM_C, NUM_SUBCARRIER);

    // 変調
    qam_modulation(c, a, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

    // オーバーサンプリング
    over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

    // IFFT
    ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

    // マッピングを出力
    fp = fsopen("w", "raw_mapping_%d-QAM_%d-subcarriers(extend).dat", NUM_QAM, NUM_SUBCARRIER);
    print_map(fp, t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

    // トレリスシェーピング
    trellis_shaping(c, a, NUM_SUBCARRIER, NUM_QAM);

    // オーバーサンプリング
    over_sampling(a, f, OVER_SAMPLING_FACTOR, NUM_SUBCARRIER);

    // IFFT
    ifft(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER, f, t);

    // マッピングを出力
    fp = fsopen("w", "ts_mapping_%d-QAM_%d-subcarriers(extend).dat", NUM_QAM, NUM_SUBCARRIER);
    print_map(fp, t, OVER_SAMPLING_FACTOR * NUM_SUBCARRIER);

    // メモリ解放
    free(constellation);
    free(d);
    free(s);
    free(z);
    free(b);
    free(c);
    free(a);
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
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    double power;                               // 正規化瞬時電力
    double ccdf;                                // CCDF
    int i;                                      // ループカウンタ
    FILE *fp;                                   // 出力用ファイルポインタ

    // メモリの確保
    d = (int *)malloc(NUM_D * NUM_SUBCARRIER * sizeof(int));
    s = (int *)malloc(NUM_S * NUM_SUBCARRIER * sizeof(int));
    b = (int *)malloc(NUM_B * NUM_SUBCARRIER * sizeof(int));
    z = (int *)malloc(NUM_Z * NUM_SUBCARRIER * sizeof(int));
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // コンステレーションを初期化
    construct_constellation(NUM_QAM, MAPPING_TYPE);

    // 出力ファイルを開く
    fp = fsopen("w", "ccdf_normalized_%d-QAM_%d-subcarriers(extend).dat", NUM_QAM, NUM_SUBCARRIER);

    for (power = 2.0; power < 9.5; power += 0.5) {

        // CCDFを初期化
        ccdf = 0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(d, NUM_D * NUM_SUBCARRIER);

            // 信号を分離
            demultiplexer(d, s, b, NUM_D, NUM_S, NUM_B, NUM_SUBCARRIER);

            // 符号化
            inverse_parity_check_coding(s, z, NUM_S * NUM_SUBCARRIER);

            // 信号を合成
            multiplexer(z, b, c, NUM_Z, NUM_B, NUM_C, NUM_SUBCARRIER);

            // 変調
            qam_modulation(c, a, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

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
    free(constellation);
    free(d);
    free(s);
    free(b);
    free(z);
    free(c);
    free(a);
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
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 符号化率を計算
    rate = (double)(NUM_D) / (double)(NUM_Z + NUM_B);

    // コンステレーションを初期化
    construct_constellation(NUM_QAM, MAPPING_TYPE);

    // 出力ファイルを開く
    fp = fsopen("w", "ber_%d-QAM_%d-subcarriers(extend).dat", NUM_QAM, NUM_SUBCARRIER);

    for (ebn0 = 6; ebn0 < 26; ebn0++) {
        // E_b/N_0を計算
        snr = pow(10, ((double)ebn0 / 10.0)) * rate * log2(NUM_QAM);
        sigma = sqrt(((double)(NUM_QAM - 1) * 2.0 / 3.0) / (double)(2.0*snr));

        // BERを初期化
        ber = 0;

        for (i = 0; i < NUM_OFDM; i++) {
            // 信号を生成
            make_signal(d1, NUM_D * NUM_SUBCARRIER);

            // 信号を分離
            demultiplexer(d1, s1, b1, NUM_D, NUM_S, NUM_B, NUM_SUBCARRIER);

            // 符号化
            inverse_parity_check_coding(s1, z1, NUM_S * NUM_SUBCARRIER);

            // 信号を合成
            multiplexer(z1, b1, c1, NUM_Z, NUM_B, NUM_C, NUM_SUBCARRIER);

            // 変調
            qam_modulation(c1, a, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

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
            qam_demodulation(a, c2, NUM_SUBCARRIER, NUM_QAM, MAPPING_TYPE);

            // 信号を分離
            demultiplexer(c2, z2, b2, NUM_C, NUM_Z, NUM_B, NUM_SUBCARRIER);

            // 復号
            parity_check_decoding(z2, s2, NUM_S * NUM_SUBCARRIER);

            // 信号を合成
            multiplexer(s2, b2, d2, NUM_S, NUM_B, NUM_D, NUM_SUBCARRIER);

            // BERを計算
            ber += count_be(d1, d2, NUM_D * NUM_SUBCARRIER) / (double)(NUM_D * NUM_SUBCARRIER);

            // 進捗を出力
            fprintf(stderr, "Eb/N0 = %d, trial = %d, BER = %e   \r", ebn0, i, ber / (double)(i+1));
        }

        // BERを計算
        ber = ber / (double)(NUM_OFDM);

        // ファイル出力
        fprintf(fp, "%d %e\n", ebn0, ber);

        // 改行
        printf("\n");
    }

    // メモリ解放
    free(constellation);
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
        // マッピングを出力
        printf("Make mapping graph.\n");
        run_mapping();
    } else if (strcmp(argv[1], "ccdf") == 0) {
        printf("Make normalized power - CCDF graph.\n");
        run_calc_normalized_ccdf();
    } else if (strcmp(argv[1], "ber") == 0) {
        printf("Make Eb/N0 - BER graph.\n");
        run_calc_ber();
    } else {
        printf("Invalid argument\n");
        exit(-1);
    }

    return 0;
}
