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
#define NUM_OFDM 30000                                              // OFDMシンボルを送る回数
#define NUM_QAM 64                                                 // QAMのコンステレーション数
#define NUM_D 6                                                     // dのビット数
#define NUM_SUBCARRIER 64                                           // サブキャリア数
#define OVER_SAMPLING_FACTOR 8                                      // オーバーサンプリング係数
#define MAPPING_TYPE 1                                           // トレリスシェーピングのマッピングタイプ


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

    // コンステレーションを初期化
    construct_constellation(NUM_QAM, MAPPING_TYPE);

    // 出力ファイルを開く
    fp = fsopen("w", "ber_%d-QAM_%d-subcarriers(extend).dat", NUM_QAM, NUM_SUBCARRIER);

    for (ebn0 = 6; ebn0 < 21; ebn0++) {
        // E_b/N_0を計算
        snr = pow(10, ((double)ebn0 / 10.0)) * rate * log2(NUM_QAM);
        sigma = sqrt(((double)(NUM_QAM - 1) * 2.0 / 3.0) / (double)(2.0*snr));

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
    if (strcmp(argv[1], "ber") == 0) {
        printf("Make Eb/N0 - BER graph.\n");
        run_calc_ber();
    } else {
        printf("Invalid argument\n");
        exit(-1);
    }

    return 0;
}
