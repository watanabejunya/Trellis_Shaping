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


int main(void)
{
    int *c;                                     // 符号語
    complex *a;                                 // OFDMシンボル
    complex *a_caf;                             // CAF後のOFDMシンボル
    fftw_complex *f;                            // FFT用(周波数領域)
    fftw_complex *t;                            // FFT用(時間領域)
    FILE *fp;                                   // 出力用ファイルポインタ
    const int num_state = 4;                                                            // 状態数
    const int num_trans = 2;                                                            // 遷移数
    const int next_state[num_state][num_trans] = {{0,1}, {2,3}, {0,1}, {2,3}};          // 次状態
    const int output1[num_state][num_trans] = {{0,1}, {0,1}, {1,0}, {1,0}};             // 1ビット目の出力
    const int output2[num_state][num_trans] = {{0,1}, {1,0}, {1,0}, {0,1}};             // 2ビット目の出力
    const double infty = (double)NUM_SUBCARRIER * 10000000;                             // 無限
    int tmp_index;                                                                      // 合成した信号
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    static complex *constellation;                                                      // 変調信号のテーブル
    static complex **delta_table1;                                                      // δテーブル
    static node **nodes;                                                                // ノード
    static int memory_flag;                                                             // メモリ管理フラグ
    int sub, state, input;                                                              // ループカウンタ
    static double **metric_ts, **metric_tscaf;                                          // メトリック
    int i, j, l, k;                                                                     // ループカウンタ

    // メモリの確保
    c = (int *)malloc(NUM_C * NUM_SUBCARRIER * sizeof(int));
    a = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    a_caf = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
    f = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));
    t = (fftw_complex *)fftw_malloc(OVER_SAMPLING_FACTOR * NUM_SUBCARRIER * sizeof(fftw_complex));

    // 乱数の初期化
    srandom((unsigned)time(NULL));

    // 出力ファイルを開く
    fp = fsopen("w", "./Result/metrics_%d-QAM_%d-subs(TS_CAF_U1).dat", NUM_QAM, NUM_SUBCARRIER);

    for (k = 0; k < 1; k++) {
        // 信号を生成
        make_signal(c, NUM_C * NUM_SUBCARRIER);

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

        if (memory_flag == 0) {
            // メモリの確保
            constellation = (complex *)malloc(NUM_QAM * sizeof(complex));
            nodes = (node **)malloc((NUM_SUBCARRIER+1) * sizeof(node *));
            delta_table1 = (complex **)malloc(NUM_QAM * sizeof(complex *));
            metric_ts = (double **)malloc((NUM_SUBCARRIER+1) * num_state * sizeof(double *));
            metric_tscaf = (double **)malloc((NUM_SUBCARRIER+1) * num_state * sizeof(double *));

            for (i = 0; i <= NUM_SUBCARRIER; i++) {
                nodes[i] = (node *)malloc(num_state * sizeof(node));
                metric_ts[i] = (double *)malloc(num_state * sizeof(double));
                metric_tscaf[i] = (double *)malloc(num_state * sizeof(double));

                for (j = 0; j < num_state; j++) {
                    nodes[i][j].a_index = (int *)malloc(NUM_SUBCARRIER * sizeof(int));
                    nodes[i][j].autocor = (complex *)malloc(NUM_SUBCARRIER * sizeof(complex));
                }
            }
            for (i = 0; i < NUM_QAM; i++) {
                delta_table1[i] = (complex *)malloc(NUM_QAM * sizeof(complex));
            }

            // 信号点配置を初期化
            construct_constellation(constellation, NUM_QAM, 1);

            // フラグをチェック
            memory_flag = 1;
        }

        // 初期化
        for (i = 0; i <= NUM_SUBCARRIER; i++) {
            for (j = 0; j < num_state; j++) {
                nodes[i][j].metric = infty;
                metric_ts[i][j] = 0;
                metric_tscaf[i][j] = 0;
            }
        }
        nodes[0][0].metric = 0;

        // ビタビアルゴリズム開始
        for (i = 1; i <= NUM_SUBCARRIER; i++) {
            // 決定される信号位置
            sub = i - 1;

            // ノードメトリックの計算
            for (state = 0; state < num_state; state++) {
                if (nodes[i-1][state].metric > infty - 1.0) continue;

                for (input = 0; input < num_trans; input++) {
                    // 仮の変調信号を求める
                    if (NUM_QAM == 16) {
                        tmp_index = (c[4*sub] << 3) + (c[4*sub+1] << 2) + ((c[4*sub+2] ^ output1[state][input]) << 1) + (c[4*sub+3] ^ output2[state][input]);
                    } else if (NUM_QAM == 64) {
                        tmp_index = (c[6*sub] << 5) + (c[6*sub+1] << 4) + (c[6*sub+2] << 3) + (c[6*sub+3] << 2) + ((c[6*sub+4] ^ output1[state][input]) << 1) + (c[6*sub+5] ^ output2[state][input]);
                    } else if (NUM_QAM == 256) {
                        tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + ((c[8*sub+6] ^ output1[state][input]) << 1) + (c[8*sub+7] ^ output2[state][input]);
                    }

                    // ブランチメトリック
                    branch_metric = pow(cabs(a_caf[sub] - constellation[tmp_index]), 2.0);

                    // パスの選択
                    if (nodes[i-1][state].metric + branch_metric < nodes[i][next_state[state][input]].metric) {
                        // メトリックの更新
                        nodes[i][next_state[state][input]].metric = nodes[i-1][state].metric + branch_metric;

                        // 前状態を保存
                        nodes[i][next_state[state][input]].pre_state = state;
                        nodes[i][next_state[state][input]].a_index[sub] = tmp_index;
                    }
                }
            }

            // 次状態に情報を渡す
            for (state = 0; state < num_state; state++) {
                if (nodes[i][state].metric > infty - 1.0) continue;

                // 変調信号を更新する
                for (l = 0; l < sub; l++) {
                    nodes[i][state].a_index[l] = nodes[i-1][nodes[i][state].pre_state].a_index[l];
                }

                // 最尤の信号をaに入力
                for (l = 0; l < i; l++) {
                    a[l] = constellation[nodes[i][state].a_index[l]];
                }

                // TS CAFのメトリック
                metric_tscaf[i][state] = nodes[i][state].metric;

                // TSのブランチメトリック
                metric_ts[i][state] = sum_square_autocor(a, i);

                fprintf(fp, "%lf %lf\n", metric_tscaf[i][state] - metric_tscaf[i-1][nodes[i][state].pre_state], metric_ts[i][state] - metric_ts[i-1][nodes[i][state].pre_state]);
            }
        }

        // 最尤の最終状態を見つける
        likely_state = 0;
        for (state = 0; state < num_state; state++) {
            if (nodes[NUM_SUBCARRIER][state].metric < nodes[NUM_SUBCARRIER][likely_state].metric) {
                likely_state = state;
            }
        }

        // 最尤の信号をaに入力
        for (i = 0; i < NUM_SUBCARRIER; i++) {
            a[i] = constellation[nodes[NUM_SUBCARRIER][likely_state].a_index[i]];
        }
    }
}
