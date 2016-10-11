#ifndef LIB_TS_H
#define LIB_TS_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "lib_base.h"

#ifndef TYPE_COMPLEX
#define TYPE_COMPLEX
#undef complex
typedef _Complex double complex;
#endif

#ifndef WINDOW_SIZE
#define WINDOW_SIZE 64
#endif

// TSのノード
typedef struct {
    int *a_index;
    complex *autocor;
    double metric;
    int pre_state;
} node;


void demultiplexer(int *d, int *s, int *b, int num_d, int num_s, int num_b, int n);
void multiplexer(int *s, int *b, int *d, int num_s, int num_b, int num_d, int n);
void convolutional_encoding2(int *u, int *c, int n);
void convolutional_encoding3(int *u, int *c, int n);
void parity_check_decoding2(int *c, int *u, int n);
void parity_check_decoding3(int *c, int *u, int n);
void inverse_parity_check_encoding2(int *u, int *c, int n);
void inverse_parity_check_encoding3(int *u, int *c, int n);
void construct_constellation(complex *constellation, int m, int type);
void cbts_qam_modulation(int *c, complex *a, int n, int m_m, int m_l, int type);
void ts_qam_modulation(int *c, complex *a, int n, int m, int type);
void ts_qam_demodulation(complex *a, int *c, int n, int m, int type);
void trellis_shaping(int *c, complex *a, int n, int num_qam, int type);
void trellis_shaping_caf(int *c, complex *a_caf, complex *a, int n, int m, int type);
void trellis_shaping_caf2(int *c, complex *a_caf, complex *a, int n, int m);
void trellis_shaping_caf3(int *c, complex *a_caf, complex *a, int n, int m);


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
        fprintf(stderr, "%d and %d bits cannot be multiplexed to %d bits.\n", num_s, num_b, num_d);
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
void convolutional_encoding2 (int *u, int *c, int n) {
    int i;                              // ループカウンタ
    int register11, register12;         // シフトレジスタ
    int register21, register22;         // シフトレジスタ


    // シフトレジスタを初期化
    register11 = 0;
    register12 = 0;
    register21 = 0;
    register22 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[4*i] = u[2*i] ^ register12;
        c[4*i+1] = u[2*i] ^ register11 ^ register12;
        c[4*i+2] = u[2*i+1] ^ register22;
        c[4*i+3] = u[2*i+1] ^ register21 ^ register22;

        // レジスタをシフト
        register12 = register11;
        register11 = u[2*i];
        register22 = register21;
        register21 = u[2*i+1];
    }
}


// 畳み込み符号化(拘束長3)
void convolutional_encoding3 (int *u, int *c, int n) {
    int i;                              // ループカウンタ
    int register11, register12;           // シフトレジスタ
    int register21, register22;         // シフトレジスタ
    int register31, register32;         // シフトレジスタ


    // シフトレジスタを初期化
    register11 = 0;
    register12 = 0;
    register21 = 0;
    register22 = 0;
    register31 = 0;
    register32 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[6*i] = u[3*i] ^ register12;
        c[6*i+1] = u[3*i] ^ register11 ^ register12;
        c[6*i+2] = u[3*i+1] ^ register22;
        c[6*i+3] = u[3*i+1] ^ register21 ^ register22;
        c[6*i+4] = u[3*i+2] ^ register32;
        c[6*i+5] = u[3*i+2] ^ register31 ^ register32;

        // レジスタをシフト
        register12 = register11;
        register11 = u[3*i];
        register22 = register21;
        register21 = u[3*i+1];
        register32 = register31;
        register31 = u[3*i+2];
    }
}


// パリティ検査行列による復号(拘束長3)
void parity_check_decoding2 (int *c, int *u, int n) {
    int i;                                                      // ループカウンタ
    int register11, register12, register13, register14;         // シフトレジスタ
    int register21, register22, register23, register24;         // シフトレジスタ

    // シフトレジスタを初期化
    register11 = 0;
    register12 = 0;
    register13 = 0;
    register14 = 0;
    register21 = 0;
    register22 = 0;
    register23 = 0;
    register24 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        u[2*i] = c[4*i] ^ register11 ^ register12 ^ c[4*i+1] ^ register14;
        u[2*i+1] = c[4*i+2] ^ register21 ^ register22 ^ c[4*i+3] ^ register24;

        // レジスタをシフト
        register14 = register13;
        register13 = c[4*i+1];
        register12 = register11;
        register11 = c[4*i];
        register24 = register23;
        register23 = c[4*i+3];
        register22 = register21;
        register21 = c[4*i+2];
    }
}


// パリティ検査行列による復号(拘束長3)
void parity_check_decoding3 (int *c, int *u, int n) {
    int i;                                                      // ループカウンタ
    int register11, register12, register13, register14;
    int register21, register22, register23, register24;          // シフトレジスタ
    int register31, register32, register33, register34;          // シフトレジスタ

    // シフトレジスタを初期化
    register11 = 0;
    register12 = 0;
    register13 = 0;
    register14 = 0;
    register21 = 0;
    register22 = 0;
    register23 = 0;
    register24 = 0;
    register31 = 0;
    register32 = 0;
    register33 = 0;
    register34 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        u[3*i] = c[6*i] ^ register11 ^ register12 ^ c[6*i+1] ^ register14;
        u[3*i+1] = c[6*i+2] ^ register21 ^ register22 ^ c[6*i+3] ^ register24;
        u[3*i+2] = c[6*i+4] ^ register31 ^ register32 ^ c[6*i+5] ^ register34;

        // レジスタをシフト
        register14 = register13;
        register13 = c[6*i+1];
        register12 = register11;
        register11 = c[6*i];
        register24 = register23;
        register23 = c[6*i+3];
        register22 = register21;
        register21 = c[6*i+2];
        register34 = register33;
        register33 = c[6*i+5];
        register32 = register31;
        register31 = c[6*i+4];
    }
}


// 並列パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_encoding2 (int *u, int *c, int n) {
    int i;                                  // ループカウンタ
    int register1, register2;               // シフトレジスタ

    // シフトレジスタを初期化
    register1 = 0;
    register2 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[4*i] = register1;
        c[4*i+1] = u[2*i] ^ register1;

        c[4*i+2] = register2;
        c[4*i+3] = u[2*i+1] ^ register2;

        // レジスタをシフト
        register1 = u[2*i];
        register2 = u[2*i+1];
    }
}


// 並列パリティ検査行列の左逆行列による符号化(拘束長3)
void inverse_parity_check_encoding3 (int *u, int *c, int n) {
    int i;                                      // ループカウンタ
    int register1, register2, register3;        // シフトレジスタ

    // シフトレジスタを初期化
    register1 = 0;
    register2 = 0;
    register3 = 0;

    // 畳み込み
    for (i = 0; i < n; i++) {
        c[6*i] = register1;
        c[6*i+1] = u[3*i] ^ register1;

        c[6*i+2] = register2;
        c[6*i+3] = u[3*i+1] ^ register2;

        c[6*i+4] = register3;
        c[6*i+5] = u[3*i+2] ^ register3;

        // レジスタをシフト
        register1 = u[3*i];
        register2 = u[3*i+1];
        register3 = u[3*i+2];
    }
}


// コンステレーションを初期化
void construct_constellation (complex *constellation, int m, int type) {
    // コンステレーション数のチェック
    if (check_power_of_4(m) == 0 && m <= 256) {
        printf("invalid number of constellation.\n");
        exit(-1);
    }

    // タイプのチェック
    if (type != 1 && type != 2) {
        printf("invalid argument for type (type = 1 or 2).\n");
        exit(-1);
    }

    const int num_bit = (int)log2(m);
    int bits[num_bit];
    int i;

    for (i = 0; i < m; i++) {
        convert_decimal_into_binary(i, bits, num_bit);

        if (m == 4) {
            constellation[i] = map_qpsk(bits[0], bits[1]);
        } else if (m == 16) {
            if (type == 1) {
                constellation[i] = map_16qam_type1(bits[0], bits[1], bits[2], bits[3]);
            } else {
                constellation[i] = map_16qam_type2(bits[0], bits[1], bits[2], bits[3]);
            }
        } else if (m == 64) {
            if (type == 1) {
                constellation[i] = map_64qam_type1(bits[0], bits[1], bits[2], bits[3], bits[4], bits[5]);
            } else {
                constellation[i] = map_64qam_type2(bits[0], bits[1], bits[2], bits[3], bits[4], bits[5]);
            }
        } else if (m == 256) {
            if (type == 1) {
                constellation[i] = map_256qam_type1(bits[0], bits[1], bits[2], bits[3], bits[4], bits[5], bits[6], bits[7]);
            } else {
                constellation[i] = map_256qam_type2(bits[0], bits[1], bits[2], bits[3], bits[4], bits[5], bits[6], bits[7]);
            }
        }
    }
}

// トレリスシェイピングの変調
void cbts_qam_modulation (int *c, complex *a, int n, int m_m, int m_l, int type) {
    // コンステレーション数のチェック
    if (check_power_of_4(m_m) == 0 || m_m > 64) {
        printf("invalid number of constellation!!.\n");
        exit(-1);
    }

    // コンステレーション数のチェック
    if (check_power_of_4(m_l) == 0 || m_m * m_l > 256) {
        printf("invalid number of constellation.\n");
        exit(-1);
    }

    // マッピングタイプのチェック
    if (type != 1 && type != 2) {
        printf("invalid argument for type (type = 1 or 2).\n");
        exit(-1);
    }

    const int num_code_bit = (int)log2(m_m * m_l);      // 符号ビット数
    const int num_info_bit = (int)log2(m_m);            // 情報ビット数
    const double scaling_size = sqrt(m_l);              // スケーリングサイズ
    int bits[num_info_bit];                             // ビット列
    int i, j;                                           // ループカウンタ

    for (i = 0; i < n; i++) {
        for (j = 0; j < num_info_bit; j++) {
            bits[j] = c[num_code_bit * i + j];
        }

        if (m_m == 4) {
            a[i] = map_qpsk(bits[0], bits[1]) * scaling_size;
        } if (m_m == 16) {
            if (type == 1) {
                a[i] = map_16qam_type1(bits[0], bits[1], bits[2], bits[3]) * scaling_size;
            } else {
                a[i] = map_16qam_type2(bits[0], bits[1], bits[2], bits[3]) * scaling_size;
            }
        } else if (m_m == 64) {
            if (type == 1) {
                a[i] = map_64qam_type1(bits[0], bits[1], bits[2], bits[3], bits[4], bits[5]) * scaling_size;
            } else {
                a[i] = map_64qam_type2(bits[0], bits[1], bits[2], bits[3], bits[4], bits[5]) * scaling_size;
            }
        } else if (m_m == 256) {
            if (type == 1) {
                a[i] = map_256qam_type1(bits[0], bits[1], bits[2], bits[3], bits[4], bits[5], bits[6], bits[7]) * scaling_size;
            } else {
                a[i] = map_256qam_type2(bits[0], bits[1], bits[2], bits[3], bits[4], bits[5], bits[6], bits[7]) * scaling_size;
            }
        }
    }
}


// トレリスシェイピングのQAM変調
void ts_qam_modulation(int *c, complex *a, int n, int m, int type)
{
    // コンステレーション数のチェック
    if (check_power_of_4(m) == 0 || m > 256) {
        printf("invalid number of constellation.\n");
        exit(-1);
    }

    for (int i = 0; i < n; i++) {
        if (m == 4) {
            a[i] = map_qpsk(c[2*i], c[2*i+1]);
        } else if (m == 16) {
            if (type == 1) {
                a[i] = map_16qam_type1(c[4*i], c[4*i+1], c[4*i+2], c[4*i+3]);
            } else {
                a[i] = map_16qam_type2(c[4*i], c[4*i+1], c[4*i+2], c[4*i+3]);
            }
        } else if (m == 64) {
            if (type == 1) {
                a[i] = map_64qam_type1(c[6*i], c[6*i+1], c[6*i+2], c[6*i+3], c[6*i+4], c[6*i+5]);
            } else {
                a[i] = map_64qam_type2(c[6*i], c[6*i+1], c[6*i+2], c[6*i+3], c[6*i+4], c[6*i+5]);
            }
        } else if (m == 256) {
            if (type == 1) {
                a[i] = map_256qam_type1(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3], c[8*i+4], c[8*i+5], c[8*i+6], c[8*i+7]);
            } else {
                a[i] = map_256qam_type2(c[8*i], c[8*i+1], c[8*i+2], c[8*i+3], c[8*i+4], c[8*i+5], c[8*i+6], c[8*i+7]);
            }
        }
    }
}


void ts_qam_demodulation(complex *a, int *c, int n, int m, int type)
{
    const int num_bit = (int)log2(m);             // ビット数
    int i;                                        // ループカウンタ

    // コンステレーション数のチェック
    if (m != 16 && m != 64 && m != 256) {
        printf("invalid number of constellation.\n");
        exit(-1);
    }

    // タイプのチェック
    if (type != 1 && type != 2) {
        printf("invalid argument for type (type = 1 or 2).\n");
        exit(-1);
    }

    // 復調
    for (i = 0; i < n; i++) {
        if (m == 16) {
            if (type == 1) {
                unmap_16qam_type1(a[i], c + i * num_bit);
            } else {
                unmap_16qam_type2(a[i], c + i * num_bit);
            }
        } else if (m == 64) {
            if (type == 1) {
                unmap_64qam_type1(a[i], c + i * num_bit);
            } else {
                unmap_64qam_type2(a[i], c + i * num_bit);
            }
        } else if (m == 256) {
            if (type == 1) {
                unmap_256qam_type1(a[i], c + i * num_bit);
            } else {
                unmap_256qam_type2(a[i], c + i * num_bit);
            }
        }
    }
}


// トレリスシェーピング
void trellis_shaping (int *c, complex *a, int n, int m, int type)
{
    const int num_state = 4;                                                            // 状態数
    const int num_trans = 2;                                                            // 遷移数
    const int next_state[num_state][num_trans] = {{0,1}, {2,3}, {0,1}, {2,3}};          // 次状態
    const int output1[num_state][num_trans] = {{0,1}, {0,1}, {1,0}, {1,0}};             // 1ビット目の出力
    const int output2[num_state][num_trans] = {{0,1}, {1,0}, {1,0}, {0,1}};             // 2ビット目の出力
    const double infty = (double)n * 10000000;                                          // 無限
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    int tmp_index;                                                                      // 合成した信号
    static complex *constellation;                                                      // 変調信号のテーブル
    static complex **delta_table1;                                                      // δテーブル
    static double **delta_table2;                                                       // |δ^2|テーブル
    static node **nodes;                                                                // ノード
    static int memory_flag;                                                             // メモリ管理フラグ
    int sub, state, input;                                                              // ループカウンタ
    int min_size;                                                                       // メトリックを計算するサイズ
    int i, j, l;                                                                        // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (check_power_of_4(m) == 0 || m > 256) {
            printf("invalid number of constellation.\n");
            exit(-1);
        }
        // マッピングタイプのチェック
        if (type != 1 && type != 2) {
            printf("invalid argument for type (type = 1 or 2).\n");
            exit(-1);
        }

        // メモリの確保
        constellation = (complex *)malloc(m * sizeof(complex));
        nodes = (node **)malloc((n+1) * sizeof(node *));
        delta_table1 = (complex **)malloc(m * sizeof(complex *));
        delta_table2 = (double **)malloc(m * sizeof(double *));

        for (i = 0; i <= n; i++) {
            nodes[i] = (node *)malloc(num_state * sizeof(node));
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = (int *)malloc(n * sizeof(int));
                nodes[i][j].autocor = (complex *)malloc(n * sizeof(complex));
            }
        }
        for (i = 0; i < m; i++) {
            delta_table1[i] = (complex *)malloc(m * sizeof(complex));
            delta_table2[i] = (double *)malloc(m * sizeof(double));
        }

        // 信号点配置を初期化
        construct_constellation(constellation, m, type);

        // δのテーブルを用意する
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {
                delta_table1[i][j] = constellation[i] * conj(constellation[j]);
                delta_table2[i][j] = pow(cabs(delta_table1[i][j]), 2.0);
            }
        }

        // フラグをチェック
        memory_flag = 1;
    }

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0.0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // サブキャリアの位置
        sub = i - 1;

        // 窓サイズを計算
        min_size = ((i - 1) < WINDOW_SIZE) ? i - 1 : WINDOW_SIZE;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 16) {
                    tmp_index = ((c[4*sub] ^ output1[state][input]) << 3) + ((c[4*sub+1] ^ output2[state][input]) << 2) + (c[4*sub+2] << 1) + c[4*sub+3];
                } else if (m == 64) {
                    tmp_index = ((c[6*sub] ^ output1[state][input]) << 5) + ((c[6*sub+1] ^ output2[state][input]) << 4) + (c[6*sub+2] << 3) + (c[6*sub+3] << 2) + (c[6*sub+4] << 1) + c[6*sub+5];
                } else if (m == 256){
                    tmp_index = ((c[8*sub] ^ output1[state][input]) << 7) + ((c[8*sub+1] ^ output2[state][input]) << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + (c[8*sub+6] << 1) + c[8*sub+7];
                }

                // ブランチメトリックを求める
                branch_metric = 0;

                // 第2項
                for (l = 1; l < min_size; l++) {
                    branch_metric = branch_metric + 2.0 * creal(conj(nodes[i-1][state].autocor[l]) * delta_table1[tmp_index][nodes[i-1][state].a_index[i-1-l]]);
                }

                // 第3項
                if (type == 2) {
                    for (l = 1; l <= min_size; l++) {
                        branch_metric =  branch_metric + delta_table2[tmp_index][nodes[i-1][state].a_index[i-1-l]];
                    }
                }

                // パスの選択
                if (nodes[i-1][state].metric + branch_metric < nodes[i][next_state[state][input]].metric) {
                    // メトリックの更新
                    nodes[i][next_state[state][input]].metric = nodes[i-1][state].metric + branch_metric;

                    // 前状態を保存
                    nodes[i][next_state[state][input]].pre_state = state;
                    nodes[i][next_state[state][input]].a_index[i-1] = tmp_index;
                }
            }
        }

        // 次状態に情報を渡す
        for (state = 0; state < num_state; state++) {
            if (nodes[i][state].metric > infty - 1.0) continue;

            // 変調信号を更新する
            for (l = 0; l <= i-2; l++) {
                nodes[i][state].a_index[l] = nodes[i-1][nodes[i][state].pre_state].a_index[l];
            }

            // 自己相関を更新する
            for (l = 1; l <= i-2; l++) {
                nodes[i][state].autocor[l] = nodes[i-1][nodes[i][state].pre_state].autocor[l] + delta_table1[nodes[i][state].a_index[i-1]][nodes[i][state].a_index[i-1-l]];
            }
            nodes[i][state].autocor[i-1] = delta_table1[nodes[i][state].a_index[i-1]][nodes[i][state].a_index[0]];
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


// トレリスシェーピング
void trellis_shaping_caf (int *c, complex *a_caf, complex *a, int n, int m, int type) {
    const int num_state = 4;                                                            // 状態数
    const int num_trans = 2;                                                            // 遷移数
    const int next_state[num_state][num_trans] = {{0,1}, {2,3}, {0,1}, {2,3}};          // 次状態
    const int output1[num_state][num_trans] = {{0,1}, {0,1}, {1,0}, {1,0}};             // 1ビット目の出力
    const int output2[num_state][num_trans] = {{0,1}, {1,0}, {1,0}, {0,1}};             // 2ビット目の出力
    const double infty = (double)n * 10000000;                                          // 無限
    int tmp_index;                                                                      // 合成した信号
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    static complex *constellation;                                                      // 変調信号のテーブル
    static node **nodes;                                                                // ノード
    static int memory_flag;                                                             // メモリ管理フラグ
    int sub, state, input;                                                              // ループカウンタ
    int i, j, l;                                                                        // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (m != 16 && m != 64 && m != 256) {
            printf("invalid number of constellation.\n");
            exit(-1);
        }
        // マッピングタイプのチェック
        if (type != 1 && type != 2) {
            printf("invalid argument for type (type = 1 or 2).\n");
            exit(-1);
        }

        // メモリの確保
        constellation = (complex *)malloc(m * sizeof(complex));
        nodes = (node **)malloc((n+1) * sizeof(node *));

        for (i = 0; i <= n; i++) {
            nodes[i] = (node *)malloc(num_state * sizeof(node));
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = (int *)malloc(n * sizeof(int));
            }
        }

        // 信号点配置を初期化
        construct_constellation(constellation, m, type);

        // フラグをチェック
        memory_flag = 1;
    }

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 16) {
                    tmp_index = (c[4*sub] << 3) + (c[4*sub+1] << 2) + ((c[4*sub+2] ^ output1[state][input]) << 1) + (c[4*sub+3] ^ output2[state][input]);
                } else if (m == 64) {
                    tmp_index = (c[6*sub] << 5) + (c[6*sub+1] << 4) + (c[6*sub+2] << 3) + (c[6*sub+3] << 2) + ((c[6*sub+4] ^ output1[state][input]) << 1) + (c[6*sub+5] ^ output2[state][input]);
                } else if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + ((c[8*sub+6] ^ output1[state][input]) << 1) + (c[8*sub+7] ^ output2[state][input]);
                }

                // ブランチメトリックを求める
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
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


// トレリスシェーピング
void trellis_shaping_caf2 (int *c, complex *a_caf, complex *a, int n, int m) {
    const int num_state = 4;                                                            // 状態数
    const int num_trans = 2;                                                            // 遷移数
    const int next_state[num_state][num_trans] = {{0,1}, {2,3}, {0,1}, {2,3}};          // 次状態
    const int output1[num_state][num_trans] = {{0,1}, {0,1}, {1,0}, {1,0}};             // 1ビット目の出力
    const int output2[num_state][num_trans] = {{0,1}, {1,0}, {1,0}, {0,1}};             // 2ビット目の出力
    const double infty = (double)n * 10000000;                                          // 無限
    int tmp_index;                                                                      // 合成した信号
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    static complex *constellation;                                                      // 変調信号のテーブル
    static node **nodes;                                                                // ノード
    static int memory_flag;                                                             // メモリ管理フラグ
    int sub, state, input;                                                              // ループカウンタ
    int i, j, l;                                                                        // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (m != 64 && m != 256) {
            printf("invalid number of constellation.\n");
            exit(-1);
        }

        // メモリの確保
        constellation = (complex *)malloc(m * sizeof(complex));
        nodes = (node **)malloc((n+1) * sizeof(node *));

        for (i = 0; i <= n; i++) {
            nodes[i] = (node *)malloc(num_state * sizeof(node));
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = (int *)malloc(n * sizeof(int));
            }
        }

        // 信号点配置を初期化
        construct_constellation(constellation, m, 1);

        // フラグをチェック
        memory_flag = 1;
    }

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 64) {
                    tmp_index = (c[6*sub] << 5) + (c[6*sub+1] << 4) + ((c[6*sub+2] ^ output1[state][input]) << 3) + ((c[6*sub+3] ^ output2[state][input]) << 2) + (c[6*sub+4] << 1) + c[6*sub+5];
                } else if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + ((c[8*sub+4] ^ output1[state][input]) << 3) + ((c[8*sub+5] ^ output2[state][input]) << 2) + (c[8*sub+6] << 1) + c[8*sub+7];
                }

                // ブランチメトリックを求める
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
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }

    ts_qam_demodulation(a, c, n, m, 1);

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 16) {
                    tmp_index = (c[4*sub] << 3) + (c[4*sub+1] << 2) + ((c[4*sub+2] ^ output1[state][input]) << 1) + (c[4*sub+3] ^ output2[state][input]);
                } else if (m == 64) {
                    tmp_index = (c[6*sub] << 5) + (c[6*sub+1] << 4) + (c[6*sub+2] << 3) + (c[6*sub+3] << 2) + ((c[6*sub+4] ^ output1[state][input]) << 1) + (c[6*sub+5] ^ output2[state][input]);
                } else if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + ((c[8*sub+6] ^ output1[state][input]) << 1) + (c[8*sub+7] ^ output2[state][input]);
                }

                // ブランチメトリックを求める
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
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


// トレリスシェーピング
void trellis_shaping_caf3 (int *c, complex *a_caf, complex *a, int n, int m) {
    const int num_state = 4;                                                            // 状態数
    const int num_trans = 2;                                                            // 遷移数
    const int next_state[num_state][num_trans] = {{0,1}, {2,3}, {0,1}, {2,3}};          // 次状態
    const int output1[num_state][num_trans] = {{0,1}, {0,1}, {1,0}, {1,0}};             // 1ビット目の出力
    const int output2[num_state][num_trans] = {{0,1}, {1,0}, {1,0}, {0,1}};             // 2ビット目の出力
    const double infty = (double)n * 10000000;                                          // 無限
    int tmp_index;                                                                      // 合成した信号
    double branch_metric;                                                               // ブランチメトリック
    int likely_state;                                                                   // トレースバックの最尤状態
    static complex *constellation;                                                      // 変調信号のテーブル
    static node **nodes;                                                                // ノード
    static int memory_flag;                                                             // メモリ管理フラグ
    int sub, state, input;                                                              // ループカウンタ
    int i, j, l;                                                                        // ループカウンタ

    if (memory_flag == 0) {
        // コンステレーション数のチェック
        if (m != 256) {
            printf("invalid number of constellation.\n");
            exit(-1);
        }

        // メモリの確保
        constellation = (complex *)malloc(m * sizeof(complex));
        nodes = (node **)malloc((n+1) * sizeof(node *));

        for (i = 0; i <= n; i++) {
            nodes[i] = (node *)malloc(num_state * sizeof(node));
            for (j = 0; j < num_state; j++) {
                nodes[i][j].a_index = (int *)malloc(n * sizeof(int));
            }
        }

        // 信号点配置を初期化
        construct_constellation(constellation, m, 1);

        // フラグをチェック
        memory_flag = 1;
    }

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + ((c[8*sub+2] ^ output1[state][input]) << 5) + ((c[8*sub+3] ^ output2[state][input]) << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + (c[8*sub+6] << 1) + c[8*sub+7];
                }

                // ブランチメトリックを求める
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
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }

    ts_qam_demodulation(a, c, n, m, 1);

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + ((c[8*sub+4] ^ output1[state][input]) << 3) + ((c[8*sub+5] ^ output2[state][input]) << 2) + (c[8*sub+6] << 1) + c[8*sub+7];
                }

                // ブランチメトリックを求める
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
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }

    ts_qam_demodulation(a, c, n, m, 1);

    // 初期化
    for (i = 0; i <= n; i++) {
        for (j = 0; j < num_state; j++) {
            nodes[i][j].metric = infty;
        }
    }
    nodes[0][0].metric = 0;

    // ビタビアルゴリズム開始
    for (i = 1; i <= n; i++) {
        // 決定される信号位置
        sub = i - 1;

        // ノードメトリックの計算
        for (state = 0; state < num_state; state++) {
            if (nodes[i-1][state].metric > infty - 1.0) continue;

            for (input = 0; input < num_trans; input++) {
                // 仮の変調信号を求める
                if (m == 256) {
                    tmp_index = (c[8*sub] << 7) + (c[8*sub+1] << 6) + (c[8*sub+2] << 5) + (c[8*sub+3] << 4) + (c[8*sub+4] << 3) + (c[8*sub+5] << 2) + ((c[8*sub+6] ^ output1[state][input]) << 1) + (c[8*sub+7] ^ output2[state][input]);
                }

                // ブランチメトリックを求める
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
        }
    }

    // 最尤の最終状態を見つける
    likely_state = 0;
    for (state = 0; state < num_state; state++) {
        if (nodes[n][state].metric < nodes[n][likely_state].metric) {
            likely_state = state;
        }
    }

    // 最尤の信号をaに入力
    for (i = 0; i < n; i++) {
        a[i] = constellation[nodes[n][likely_state].a_index[i]];
    }
}


#endif
