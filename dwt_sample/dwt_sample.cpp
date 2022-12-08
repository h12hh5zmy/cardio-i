#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "wavelib.h"
#include "HolterDataManager.h"
#include "HolterData.h"

int main(int argc, char *argv[])
{
    std::vector<HolterData> holterDatas = MsgPackParse(argv[1]);
    int N = 3750;//1ブロックのフレーム数
    int block_num = int(holterDatas.size()/3750);//holterDatasブロック数

    wave_object obj;
    wt_object wt;
    double *inp,*out,*diff;
    int i,J;

    char *name = "sym4";
    obj = wave_init(name);// Initialize the wavelet
    inp = (double *) malloc(sizeof(double) * N);
    out = (double *) malloc(sizeof(double) * N);

    for (int block=0;block<block_num;block++) {
        for (i = 0; i < N; ++i) {
            inp[i] = holterDatas.at(i + block * N).ECG1 * 5.0 / 1000;
        }

        J = 9;//9レベルに離散ウェーブレット変換する

        wt = wt_init(obj, "dwt", N, J);// Initialize the wavelet transform object
        setDWTExtension(wt, "sym");// Options are "per" and "sym". Symmetric is the default option
        setWTConv(wt, "direct");

        dwt(wt, inp);// Perform DWT

        //DWT output can be accessed using wt->output vector. Use wt_summary to find out how to extract appx and detail coefficients
        // wt->outputに14 14 36 65 123 240 474 942 1878と直列にレベルごとのwavelet係数が入っている
        for (i = 0; i < wt->outlength; ++i) {
            //低周波除去の場合
            if (i < 14 + 14)// 14 14 36 65 123 240 474 942 1878
            {
                wt->output[i] = 0;
            }
            /*
            //高周波除去の場合
            if (i > 14 + 14 + 36 + 65 + 123 + 240 + 474 + 942)
            {
                wt->output[i] = 0;
            }*/
        }
        
        idwt(wt, out);// Perform IDWT (if needed)
        
        std::ofstream output("reconstruction_" + std::to_string(block) + ".csv");

        for (i = 0; i < wt->siglength; ++i) {
            output << out[i] << std::endl;
            //diff[i] = out[i] - inp[i];
        }
        output.close();
    }
    wave_free(obj);
    wt_free(wt);

    free(inp);
    free(out);
}


