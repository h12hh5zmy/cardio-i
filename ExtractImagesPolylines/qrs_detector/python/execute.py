import pandas as pd
import numpy as np
import glob,pickle,time,re,MFERManager
import matplotlib.pyplot as plt
import msgpack as mp
import os, sys
import shutil

# https://github.com/c-labpl/qrs_detector
from MyQRSDetectorOffline import QRSDetectorOffline

to_plot = False

#フォルダ名
datafolder = ""


#mferの場合
mfer = "/home/cloud/Dropbox/HolterDeep/output.mwf" #datafolder + "ID60069_63h.MWF"
import glob

for mfer in glob.glob(sys.argv[1]):
    print(mfer)
    #mfer = sys.argv[1]
    mferdata = MFERManager.MFERData(mfer)
    df = mferdata.df.copy()

    #csvの場合
    #csv="af9.csv"
    #df=pd.read_csv(datafolder +csv)

    df=df.loc[:,['block','timestamp_msec','ECG1','ECG2','STATUS','MOTION']]

    #mferの場合
    df["ECG1_mV"] = df["ECG1"] * 5.0 / 1000
    df["ECG2_mV"] = df["ECG2"] * 5.0 / 1000

    #csvの場合
    #df=df.rename(columns={'ECG1':'ECG1_mV'})

    df = df.sort_values(by="timestamp_msec")
    df = df.reset_index()
    start = time.time()
    rr_indices = []
    block_std = []
    peaks_num=pd.DataFrame(columns=['block', 'peak num', 'flag(below 10)', 'flag(below 20)'])
    peak_data=[]

    #ーーーーーー定数設定ーーーーーーーーーー
    #標準偏差の閾値
    threshold_std_dv=3
    #移動平均の変化最大値の閾値
    threshold_max_min=5
    #微分の2乗移動平均の０判定閾値
    threshold_zero=0.02
    #各ブロックごとの微分の2乗移動平均の０以上が連続して続く最大回数の閾値
    threshold_cont_drv_big=200
    #各ブロックごとの微分の2乗移動平均の０が連続して続く最大回数の閾値(6秒x125=750)
    threshold_cont_drv_zero=500


    i = 0
    #保存先の削除
    figDir = "figure/"
    if os.path.exists(datafolder+figDir)==True:
        shutil.rmtree(datafolder+figDir)
    convolve=pd.DataFrame()
    start_block=df.iloc[0]["block"]
    for block in df["block"].unique():
        start_index = block * 3750
        indices = df[df["block"] == block].index.tolist()

        if block != start_block:
            indices.extend(list(range(indices[0] - 20, indices[0])))
            indices.sort()
        if block != df.iloc[df.shape[0] - 1]["block"]:
            indices.extend(list(range(indices[-1] + 1, indices[-1] + 1 + 20)))
            indices.sort()

        dfdata = df.loc[indices].sort_values(by="timestamp_msec")
        start_time = dfdata.iloc[0]["timestamp_msec"]

        ecg1_data = np.array(dfdata.loc[:, ("timestamp_msec", "ECG1_mV")]).astype(np.float64)
        qrs_detector = QRSDetectorOffline(ecg_data=ecg1_data, verbose=False,
                                          log_data=False, plot_data=False, show_plot=False)
        qrs_detector.signal_frequency = 125  # Set ECG device frequency in samples per second here.
        max_val = dfdata["ECG1_mV"].max()

        qrs_detector.detect(verbose=False, log_data=False, plot_data=to_plot, show_plot=to_plot)

        #print("block:" + str(block) + " " + str(len(qrs_detector.qrs_peaks_indices)))

        mean_convolve=pd.DataFrame(pd.Series(qrs_detector.convolve_ecg_measurements.ravel()).describe()).at['mean',0]
        std_convolve = pd.DataFrame(pd.Series(qrs_detector.convolve_ecg_measurements.ravel()).describe()).at['std', 0]
        max_convolve = pd.DataFrame(pd.Series(qrs_detector.convolve_ecg_measurements.ravel()).describe()).at['max', 0]
        min_convolve = pd.DataFrame(pd.Series(qrs_detector.convolve_ecg_measurements.ravel()).describe()).at['min', 0]
        dt_convolve_drv=pd.DataFrame(qrs_detector.integrated_ecg_measurements)
        cont_drv_big=0
        cont_drv_zero=0
        for x in range(0, len(qrs_detector.integrated_ecg_measurements) - 1, 2):
            for y in range(0, len(qrs_detector.integrated_ecg_measurements) - x, 1):
                if qrs_detector.integrated_ecg_measurements[x + y] < threshold_zero or y > threshold_cont_drv_big:
                    break
            if threshold_cont_drv_big < y:
                cont_drv_big = y
                break
            elif cont_drv_big < y:
                cont_drv_big = y

        for x in range(0, len(qrs_detector.integrated_ecg_measurements) - 1, 5):
            for y in range(0, len(qrs_detector.integrated_ecg_measurements) - x, 10):
                if qrs_detector.integrated_ecg_measurements[x + y] > threshold_zero or y > threshold_cont_drv_zero:
                    break
            if threshold_cont_drv_zero < y:
                cont_drv_zero = y
                break
            elif cont_drv_zero < y:
                cont_drv_zero = y

        # 印字用
        for k in range(len(qrs_detector.qrs_peaks_indices)):
            if dfdata.iloc[qrs_detector.qrs_peaks_indices[k]]["block"] == block:
                rr_indices.append(indices[0] + qrs_detector.qrs_peaks_indices[k])
                i += 1
        #202109追加フィルタ
        if std_convolve>threshold_std_dv or (max_convolve-min_convolve)>threshold_max_min or cont_drv_big>threshold_cont_drv_big or cont_drv_zero>threshold_cont_drv_zero:
            block_std.append(block)
        #データ確認出力用
        peak_data=pd.DataFrame([(block,len(qrs_detector.qrs_peaks_indices),mean_convolve,std_convolve,max_convolve,min_convolve,cont_drv_big,cont_drv_zero)],
                               columns=['block', 'peak num','mean','std','max','min','cont_drv_big','cont_drv_zero'])
        peaks_num=peaks_num.append(peak_data,ignore_index=True)


        # if block % 500 == 0:
        #    print(block, time.time() - start)
    for idx in rr_indices:
        print(idx)
    with open(mfer.replace(".mwf","_ignores.txt").replace(".MWF","_ignores.txt"),"w") as f:
        for block in block_std:
            f.write(str(block)+"\n")
    """
    #dataFrameにラベルづけ
    df["label"] = ""
    df.loc[rr_indices, "label"] = "R"
    for i in range(len(block_std)):
     df.loc[df["block"]==block_std[i], "label"] = "F"

    #配列整理
    values = []
    idx = 0
    for block, bgroup in df.groupby("block"):
        values.extend([idx] * bgroup.shape[0])
        idx += 1

    df["idx"] = values
    df["step"] = np.concatenate(np.array([list(range(bgroup.shape[0]))] * len(df["block"].unique())))
    dfresult = df.reindex(columns=["step", "fileid", "block", "ECG1_mV", "ECG2_mV", "label"])
    dfresult.to_pickle(mfer+"_result.pkl")
    #dfresult.to_csv(datafolder+"df.CSV")
    """
