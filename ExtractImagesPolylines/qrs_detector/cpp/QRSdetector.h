#ifndef QRSDETECTOR_H
#define QRSDETECTOR_H
#pragma once
#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <iterator>
#include <algorithm>
#include <list>
#include <numeric>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>


float calc_2ndMax(std::vector<float> scores);

class QRSdetector
{
private:
    static void logical_and(const bool array1[], const bool array2[], size_t array_length, bool out_array[]);
    static void logical_lt(const float array1[], const float array2[], size_t array_length, bool out_array[]);
    void filter_peaks(std::vector<float>& qrs_act_peaks_values);
    double highfreq_threshold = 0.2;
    double lowfreq_threshold  = 1.5;

public:
    explicit QRSdetector(std::vector<float> _ecg_measurements)
    {
        ecg_measurements = std::move(_ecg_measurements);
    }

    std::vector<float> ecg_measurements;
    std::vector<float> filtered_ecg_measurements;

    std::vector<int> detected_peaks_indices;
    std::vector<float> detected_peaks_values;
    std::vector<float> integrated_ecg_measurements;
    std::vector<int> qrs_peaks_indices;
    std::vector<int> noise_peaks_indices;
    std::vector<float> convolve_ecg_measurements;

    float signal_frequency = 125;
    float filter_lowcut = 0.0f;
    float filter_highcut = 60.0f;
    int filter_order = 1;
    //detect_peaks
    int integration_window = 15;
    int integration_window_detectPeak = 30;

    //thresholds for findpeaks_revA
    float findpeaks_absolute_limit = 0.01f;
    float findpeaks_relative_limit = 0.2f;
    float min_wave_range= 0.1f; //minimal wave range

    size_t findpeaks_spacing = 25;
    //detect_qrs_revA
    int refractory_period = 30;
    float threshold_coeff = 0.4f;

    float noise_peak_value = 0;
    float qrs_peak_value = 0;
    float threshold_value = 0;
    float threshold_value_change_coeff = 0.5;
    // parameters (for signal_frequency:125Hz)
    float qrs_peak_filtering_factor = 0.125f;
    float qrs_noise_diff_weight=0.1f;
    float noise_peak_filtering_factor = 0.1f;// 0.125;
    //
    //標準偏差の閾値
    float threshold_std_dv = 3;
    //#移動平均の変化最大値の閾値
    float threshold_max_min = 5;
    //#微分の2乗移動平均の０判定閾値
    float threshold_zero = 0.02;//validate()で微分の2乗移動平均のピーク値の10%に調整
    float threshold_zero_coeff = 0.1f;
    //#各ブロックごとの微分の2乗移動平均の0以上が連続して続く最大回数の閾値
    int threshold_cont_drv_big = 100;
    //#各ブロックごとの微分の2乗移動平均の0が連続して続く最大回数の閾値(4秒x125)
    int threshold_cont_drv_zero = 15*125;//もっと長く。15秒x125など

    void detect_peaks();
    static std::vector<int> findpeaks(const float integrated_ecg_measurements[], size_t measurements_length, int spacing, float limit);
    std::vector<int> findpeaks_revA(const float integrated_ecg_measurements[], size_t measurements_length, size_t spacing, float coeff_limit, float absolute_limit) const;
    void detect_qrs();
    void detect_qrs_revA();
    void adjust_peak_indices();
    bool validate();
    bool validate_fft(double& out_lowfreq_Amp, double& out_highfreq_Amp);
    static std::vector<int> shrink_indices(const std::vector<int>& indices, int threshold=25);
    static void bandpass_filter(float data[], size_t data_size, float lowcut, float highcut, int signal_freq, int filter_order, float out_result[]);
    static void linear_filter(const float data[],size_t data_size, const float numerator[], size_t n_size, const float denominator[], size_t d_size, float out_result[]);

};
#endif