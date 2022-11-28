#include "QRSdetector.h"

void QRSdetector::bandpass_filter(float data[], size_t data_size, float lowcut, float highcut, int signal_freq, int filter_order, float out_result[])
{
    /*
	float nyquist_freq = 0.5f*signal_freq;
	float low = lowcut / nyquist_freq;
	float high = highcut / nyquist_freq;
     */
    //pre-calculated by scipy.signal.butter
    /*
      from scipy.signal import butter, iirnotch, lfilter
      ## A low pass filter allows frequencies lower than a cut-off value
      def butter_lowpass(cutoff, fs, order=5):
      nyq = 0.5*fs
      normal_cutoff = cutoff/nyq
      b, a = butter(order, normal_cutoff, btype='low', analog=False, output='ba')
      return b, a

      highcut = 25;
      signal_freq = 125;//Hz
      order = 5;
      num,den = butter_lowpass(25, 125, 5)
      float numerator[6] = { 0.00276693, 0.01383464, 0.02766928, 0.02766928, 0.01383464, 0.00276693};
      float denominator[6] = {1.        , -2.57441532,  2.99386863, -1.85070748,  0.60066675, -0.08087089};


      float numerator[3] = { 0.28363068, 0., -0.28363068};
      float denominator[3] = {1.0, -1.43273864,  0.43273864};

      linear_filter(data, data_size, numerator, 3, denominator, 3, out_result);
      */

    //float numerator[6] ={0.81587532, 4.0793766 , 8.1587532 , 8.1587532 , 4.0793766 , 0.81587532};
    //float denominator[6] = {1. , 4.5934214 , 8.45511522, 7.79491832, 3.59890277, 0.66565254};

    float numerator[3]={ 0.94080927,  0.        , -0.94080927};
    float denominator[3]={ 1.        , -0.11838136, -0.88161855};
    linear_filter(data, data_size, numerator, 3, denominator, 3, out_result);

    /*
    float numerator[5] = {0.92117099, -3.68468397,  5.52702596, -3.68468397,  0.92117099};
    float denominator[5] = { 1.        , -3.83582554,  5.52081914, -3.53353522,  0.848556 };
    linear_filter(data, data_size, numerator, 5, denominator, 5, out_result);
     */
}

void QRSdetector::linear_filter(const float data[],size_t data_size, const float numerator[], size_t n_size, const float denominator[], size_t d_size, float out_result[])
{
    for (size_t k=0;k<data_size;k++)
    {
        float s = 0;
        for (size_t m=0;m<n_size;m++)
        {
            if (k>=m)
            {
                s+=numerator[m]*data[k-m];
            }
        }
        for (size_t n=1;n<d_size;n++)
        {
            if (k>=n)
            {
                s-=denominator[n]*out_result[k-n];
            }
        }
        out_result[k] = s/denominator[0];
    }
}

void QRSdetector::detect_peaks()
{
    float _filtered_ecg_measurements[ecg_measurements.size()];
    bandpass_filter(ecg_measurements.data(), ecg_measurements.size(), 0,60,125,1,_filtered_ecg_measurements);

    //移動平均
    float sum;
    convolve_ecg_measurements.clear();
    for (int i=0;i<ecg_measurements.size();i++) {
        sum = 0;
        int count = 0;
        for (int k=0;k<integration_window_detectPeak;k++)
        {
            if(i-k>=0)
            {
                sum+= _filtered_ecg_measurements[i-k];
                count++;
            }
        }
        convolve_ecg_measurements.push_back(sum/count);
    }
    for (int i=0;i< (int)ecg_measurements.size() - (integration_window_detectPeak - 1) / 2; i++)
    {
        convolve_ecg_measurements[i] = convolve_ecg_measurements[i + (integration_window_detectPeak - 1) / 2];
    }
    for (unsigned long i= (int)ecg_measurements.size() - (integration_window_detectPeak - 1) / 2; i < ecg_measurements.size(); i++)
    {
        convolve_ecg_measurements[i] = convolve_ecg_measurements[(int)ecg_measurements.size()-1];
    }

    for (int i=0;i<5;i++)
    {
        _filtered_ecg_measurements[i] = _filtered_ecg_measurements[5];
    }
    for (int k=0;k<ecg_measurements.size();k++)
    {
        filtered_ecg_measurements.push_back(_filtered_ecg_measurements[k]);
    }
    float squared_ecg_measurements[ecg_measurements.size()];
    for(int i=0;i<(int)ecg_measurements.size()-1;i++)
    {
        squared_ecg_measurements[i]=(filtered_ecg_measurements[i+1]-filtered_ecg_measurements[i])*(filtered_ecg_measurements[i+1]-filtered_ecg_measurements[i]);
    }
    squared_ecg_measurements[(int)ecg_measurements.size()-1] = 0;
    float _integrated_ecg_measurements[ecg_measurements.size()];
    for (int i=0;i<ecg_measurements.size();i++)
    {
        float s = 0;
        for (int k=0;k<integration_window;k++)
        {
            if(i-k>=0)
            {
                s+= squared_ecg_measurements[i-k];
            }
        }
        _integrated_ecg_measurements[i] = s;
    }
    for (int i=0;i< static_cast<int>(ecg_measurements.size()) - (integration_window - 1) / 2; i++)
    {
        _integrated_ecg_measurements[i] = _integrated_ecg_measurements[i + (integration_window - 1) / 2];
    }
    for (int i= static_cast<int>(ecg_measurements.size()) - (integration_window - 1) / 2; i < ecg_measurements.size(); i++)
    {
        _integrated_ecg_measurements[i] = _integrated_ecg_measurements[(int)ecg_measurements.size()-1];
    }

    detected_peaks_indices = findpeaks_revA(_integrated_ecg_measurements,ecg_measurements.size(), findpeaks_spacing, findpeaks_relative_limit, findpeaks_absolute_limit);
    adjust_peak_indices();

    for (float val : _integrated_ecg_measurements)
    {
        integrated_ecg_measurements.push_back(val);
    }
    detected_peaks_values = std::vector<float>();
    for (int idx : detected_peaks_indices)
    {
        detected_peaks_values.push_back(_integrated_ecg_measurements[idx]);
    }
}
void QRSdetector::adjust_peak_indices()
{
    std::vector<int> corrected_rr_indices;
    for(int idx:detected_peaks_indices)
    {
        std::vector<float> target;
        if(idx<20)
        {
            for (int k=0;k<idx+20;k++)
            {
                target.push_back(std::abs(ecg_measurements[k]-convolve_ecg_measurements[k]));
            }
            auto maxIt = std::max_element(target.begin(), target.end());
            unsigned long max_idx= std::distance(target.begin(), maxIt);
            corrected_rr_indices.push_back(max_idx);
        }
        else {
            if (idx+20<ecg_measurements.size())
            {
                for (int k = idx-20; k < idx + 20; k++) {
                    target.push_back(std::abs(ecg_measurements[k] - convolve_ecg_measurements[k]));
                }
            }
            else
            {
                for (int k = idx-20; k < ecg_measurements.size(); k++) {
                    target.push_back(std::abs(ecg_measurements[k] - convolve_ecg_measurements[k]));
                }
            }
            auto maxIt = std::max_element(target.begin(), target.end());
            auto max_idx= std::distance(target.begin(), maxIt);
            corrected_rr_indices.push_back(idx-20+max_idx);
        }
    }

    detected_peaks_indices=corrected_rr_indices;
}
/*
void QRSdetector::show_data(float data[],size_t length)
{
	for (int i=0;i<length;i++)
	{
		std::cout << data[i] << ",";
	}
	std::cout << std::endl;
}*/

bool QRSdetector::validate()
{
    if (convolve_ecg_measurements.size()==0)
    {
        return false;
    }
    const auto max_convolve = *std::max_element(std::begin(convolve_ecg_measurements), std::end(convolve_ecg_measurements));
    const auto min_convolve = *std::min_element(std::begin(convolve_ecg_measurements), std::end(convolve_ecg_measurements));
    if(max_convolve-min_convolve<min_wave_range)//wave range is too short
    {
        return false;
    }
    const auto ave_convolve = std::accumulate(std::begin(convolve_ecg_measurements), std::end(convolve_ecg_measurements), 0.0) / convolve_ecg_measurements.size();
    const auto var_convolve = std::accumulate(std::begin(convolve_ecg_measurements), std::end(convolve_ecg_measurements), 0.0, [ave_convolve](double sum, const auto& e){
        const auto temp = e - ave_convolve;
        return sum + temp * temp;
    }) / convolve_ecg_measurements.size();
    const auto std_convolve = std::sqrt(var_convolve);

    threshold_zero = calc_2ndMax(detected_peaks_values)*threshold_zero_coeff;//

    size_t cont_drv_big = 0;
    size_t cont_drv_zero = 0;
    size_t y;
    for (size_t x=0;x<(int)integrated_ecg_measurements.size()-1;x=x+2)
    {
        y = 0;
        while(y<integrated_ecg_measurements.size()-x)
        {
            if (integrated_ecg_measurements[x+y]<threshold_zero || y>threshold_cont_drv_big)
            {
                break;
            }
            y++;
        }
        if (threshold_cont_drv_big<y)
        {
            cont_drv_big = y;
            break;
        }
        else if(cont_drv_big<y)
        {
            cont_drv_big = y;
        }
    }

    for (size_t x=0;x<(int)integrated_ecg_measurements.size()-1;x=x+5)
    {
        y=0;
        while(y<integrated_ecg_measurements.size()-x) {
            if (integrated_ecg_measurements[x+y]>threshold_zero || y>threshold_cont_drv_zero)
            {
                break;
            }
            y = y + 10;
        }
        if (threshold_cont_drv_zero<y)
        {
            cont_drv_zero = y;
            break;
        }
        else if(cont_drv_zero<y)
        {
            cont_drv_zero = y;
        }
    }
    if (std_convolve>threshold_std_dv || (max_convolve-min_convolve)>threshold_max_min || cont_drv_big>threshold_cont_drv_big || cont_drv_zero>threshold_cont_drv_zero)
    {
        return false;
    }
    return true;
}

float calc_median(std::vector<float> vals)
{
    if (vals.size()==0)
    {
        return 1.0f;
    }

    std::sort(vals.begin(), vals.end());
    if (vals.size()%2==0)
    {
        return (vals.at(vals.size()/2-1) + vals.at(vals.size()/2))/2;
    }
    else
    {
        return vals.at(vals.size()/2);
    }
}

bool QRSdetector::validate_fft(double& lowfreq_Amp, double& highfreq_Amp)
{
    int BLOCK_LENGTH = ecg_measurements.size();
    std::vector<float> peak_values;
    for (size_t idx=0;idx<qrs_peaks_indices.size();idx++)
    {
        peak_values.push_back(integrated_ecg_measurements.at(qrs_peaks_indices.at(idx)));
    }
    float med_peak_val = calc_median(peak_values);

    double data[BLOCK_LENGTH];
    double F[BLOCK_LENGTH*2];
    double Amp[BLOCK_LENGTH*2];
    gsl_fft_real_wavetable  * real;
    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_workspace *w;

    for (size_t k=0;k<BLOCK_LENGTH;k++)
    {
        data[k]=ecg_measurements.at(k)/std::sqrt(med_peak_val);
    }
    w = gsl_fft_real_workspace_alloc(BLOCK_LENGTH);
    real = gsl_fft_real_wavetable_alloc(BLOCK_LENGTH);
    int stride = 1;
    gsl_fft_real_transform(data, stride, BLOCK_LENGTH, real, w);
    gsl_fft_halfcomplex_unpack(data,F,stride,BLOCK_LENGTH);

    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(w);

    for (int i=0;i<BLOCK_LENGTH*2;i++)
    {
        F[i]=F[i]/(BLOCK_LENGTH/2);
    }

    lowfreq_Amp= 0.0;
    highfreq_Amp = 0.0;
    double freq_width = signal_frequency/BLOCK_LENGTH;
    for (int i=0;i<BLOCK_LENGTH/2;i++)
    {
        Amp[i] = sqrt(F[2*i]*F[2*i]+F[2*i+1]*F[2*i+1]);
        if (freq_width*i<0.5)
        {
            lowfreq_Amp += Amp[i];
        }
    }
    double half1 = 0.0;
    double half2 = 0.0;
    for (int i=0;i<BLOCK_LENGTH/4;i++)
    {
        half1 += Amp[i];
        half2 += Amp[i+BLOCK_LENGTH/4];
    }
    if (half1>0)
    {
        highfreq_Amp = half2/half1;
    }
    else
    {
        return false;
    }

    if (lowfreq_Amp>lowfreq_threshold)
    {
        return false;
    }
    if (highfreq_Amp>highfreq_threshold)
    {
        return false;
    }
    return true;
}


float calc_2ndMax(std::vector<float> scores)
{
    size_t size = scores.size();

    if (size == 0)
    {
        return 1.0f;
    }
    else
    {
        if (size==1)
        {
            return scores[0];
        }
        std::sort(scores.begin(), scores.end());
        return scores[scores.size()-2];
    }
}

std::vector<int> QRSdetector::findpeaks_revA(const float integrated_ecg_measurements[], size_t measurements_length,  size_t spacing=1, float coeff_limit = 0, float absolute_limit = 0.01f) const
{
    float x[measurements_length+2*spacing];

    for (int i=0;i<spacing;i++)
    {
        x[i]=integrated_ecg_measurements[0]-1e-6f;
        x[measurements_length+2*spacing-1-i]=integrated_ecg_measurements[measurements_length-1]-1e-6f;
    }
    for (int i=0;i<measurements_length;i++)
    {
        x[spacing+i]=integrated_ecg_measurements[i];
    }
    bool peak_candidate[measurements_length];
    for (int i=0;i<measurements_length;i++)
    {
        peak_candidate[i] = true;
    }

    size_t start;

    float h_b[measurements_length];
    float h_c[measurements_length];
    float h_a[measurements_length];
    bool h_t1[measurements_length];
    bool h_t2[measurements_length];
    for (size_t s=0;s<spacing;s++)
    {
        start = spacing-s-1;
        for (size_t i=0;i<measurements_length;i++)
        {
            h_b[i]=x[start+i];
        }
        start = spacing;
        for (size_t i=0;i<measurements_length;i++)
        {
            h_c[i]=x[start+i];
        }
        start=spacing+s+1;
        for (size_t i=0;i<measurements_length;i++)
        {
            h_a[i]=x[start+i];
        }
        logical_lt(h_b,h_c,measurements_length,h_t1);
        logical_lt(h_a,h_c,measurements_length,h_t2);
        logical_and(h_t1,h_t2,measurements_length,h_t1);
        logical_and(peak_candidate,h_t1,measurements_length,peak_candidate);
    }
    /*
    std::vector<int> peaks;
    float peak_sum=0;
    for (int i=0;i<measurements_length;i++)
    {
        if (peak_candidate[i]&&integrated_ecg_measurements[i]>limit)
        {
            peaks.push_back(i);
            //peak_sum+=integrated_ecg_measurements[i];
        }
    }
    //return peaks;
    */
    std::vector<int> peaks2;
    for (int i=0;i<measurements_length;i++)
    {
        if (peak_candidate[i])
        {
            peaks2.push_back(i);
        }
    }

    std::vector<int> results;
    int half_window = 3;

    std::vector<float> peak_values;
    for (int k=0;k<peaks2.size();k++)
    {
        if (k<half_window)
        {
            if (k+half_window<peaks2.size())
            {
                for (int i=0;i<k+half_window;i++)
                {
                    peak_values.push_back(integrated_ecg_measurements[peaks2[i]]);
                }
            }
            else{
                peak_values.push_back(1.0f);
            }
        }
        else
        {
            if (k+half_window<peaks2.size())
            {
                for (int i=k-half_window;i<k+half_window;i++)
                {
                    peak_values.push_back(integrated_ecg_measurements[peaks2[i]]);
                }
            }
            else{
                for (int i=k-half_window;i<peaks2.size();i++)
                {
                    peak_values.push_back(integrated_ecg_measurements[peaks2[i]]);
                }
            }
        }
        float max2 = calc_2ndMax(peak_values);
        if (integrated_ecg_measurements[peaks2[k]]>max2*coeff_limit && integrated_ecg_measurements[peaks2[k]]>absolute_limit)
        {
            results.push_back(peaks2[k]);
        }
        peak_values.clear();
    }
    return results;
}
std::vector<int> QRSdetector::findpeaks(const float integrated_ecg_measurements[], size_t measurements_length,  int spacing=1, float limit = 0)
{
    float x[measurements_length+2*spacing];

    for (int i=0;i<spacing;i++)
    {
        x[i]=integrated_ecg_measurements[0]-1e-6f;
        x[measurements_length+2*spacing-1-i]=integrated_ecg_measurements[measurements_length-1]-1e-6f;
    }
    for (int i=0;i<measurements_length;i++)
    {
        x[spacing+i]=integrated_ecg_measurements[i];
    }
    bool peak_candidate[measurements_length];
    for (int i=0;i<measurements_length;i++)
    {
        peak_candidate[i] = true;
    }

    size_t start;

    float h_b[measurements_length];
    float h_c[measurements_length];
    float h_a[measurements_length];
    bool h_t1[measurements_length];
    bool h_t2[measurements_length];
    for (size_t s=0;s<spacing;s++)
    {
        start = spacing-s-1;
        for (size_t i=0;i<measurements_length;i++)
        {
            h_b[i]=x[start+i];
        }
        start = spacing;
        for (size_t i=0;i<measurements_length;i++)
        {
            h_c[i]=x[start+i];
        }
        start=spacing+s+1;
        for (size_t i=0;i<measurements_length;i++)
        {
            h_a[i]=x[start+i];
        }
        logical_lt(h_b,h_c,measurements_length,h_t1);
        logical_lt(h_a,h_c,measurements_length,h_t2);
        logical_and(h_t1,h_t2,measurements_length,h_t1);
        logical_and(peak_candidate,h_t1,measurements_length,peak_candidate);
    }

    std::vector<int> peaks;
    for (int i=0;i<measurements_length;i++)
    {
        if (peak_candidate[i]&&integrated_ecg_measurements[i]>limit)
        {
            peaks.push_back(i);
        }
    }
    return peaks;
}

void QRSdetector::detect_qrs()
{
    int last_qrs_index = 0;
    qrs_peaks_indices = std::vector<int>();
    for (int i=0;i<detected_peaks_indices.size();i++)
    {
        int detected_peak_index = detected_peaks_indices[i];
        float detected_peak_value = detected_peaks_values[i];
        if (!qrs_peaks_indices.empty())
        {
            last_qrs_index = qrs_peaks_indices[(int)qrs_peaks_indices.size()-1];
        }
        if (detected_peak_index - last_qrs_index > refractory_period || qrs_peaks_indices.empty())
        {
            if (detected_peak_value > threshold_value)
            {
                qrs_peaks_indices.push_back(detected_peak_index);
                qrs_peak_value = qrs_peak_filtering_factor * detected_peak_value + (1-qrs_peak_filtering_factor)*qrs_peak_value;
            }
            else
            {
                noise_peaks_indices.push_back(detected_peak_index);
                noise_peak_value = noise_peak_filtering_factor * detected_peak_value + (1-noise_peak_filtering_factor) * noise_peak_value;
            }
            threshold_value = noise_peak_value + qrs_noise_diff_weight * (qrs_peak_value - noise_peak_value);
        }
    }

}

void QRSdetector::detect_qrs_revA()
{
    std::vector<float> qrs_act_peaks_values;
    for (int idx:detected_peaks_indices)
    {
        qrs_act_peaks_values.push_back(std::abs(ecg_measurements[idx]-convolve_ecg_measurements[idx]));
    }
    /*
    if (!qrs_act_peaks_values.empty())
    {
        threshold_value=std::accumulate( qrs_act_peaks_values.begin(), qrs_act_peaks_values.end(), 0.0f) / static_cast<float>(qrs_act_peaks_values.size())*threshold_coeff;
    }*/
    if (qrs_act_peaks_values.size()==0)
    {
        return;
    }

    filter_peaks(qrs_act_peaks_values);
    for (int k=0;k<(int)qrs_peaks_indices.size()-1;k++)
    {
        if (qrs_peaks_indices[k+1]-qrs_peaks_indices[k]>2*signal_frequency)//RR>2.0
        {
            threshold_coeff = threshold_coeff * threshold_value_change_coeff;
            //threshold_value = threshold_value*threshold_value_change_coeff;
            filter_peaks(qrs_act_peaks_values);
            break;
        }
    }
}

void QRSdetector::filter_peaks(std::vector<float>& qrs_act_peaks_values)
{
    int last_qrs_index = 0;
    qrs_peaks_indices = std::vector<int>();
    if (qrs_act_peaks_values.size()>5)
    {
        float tmp_threshold = std::accumulate( qrs_act_peaks_values.begin(), qrs_act_peaks_values.begin()+5, 0.0f) / 5 * threshold_coeff;
        for (int i=0;i<detected_peaks_indices.size();i++)
        {
            int detected_peak_index = detected_peaks_indices[i];
            if (i>2 && i+3<qrs_act_peaks_values.size())
            {
                tmp_threshold=std::accumulate( qrs_act_peaks_values.begin()+i-2, qrs_act_peaks_values.begin()+i+3, 0.0f) / 5*threshold_coeff;
            }
            if (!qrs_peaks_indices.empty())
            {
                last_qrs_index = qrs_peaks_indices[(int)qrs_peaks_indices.size()-1];
            }
            if (detected_peak_index - last_qrs_index > refractory_period || qrs_peaks_indices.empty())
            {
                if (std::abs(qrs_act_peaks_values[i]) > std::abs(tmp_threshold))
                {
                    qrs_peaks_indices.push_back(detected_peak_index);
                }
            }
        }
    }
    else if (qrs_act_peaks_values.size()>0)
    {
        float tmp_threshold = std::accumulate( qrs_act_peaks_values.begin(), qrs_act_peaks_values.end(), 0.0f) / qrs_act_peaks_values.size() * threshold_coeff;
        for (int i=0;i<detected_peaks_indices.size();i++)
        {
            int detected_peak_index = detected_peaks_indices[i];
            if (!qrs_peaks_indices.empty())
            {
                last_qrs_index = qrs_peaks_indices[(int)qrs_peaks_indices.size()-1];
            }
            if (detected_peak_index - last_qrs_index > refractory_period || qrs_peaks_indices.empty())
            {
                if (std::abs(qrs_act_peaks_values[i]) > std::abs(tmp_threshold))
                {
                    qrs_peaks_indices.push_back(detected_peak_index);
                }
            }
        }
    }
}

void QRSdetector::logical_and(const bool array1[], const bool array2[], size_t array_length, bool out_array[])
{
    for (size_t i=0;i<array_length;i++)
    {
        out_array[i]=array1[i] && array2[i];
    }
}
void QRSdetector::logical_lt(const float array1[], const float array2[], size_t array_length, bool out_array[])
{
    for (size_t i=0;i<array_length;i++)
    {
        out_array[i] = array1[i] <= array2[i];
    }
}
std::vector<int> QRSdetector::shrink_indices(const std::vector<int>& indices, int threshold)//threshold = 25
{

    std::vector<int> result;
    int index = -100;
    for(auto idx:indices)
    {
        if (idx-index>threshold)
        {
            result.push_back(idx);
            index = idx;
        }
    }
    return result;
}
