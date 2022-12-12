#include <iostream>
#include <fstream>
#include "opencv2/plot.hpp"
#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include <cmath>

#include <msgpack.hpp>
#include <array>
#include <string>
#include <sstream>
#include <fstream>
#include "HolterData.h"
#include "HolterDataManager.h"

#include <numeric>
#define BLOCK_FREQUENCY 125
#define BLOCK_TIME_SEC 30
#define BLOCK_LENGTH (BLOCK_FREQUENCY * BLOCK_TIME_SEC)
#define PLOT_HEIGHT  783
#define OUTPUT_HEIGHT 264
#define OUTPUT_WIDTH 1264
#define CV_AA 16
#include <vector>
#include <memory>
#include <string>
#include <regex>

#include <iomanip>

#ifndef _WIN64
#include <sys/stat.h> // UNIXディレクトリ作成用
#endif
#include <direct.h>

std::array<double, BLOCK_LENGTH> linear_filter(std::array<double, BLOCK_LENGTH> data, const float numerator[], size_t n_size, const float denominator[], size_t d_size)
{
    std::array<double, BLOCK_LENGTH> out_result;
    for (size_t k=0;k<BLOCK_LENGTH;k++)
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
    return out_result;
}

std::array<double, BLOCK_LENGTH> bandpass_filter(std::array<double, BLOCK_LENGTH> data)
{
    float numerator[5] = {0.99868736, -3.99474945,  5.99212417, -3.99474945,  0.99868736};
    float denominator[5] = {1.        , -3.997373  ,  5.99212245, -3.9921259 ,  0.99737645};
    return linear_filter(data, numerator, 5, denominator, 5);
}


std::array<double, BLOCK_LENGTH>  normalize(std::array<double, BLOCK_LENGTH> data)
{
	double ave = std::accumulate(std::begin(data), std::end(data), 0.0) / data.size();
	std::transform(std::begin(data), std::end(data), std::begin(data), [ave](const double &e) {return e - ave; });
	double var = std::inner_product(std::begin(data), std::end(data), std::begin(data), 0.0) /data.size();
	std::transform(std::begin(data), std::end(data), std::begin(data), [var](const double &e) {return e /sqrt(var); });
	return data;
}

int conv_ycoordinate(double y, double ymin, double ymax)
{
	int yconv = 0;

	ymax=std::ceil(ymax);
	ymin=std::floor(ymin);

	if (ymax != ymin) {
		yconv = (int)((ymax-y)*PLOT_HEIGHT/(ymax-ymin));
	} else {
		// middle of the height
		yconv = (int)((ymax - y) + PLOT_HEIGHT / 2);
	}

	return yconv;

}

std::string make_datetime()
{
    time_t t = time(NULL);
    struct tm local;
    std::stringstream ss;

    localtime_s(&local, &t);

    ss << std::setw(4) << std::setfill('0') << local.tm_year + 1900;
    ss << std::setw(2) << std::setfill('0') << local.tm_mon + 1;  // [0-11]
    ss << std::setw(2) << std::setfill('0') << local.tm_mday;  // [1-31]
    ss << std::setw(2) << std::setfill('0') << local.tm_hour;
    ss << std::setw(2) << std::setfill('0') << local.tm_min;
    ss << std::setw(2) << std::setfill('0') << local.tm_sec;

	return ss.str();
}

int main(int argc,char *argv[])
{
	std::cout << argc<<std::endl;
	std::string out_dirname = argv[1];

	std::string out_dirname_suff = out_dirname + "_" + make_datetime();

#ifdef _WIN64
	_mkdir(out_dirname_suff.c_str());
#else
	mkdir(out_dirname_suff.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif
	for (int i=2;i<argc;i++)
	{
		std::cout << argv[i] << std::endl;
		std::string input_filename = argv[i];
		std::string out_dir_number;
		
		std::regex reg("([0-9]+)_Rac2.*");
		std::smatch match_results;
		std::string dir_name;
		if(std::regex_search(input_filename,match_results,reg))
		{
			dir_name = match_results[1].str();
			out_dir_number = out_dirname_suff + "/" + dir_name;
#ifdef _WIN64
			_mkdir(out_dir_number.c_str()); // WINDOWS系のディレクトリ作成
#else
			mkdir(out_dir_number.c_str(), S_IRWXU | S_IRWXG | S_IRWXO); // UNIX系のディレクトリ作成
#endif
		}

        //msgpack
		std::ifstream fin(input_filename, std::ios::in | std::ios::binary);
		if (!fin) {
			std::cout << "ファイルが開けません" << std::endl;
			return 1;
		}

		fin.seekg(0);

		std::stringstream buffer;
		buffer << fin.rdbuf();
		fin.close();

		std::string str(buffer.str());

		msgpack::object_handle oh =
			msgpack::unpack(str.data(), str.size());

		// deserialized object is valid during the msgpack::object_handle instance is alive.
		msgpack::object deserialized = oh.get();

		std::vector<HolterData> holterDatas;
		deserialized.convert_if_not_nil(holterDatas);

        //mfer
		/*
		std::string MWF_MAN;
		std::string MWF_FLT;
		int MWF_IVL;
		holterDatas = MFERParse(input_filename,&MWF_MAN,&MWF_IVL,&MWF_FLT);
		std::cout << MWF_MAN << std::endl;
		std::cout << MWF_FLT << std::endl;
		std::cout << MWF_IVL << std::endl;
		*/

		const int BLOCK_NUM = (int)(holterDatas.size() / BLOCK_LENGTH);


		cv::Mat xData, yData1, yData2, display;
		cv::Ptr<cv::plot::Plot2d> plot;
		xData.create(1, BLOCK_LENGTH, CV_64F);//1 Row, 100 columns, Double
		yData1.create(1, BLOCK_LENGTH, CV_64F);
		yData2.create(1, BLOCK_LENGTH, CV_64F);

		std::array<double, BLOCK_LENGTH> ecg1_data;
		std::array<double, BLOCK_LENGTH> ecg2_data;
		double ymin1,ymax1,ymin2,ymax2;

		cv::Mat plotMat;
		plotMat.create(PLOT_HEIGHT,BLOCK_LENGTH,CV_8U);

		int npt[] = {BLOCK_LENGTH};
		cv::Point pt1[1][BLOCK_LENGTH];
		const cv::Point *ppt[1] = {pt1[0]};

		for (int block = 0; block < BLOCK_NUM; block++)
		{
			for (int i = 0; i < BLOCK_LENGTH; i++)
			{
				ecg1_data[i] = (double)holterDatas[BLOCK_LENGTH * block + i].ECG1;
				ecg2_data[i] = (double)holterDatas[BLOCK_LENGTH * block + i].ECG2;
			}
			//drift cancellation
			//ecg1_data = bandpass_filter(ecg1_data);
			ecg1_data = normalize(ecg1_data);

			auto result = std::minmax_element(ecg1_data.begin(), ecg1_data.end());
			ymin1 = *result.first;
			ymax1 = *result.second;

			//ecg2_data = normalize(ecg2_data);
			//result = std::minmax_element(ecg2_data.begin(), ecg2_data.end());
			//ymin2 = *result.first;
			//ymax2 = *result.second;

			plotMat = cv::Scalar(0);
			for (int k = 0; k < BLOCK_LENGTH; k++)
			{
				pt1[0][k] = cv::Point(k,conv_ycoordinate(ecg1_data[k],ymin1,ymax1));
			}
			cv::polylines(plotMat,ppt,npt,1,false,cv::Scalar(255),2,CV_AA);

			cv::resize(plotMat,display,cv::Size(OUTPUT_WIDTH, OUTPUT_HEIGHT),0,0,cv::INTER_AREA);

			std::string output_filename1 = out_dir_number + "/" + std::to_string(block) + ".png";
			imwrite(output_filename1, display);
		}
	}
	return 0;


}
