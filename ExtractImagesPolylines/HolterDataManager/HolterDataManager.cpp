#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "HolterData.h"
#include "QRSdetector.h"
#include "msgpack.hpp"
#include <sstream>
#include <tuple>
#include "HolterDataManager.h"
#include <future>

int bytes2uint(char arr[], int length, int start = 0, bool littleEndian = true)
{
	int result = 0;
	if (littleEndian)
	{
		for (int k = 0; k < length; k++)
		{
			result += (unsigned char)arr[start + k] * (int)pow(256, k);
		}
	}
	else
	{
		for (int k = 0; k < length; k++)
		{
			result += (unsigned char)arr[start + k] * (int)pow(256, length - 1 - k);
		}
	}
	return result;
}

int bytes2int(char arr[], int length, int start = 0, bool littleEndian = true)
{
	int result = 0;
	if (littleEndian)
	{
		for (int k = 0; k < length; k++)
		{
			result += (unsigned char)arr[start + k] * (int)pow(256, k);
		}
	}
	else
	{
		for (int k = 0; k < length; k++)
		{
			result += (unsigned char)arr[start + k] * (int)pow(256, length - 1 - k);
		}
	}
	if (result >= pow(2, 8 * length - 1))
	{
		result -= (int)pow(2, 8 * length);
	}
	return result;
}

std::string bytes2string(char arr[], int length, int start=0)
{
    std::string result="";
    for (int k=0;k<length;k++)
    {
        result += arr[start+k];
    }
    return result;
}

int vector2uint(std::vector<char> vec, int length, int start = 0, bool littleEndian = true)
{
	int result = 0;
	if (littleEndian)
	{
		for (int k = 0; k < length; k++)
		{
			result += (unsigned char)vec[start + k] * (int)pow(256, k);
		}
	}
	else
	{
		for (int k = 0; k < length; k++)
		{
			result += (unsigned char)vec[start + k] * (int)pow(256, length - 1 - k);
		}
	}
	return result;
}

int vector2int(std::vector<char> vec, int length, int start = 0, bool littleEndian = true)
{
	int result = 0;
	if (littleEndian)
	{
		for (int k = 0; k < length; k++)
		{
			result += (unsigned char)vec[start + k] * (int)pow(256, k);
		}
	}
	else
	{
		for (int k = 0; k < length; k++)
		{
			result += (unsigned char)vec[start + k] * (int)pow(256, length - 1 - k);
		}
	}
	if (result >= pow(2, 8 * length - 1))
	{
		result -= (int)pow(2, 8 * length);
	}
	return result;
}

std::vector<HolterData> RemoveBlocks(std::vector<HolterData>& holterDatas, std::vector<int> ignore_blocks)
{
	std::vector<HolterData> resultHolterDatas;
	int block = 0;
	int current_block = -1;
	for (int idx=0;idx<holterDatas.size();idx++)
	{
		HolterData holterData = holterDatas[idx];
		if (std::find(ignore_blocks.begin(),ignore_blocks.end(),holterData.block)==ignore_blocks.end())
		{
			if (current_block == -1)
			{
				current_block = holterData.block;
			}
			if (holterData.block !=current_block)
			{
				block++;
				current_block = holterData.block;
			}
			holterData.block = block;
			resultHolterDatas.push_back(holterData);
		}
	}
	return resultHolterDatas;
}

std::vector<HolterData> ReadWaveData(const std::vector<std::vector<char> >& wave_data_vec, int SAMPLING_INTERVAL_MSEC, int BLOCK_BYTES
        , tm START_DATETIME, int start_block, int end_block)
{
    std::vector<HolterData> holterDatas;
    unsigned long long start_datetime =  (unsigned long long)mktime(&START_DATETIME)*1000;
    for (int block=start_block;block<end_block;block++)
    {
        std::vector<char> wave_data = wave_data_vec[block];
        int wave_start = BLOCK_LENGTH*block;
        for (int k = 0; k < BLOCK_LENGTH; k++)
        {
            HolterData holterData;
            holterData.block = block;
            holterData.timestamp_msec = start_datetime + (k + wave_start)*SAMPLING_INTERVAL_MSEC;
            holterData.ECG1 = vector2int(wave_data, BLOCK_BYTES, k*BLOCK_BYTES);
            holterData.ECG2 = vector2int(wave_data, BLOCK_BYTES, (k + BLOCK_LENGTH)*BLOCK_BYTES);
            holterData.STATUS = vector2uint(wave_data, BLOCK_BYTES, (k + 2 * BLOCK_LENGTH)*BLOCK_BYTES);
            if (k < BLOCK_LENGTH / 2)
            {
                holterData.MOTION = vector2uint(wave_data, 4, (3 * BLOCK_LENGTH) * BLOCK_BYTES);
            }
            else
            {
                holterData.MOTION = vector2uint(wave_data, 4, (3 * BLOCK_LENGTH) * BLOCK_BYTES + 4);
            }
            holterDatas.push_back(holterData);
        }
    }
    return holterDatas;
}


std::vector<HolterData> MFERParse(const char mferFilename[], std::string *outMWF_MAN, int* outMWF_IVL, std::string* outMWF_FLT)
{
    try{
        std::ifstream ifs(mferFilename, std::ios::in | std::ios::binary);
        if (!ifs) { return std::vector<HolterData>(); }

        int BLOCK_BYTES = 2;
        int SAMPLING_INTERVAL_MSEC;
        int MOTION_INTERVAL_SEC;
        int MOTION_BYTES;
        tm START_DATETIME;

        std::vector<HolterData> holterDatas;
        std::vector<char> wave_data;
        const int bufsize = 1024;
        char buffer[bufsize];
        bool isHeader = true;
        bool isLength = false;
        int tag = -1;
        int length = -1;
        int block = 0;
        int wave_start = 0;
        int reading_count = 0;
        bool isOK = false;

        std::vector<std::vector<char> > wave_data_vec;
        //Preamble
        ifs.read(buffer,1);
        tag = (unsigned char)buffer[0];
        if (tag == 16*4)
        {
            ifs.read(buffer, 1);
            length = (unsigned char)buffer[0];
            if (length != 32)
            {
                ifs.close();
                return std::vector<HolterData>();
            }

            ifs.read(buffer, length * sizeof(char));
            std::string s( reinterpret_cast<char const*>(buffer), 17 ) ;
            if (s!="MFR Long Term ECG")
            {
                ifs.close();
                return std::vector<HolterData>();
            }
        }

        while (!ifs.eof())
        {
            reading_count++;
            if (!isOK && reading_count >1000)
            {
                ifs.close();
                return std::vector<HolterData>();
            }
            if (isHeader)
            {
                ifs.read(buffer, 1);
                tag = (unsigned char)buffer[0];
                if (tag == 16 * 3 + 15)//MWF_ATT
                {
                    ifs.read(buffer, 1);
                    if ((unsigned char)buffer[0] == 3)
                    {
                        tag = 633;
                    }
                }

                if (tag == 16 * 8)//MFERファイル終端
                {
                    break;
                }

                if (tag == 16 + 14)//MWF_WAV
                {
                    ifs.read(buffer, 1);
                    if ((unsigned char)buffer[0] == 16 * 8 + 4)
                    {
                        //std::cout << std::hex << buffer[0] << std::endl;
                    }
                }

                //std::cout << "tag:"<<  tag << std::endl;
                isHeader = false;
                isLength = true;
                //std::cout << std::hex << tag << std::endl;
            }
            else if (isLength)
            {
                if (tag == 16 + 14)//MWF_WAV, 波形データ
                {
                    ifs.read(buffer, 4);
                    length = bytes2uint(buffer, 4, 0, false);
                }
                else
                {
                    ifs.read(buffer, 1);
                    length = (unsigned char)buffer[0];
                }
                isLength = false;

            }
            else
            {
                if (tag == 16 + 14)//MWF_WAV
                {
                    isOK = true;
                    wave_data = std::vector<char>();
                    while (bufsize < length)
                    {
                        ifs.read(buffer, bufsize * sizeof(char));
                        isHeader = true;
                        for (int idx = 0; idx < bufsize; idx++)
                        {
                            wave_data.push_back(buffer[idx]);
                        }
                        length -= bufsize;
                    }
                    ifs.read(buffer, length * sizeof(char));

                    for (int idx = 0; idx < length; idx++)
                    {
                        wave_data.push_back(buffer[idx]);
                    }
                    wave_data_vec.push_back(wave_data);
                    block++;
                    wave_start += BLOCK_LENGTH;
                }
                else
                {
                    ifs.read(buffer, length * sizeof(char));
                }
                if (tag == 4)//MWF_BLK
                {
                    if (BLOCK_LENGTH != bytes2uint(buffer, length))
                    {
                        //throw new InvalidDataException("BLOCK_LENGTH is not 3750");
                    }
                }
                if (tag == 10)//MWF_DTP
                {
                    if (bytes2uint(buffer, length) == 0)
                    {
                        BLOCK_BYTES = 2;
                    }
                    else
                    {
                        //throw new InvalidDataException("Invalid Data Type");
                    }
                }
                if (tag == 11)//MWF_IVL
                {
                    //ifs.read(buffer, length * sizeof(char));

                    int idx = (unsigned char)buffer[1];

                    if (idx >= pow(2, 7))
                    {
                        idx -= (int)pow(2, 8);
                    }

                    int val = bytes2uint(buffer, length - 2, 2);
                    SAMPLING_INTERVAL_MSEC = (int)(val * pow(10, idx) * 1000);
                }
                if (tag == 633)//MWF_ATT3
                {
                    //ifs.read(buffer, length * sizeof(char));

                    bool isMotion = false;
                    int i = 0;
                    while (i < length)
                    {
                        if ((unsigned int)buffer[i] == 9)
                        {
                            i++;
                            int sub_length = (unsigned int)buffer[i];
                            i++;
                            if (bytes2uint(buffer, sub_length, i) == (4 * 16 + 2) + 256 * (1 * 16 + 0))
                            {
                                isMotion = true;
                            }
                            i += sub_length;
                        }
                        if (isMotion)
                        {
                            if ((unsigned int)buffer[i] == 11)
                            {
                                i++;
                                int sub_length = (unsigned int)buffer[i];
                                i++;
                                //(unsigned int)buffer[i]=1単位
                                int idx = (unsigned int)buffer[i + 1];
                                if (idx >= pow(2, 7))
                                {
                                    idx -= (int)pow(2, 8);
                                }
                                MOTION_INTERVAL_SEC = (int)(bytes2uint(buffer, sub_length - 2, i + 2) * pow(10, idx));
                                i += sub_length;

                            }
                            if ((unsigned char)buffer[i] == 4)
                            {
                                i++;
                                int sub_length = (unsigned char)buffer[i];
                                i++;
                                int motion_samples_inBlock = bytes2uint(buffer, sub_length, i);
                                i += sub_length;
                            }

                            if ((unsigned char)buffer[i] == 10)
                            {
                                i++;
                                int sub_length = (unsigned char)buffer[i];
                                i++;
                                if ((unsigned char)buffer[i] == 6)
                                {
                                    MOTION_BYTES = 4;
                                }
                                if ((unsigned char)buffer[i] == 1)
                                {
                                    MOTION_BYTES = 2;
                                }
                                i += sub_length;
                            }
                        }
                    }

                }

                if (tag ==16*1+7)//MWF_MAN
                {
                    *outMWF_MAN = bytes2string(buffer,length,0);
                }
                if (tag ==16*1+1)//MWF_FLT
                {
                    *outMWF_FLT = bytes2string(buffer,length,0);
                }
                if (tag ==11)//MWF_IVL
                {
                    int unit = bytes2int(buffer,1,0);
                    int idx = bytes2int(buffer,1,1);
                    int num = bytes2int(buffer,length-2,2);
                    if (unit==0)//Hz
                    {
                        float heltz = num*pow(10,idx);
                        *outMWF_IVL = (int)(1000/heltz);
                    }
                    else if (unit==1)//sec
                    {
                        float sec = num*pow(10,idx);
                        *outMWF_IVL = (int)(1000*sec);
                    }
                }

                if (tag == 16 * 8 + 5)//MWF_TIM
                {
                    int year = bytes2uint(buffer, 2, 0);
                    int month = bytes2uint(buffer, 1, 2);
                    int day = bytes2uint(buffer, 1, 3);
                    int hour = bytes2uint(buffer, 1, 4);
                    int minute = bytes2uint(buffer, 1, 5);
                    int sec = bytes2uint(buffer, 1, 6);
                    START_DATETIME = { sec,minute, hour, day, month - 1, year - 1900 };
                }
                isHeader = true;
            }
        }
        ifs.close();

        int block_num = block;
        int thread_num = 4;

        int each_length = (int)std::floor(block_num/thread_num);
        int start_block = 0;
        int end_block = each_length;
        auto result0 = std::async(std::launch::async, ReadWaveData, wave_data_vec, SAMPLING_INTERVAL_MSEC, BLOCK_BYTES, START_DATETIME, start_block, end_block);
        start_block = end_block;
        end_block = end_block + each_length;
        auto result1 = std::async(std::launch::async, ReadWaveData, wave_data_vec, SAMPLING_INTERVAL_MSEC, BLOCK_BYTES, START_DATETIME, start_block, end_block);
        start_block = end_block;
        end_block = end_block + each_length;
        auto result2 = std::async(std::launch::async, ReadWaveData, wave_data_vec, SAMPLING_INTERVAL_MSEC, BLOCK_BYTES, START_DATETIME, start_block, end_block);
        start_block = end_block;
        end_block = block_num;
        auto result3 = std::async(std::launch::async, ReadWaveData, wave_data_vec, SAMPLING_INTERVAL_MSEC, BLOCK_BYTES, START_DATETIME, start_block, end_block);

        std::vector<HolterData> holterDatas0 = result0.get();
        std::vector<HolterData> holterDatas1 = result1.get();
        std::vector<HolterData> holterDatas2 = result2.get();
        std::vector<HolterData> holterDatas3 = result3.get();

        holterDatas.insert(holterDatas.end(), holterDatas0.begin(), holterDatas0.end());
        holterDatas.insert(holterDatas.end(), holterDatas1.begin(), holterDatas1.end());
        holterDatas.insert(holterDatas.end(), holterDatas2.begin(), holterDatas2.end());
        holterDatas.insert(holterDatas.end(), holterDatas3.begin(), holterDatas3.end());
        return holterDatas;
    }
    catch (char*	e)
    {
        return std::vector<HolterData>();
    }
}
std::vector<HolterData> MsgPackParse(const char inputfile[])
{
	std::ifstream fin(inputfile, std::ios::in | std::ios::binary);
	if (!fin) {
		std::cout << "ファイルが開けません" << std::endl;
		return std::vector<HolterData>();
	}

	fin.seekg(0);

	std::stringstream buffer;
	buffer << fin.rdbuf();
	fin.close();

	std::string str(buffer.str());
	std::vector<HolterData> holterDatas;
	try
	{

		msgpack::object_handle oh =
			msgpack::unpack(str.data(), str.size());

		// deserialized object is valid during the msgpack::object_handle instance is alive.
		msgpack::object deserialized = oh.get();
		deserialized.convert_if_not_nil(holterDatas);
	}
	catch(msgpack::type_error exception)
	{
		holterDatas = std::vector<HolterData>();
	}
	return holterDatas;

}
/*
convert to mV and save into csv.
*/
void WriteCsv(const std::vector<HolterData>& holterDatas, const std::string prefix, const int mod, const int remainder)
{
    const int bufsize = 1000000;
    char buf[bufsize];
    for (int block=0;block<holterDatas.size()/BLOCK_LENGTH;block++)
    {
        if (block%mod==remainder)
        {
            std::ofstream block_file;
            block_file.rdbuf()->pubsetbuf(buf,bufsize);
            block_file.open(prefix+std::to_string(block)+".csv", std::ios::out);
            for(int i=BLOCK_LENGTH*block;i<BLOCK_LENGTH*(block+1);i++)
            {
                block_file << holterDatas[i].ECG1*5/1000.0 << std::endl;
            }
            block_file.close();
        }
    }
}
std::tuple<std::vector<std::string>,int> Extract2CSV(const std::vector<HolterData>& holterDatas, const std::string prefix="test")
{
	int block=0;
	std::vector<int> ecg_data;
	std::vector<std::string> csv_files;
	int data_count = holterDatas.size();
    const int bufsize = 1000000;
    char buf[bufsize];
    const int thread_num = 4;//std::thread::hardware_concurrency();
    std::vector<std::thread> threads;

    for(int thread_id=0;thread_id<thread_num;thread_id++)
    {
        threads.emplace_back(std::thread(WriteCsv,holterDatas,prefix,thread_num,thread_id));
    }
    for(auto& thread : threads){
        thread.join();
    }
    for (int block=0;block<holterDatas.size()/BLOCK_LENGTH;block++)
    {
        csv_files.push_back(prefix+std::to_string(block)+".csv");
    }
	return std::forward_as_tuple(csv_files,data_count);
}

std::string Save2MsgPack(const std::vector<HolterData>& holterDatas, const std::string prefix="test")
{
	std::stringstream msgbuffer;
	msgpack::pack(msgbuffer, holterDatas);
	std::ofstream outfile;
	outfile.open(prefix+"_export.mp", std::ios::out | std::ios::binary);
	outfile << msgbuffer.rdbuf();
	outfile.close();

	return prefix+"_export.mp";
}

std::tuple<std::vector<int>, std::vector<int> > SaveQRStimestamps(const std::vector<HolterData>& holterDatas, const std::string outputcsv, std::vector<float>* out_filtered_wavedata)
{
	std::vector<float> data;
	std::vector<int> qrs_peaks_indices;
	size_t total_length = holterDatas.size();
	int each_length = 3750+125;
    std::vector<int> ignore_blocks;
	std::vector<int> block_heads;

    //std::vector<float> lowfreq_AmpsForFFTScoreFile; //fft結果ファイル出力用
    //std::vector<float> highfreq_AmpsForFFTScoreFile; //fft結果ファイル出力用
    //std::vector<int> blocksForFFTScoreFile; //fft結果ファイル出力用

	int block = -1;
	for (int i=0;i<holterDatas.size();i++)
	{
		HolterData holterData = holterDatas[i];
		if (block!=holterData.block)
		{
			block_heads.push_back(i);
			block = holterData.block;
		}
	}
	std::vector<float> wavedata;
	std::vector<float> filtered_wavedata;
	for (int k=0;k<holterDatas.size();k++)
	{
		HolterData holterData = holterDatas[k];
		wavedata.push_back(holterData.ECG1*5/1000.0f);
	}

    for (int block_head : block_heads)
    {
        block = holterDatas[block_head].block;
        //blocksForFFTScoreFile.push_back(block);
        data.clear();
        for (int k=block_head;k<block_head+each_length;k++)
        {
            if (k>total_length-1)
            {
                break;
            }
            HolterData holterData = holterDatas[k];
            data.push_back(holterData.ECG1*5.0f/1000);
        }

        QRSdetector qrsdetector(data);
        qrsdetector.detect_peaks();
        qrsdetector.detect_qrs_revA();
        double lowfreq_Amp = 0.0;
        double highfreq_Amp = 0.0;

        bool fft_valid = qrsdetector.validate_fft(lowfreq_Amp,highfreq_Amp);
        //lowfreq_AmpsForFFTScoreFile.push_back((float)lowfreq_Amp);
        //highfreq_AmpsForFFTScoreFile.push_back((float)highfreq_Amp);

        if ( !fft_valid || !qrsdetector.validate() || qrsdetector.qrs_peaks_indices.size() < 5) //信号が途切れたと思しきブロック
        {
            ignore_blocks.push_back(block);
        }

        if (filtered_wavedata.size()==0)
        {
            for (int k=0;k<3750+50;k++)
            {
                filtered_wavedata.push_back(qrsdetector.filtered_ecg_measurements[k]);
            }
        }
        else if (block_head==block_heads[(int)block_heads.size()-1]){
            for (int k=50;k<qrsdetector.filtered_ecg_measurements.size();k++)
            {
                filtered_wavedata.push_back(qrsdetector.filtered_ecg_measurements[k]);
            }
        }
        else
        {
            for (int k=50;k<3750+50;k++)
            {
                if(k<qrsdetector.filtered_ecg_measurements.size())
                {
                    filtered_wavedata.push_back(qrsdetector.filtered_ecg_measurements[k]);
                }
            }
        }

        for_each(qrsdetector.qrs_peaks_indices.begin(),qrsdetector.qrs_peaks_indices.end(),[&](int& idx){
            //peak correction
            std::vector<float> tmp;
            if(idx>=3)
            {
                if (idx+3<data.size())
                {
                    std::copy(data.begin()+idx-3,data.begin()+idx+3,std::back_inserter(tmp));
                }
                else{
                    std::copy(data.begin()+idx-3,data.begin()+data.size(),std::back_inserter(tmp));
                }
                int maxElementIndex = static_cast<int>(std::max_element(tmp.begin(),tmp.end()) - tmp.begin());
                idx = idx-3+maxElementIndex;
            }else{
                if (idx+3<data.size())
                {
                    std::copy(data.begin()+idx,data.begin()+idx+3,std::back_inserter(tmp));
                }
                else{
                    std::copy(data.begin()+idx,data.begin()+data.size(),std::back_inserter(tmp));
                }
                int maxElementIndex = static_cast<int>(std::max_element(tmp.begin(), tmp.end()) - tmp.begin());
                idx = maxElementIndex;
            }
            qrs_peaks_indices.push_back(block_head+idx);
        });
    }

	*out_filtered_wavedata = filtered_wavedata;

	qrs_peaks_indices = QRSdetector::shrink_indices(qrs_peaks_indices,25);
    /*
    std::string output_fft = outputcsv;
    std::ofstream fftscore_file(output_fft.replace(output_fft.size()-4,4,"_fft.csv"));
    for (size_t idx=0;idx<blocksForFFTScoreFile.size();idx++)
    {
        fftscore_file << blocksForFFTScoreFile.at(idx) << "," << lowfreq_AmpsForFFTScoreFile.at(idx) << "," << highfreq_AmpsForFFTScoreFile.at(idx) << std::endl;
    }
    fftscore_file.close();
    */

	std::ofstream qrs_file(outputcsv);
	unsigned long long start_timestamp = holterDatas[0].timestamp_msec;
	for_each(qrs_peaks_indices.begin(),qrs_peaks_indices.end(),[&](int& idx){
		qrs_file << start_timestamp+idx*8 << std::endl;
	});
	qrs_file.close();
	return std::forward_as_tuple(qrs_peaks_indices, ignore_blocks);
}
