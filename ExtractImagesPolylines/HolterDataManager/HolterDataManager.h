#pragma once
#include <vector>
#include "HolterData.h"
#include <tuple>
#include <thread>
#define BLOCK_LENGTH 3750

std::vector<HolterData> MFERParse(const char mferFilename[], std::string *outMWF_MAN, int* outMWF_IVL, std::string* outMWF_FLT);
std::vector<HolterData> MsgPackParse(const char mpackFilename[]);
std::tuple<std::vector<std::string>,int> Extract2CSV(const std::vector<HolterData>& holterDatas, const std::string prefix);
std::string Save2MsgPack(const std::vector<HolterData>& holterDatas, const std::string prefix);
std::tuple<std::vector<int>, std::vector<int> > SaveQRStimestamps(const std::vector<HolterData>& holterDatas, const std::string outputcsv, std::vector<float>* out_filtered_wavedata);
std::vector<HolterData> RemoveBlocks(std::vector<HolterData>& holterDatas, std::vector<int> ignore_blocks);
void WriteCsv(const std::vector<HolterData>& holterDatas, const std::string prefix, const int mod, const int remainder);
