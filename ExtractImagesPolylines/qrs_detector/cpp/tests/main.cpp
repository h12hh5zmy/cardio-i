#include "QRSdetector.h"
#include "HolterDataManager.h"

int main(int argc, char *argv[])
{
	//to HolterDataManager::SaveQRStimestamps

	//std::vector<HolterData> holterDatas = MsgPackParse("/mnt/datadisk/Dropbox/AI Cardiologist/storage_sample/550e8400-e29b-41d4-a716-446655440000/msgpack_export.mp");
  std::string outMWF_MAN;
  std::string outMWF_FLT;
  int outMWF_IVL;
  std::vector<HolterData> holterDatas = MFERParse(argv[1],&outMWF_MAN,&outMWF_IVL,&outMWF_FLT);
  //std::vector<HolterData> holterDatas = MsgPackParse("/mnt/datadisk/HolterData2/2017VT/17001/17001.mp");

	std::vector<int> qrs_peaks_indices;
	size_t total_length = holterDatas.size();
	int each_length = 3750+125;

	std::vector<int> block_heads;

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
    block = 0;
	for (int block_head : block_heads)
    {
	    std::vector<float> data;
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
        if (!qrsdetector.validate())
        {
            //std::cout <<"ignore: " << block << std::endl;
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
		block++;
	}

    std::sort(qrs_peaks_indices.begin(), qrs_peaks_indices.end() );
    qrs_peaks_indices = QRSdetector::shrink_indices(qrs_peaks_indices,25);

	int count=0;
	for (auto idx : qrs_peaks_indices)
	{
		std::cout << idx << std::endl;
		count++;
	}
	std::cout << "total:" << count << std::endl;
}
