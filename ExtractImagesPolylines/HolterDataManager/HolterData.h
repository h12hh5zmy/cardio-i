#pragma once
#include <msgpack.hpp>

class HolterData
{
public:
	int block;

	unsigned long long timestamp_msec;
	int ECG1;// unit: 5uV
	int ECG2;// unit: 5uV
	int STATUS;
	int MOTION;

	MSGPACK_DEFINE(block, timestamp_msec, ECG1, ECG2, STATUS, MOTION)
	HolterData();
	~HolterData();
};
