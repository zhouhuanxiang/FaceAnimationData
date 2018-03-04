#pragma once

#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <vector>
#include <string>
#include <fstream>
#include <istream>
using namespace std;

class ImageReader
{
public:
	ImageReader(string p)
	{
		path = p;
		char filename[200];
		sprintf(filename, "%s\\data.csv", path.c_str());
		ifstream ifs;
		ifs.open(filename);
		string line;
		while (getline(ifs, line))
		{
			istringstream iss(line);
			string tmp;
			iss >> tmp >> tmp;
			//cout << tmp << "\n";
			filenames_.push_back(tmp);
		}
		ifs.close();
	}

	~ImageReader()
	{}
	
	void GetFrame(int file_index, cv::Mat& frame)
	{
		char filename[200];
		file_index %= filenames_.size();
		sprintf(filename, "%s\\data\\%s", path.c_str(), filenames_[file_index].c_str());
		frame = cv::imread(filename, cv::IMREAD_UNCHANGED);
	}

	int GetFrameSize()
	{
		return filenames_.size();
	}

private:
	string path;
	vector<string> filenames_;
};