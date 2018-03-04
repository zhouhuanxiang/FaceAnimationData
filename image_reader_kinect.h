#pragma once

#include <string>
#include <fstream>

#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;

class ImageReaderKinect
{
public:
	ImageReaderKinect(string p)
	{
		path = p;
	}

	void GetFrame(int idx, cv::Mat& dframe)
	{
		char str[20];

		sprintf(str, "/d%d.png", idx);
		dframe = cv::imread(path + str, cv::IMREAD_UNCHANGED);
	}

	void GetFrame(int idx, cv::Mat& iframe, cv::Mat& dframe)
	{
		char str[20];
		sprintf(str, "/i%d.png", idx);
		iframe = cv::imread(path + str, cv::IMREAD_UNCHANGED);
		/*cv::imshow("infrared", iframe);
		cv::waitKey(0);*/

		sprintf(str, "/d%d.png", idx);
		dframe = cv::imread(path + str, cv::IMREAD_UNCHANGED);
	}
private:
	string path;
};