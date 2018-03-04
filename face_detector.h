#pragma once
#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

class FaceDetector
{
public:
	FaceDetector()
	{
		cv::String face_cascade_name = "../haarcascade_frontalface_alt.xml";
		if (!face_cascade.load(face_cascade_name)){ 
			printf("--(!)Error loading\n"); 
		};
	}

	bool Detect(cv::Mat& frame, cv::Rect& rect)
	{
		cv::Mat frame_gray = frame;
		//cv::equalizeHist(frame, frame_gray);//
		std::vector<cv::Rect> faces;
		face_cascade.detectMultiScale(frame, faces, 1.1, 3, CV_HAAR_FIND_BIGGEST_OBJECT | CV_HAAR_DO_ROUGH_SEARCH);

		/*rectangle(frame_gray, faces[0], cv::Scalar(255, 0, 255), 2, 8, 0);
		cv::imshow("test", frame_gray);
		cv::waitKey(0);*/

		if (faces.size())
		{
			rect = faces[0];
			return true;
		}
		else
			return false;
	}
private:
	cv::CascadeClassifier face_cascade;
};
