// GMM.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <opencv2/opencv.hpp>
#include"opencv2/core/core.hpp"
#include"opencv2/highgui/highgui.hpp"
#include"opencv2/imgproc/imgproc.hpp"
#include"opencv2/video/background_segm.hpp"

using namespace std;
using namespace cv;
int main(int argc, char** argv)
{
	std::string videoFile = "768x576.avi";

	cv::VideoCapture capture;
	capture.open(videoFile);

	if (!capture.isOpened())
	{
		std::cout << "read video failure" << std::endl;
		return -1;
	}

	/******************************************************/
	// YOUR CODE HERE :
	//声明一个类型为BackgroundSubtractorMOG2，名称为mog的 变量，调整其中的参数并查看效果
	Ptr<BackgroundSubtractorMOG2> mog = createBackgroundSubtractorMOG2(100, 36, false);
	/******************************************************/

	cv::Mat foreground;
	cv::Mat foregroundRGB;
	cv::Mat background;

	cv::Mat frame;
	long frameNo = 0;
	capture.read(frame);
	string outputVideoPath = ".\\GMM.avi";
	VideoWriter outputVideo;
	outputVideo.open(outputVideoPath, CV_FOURCC('M', 'J', 'P', 'G'), 10.0, cvSize(frame.cols / 2, frame.rows / 2), true);
	while (capture.read(frame))
	{
		++frameNo;

		std::cout << frameNo << std::endl;
		cv::Mat ff;
		resize(frame, ff, cvSize(frame.cols / 2, frame.rows / 2));
		//保存为视频
		
		/******************************************************/
		// YOUR CODE HERE :
		//使用BackgroundSubtractorMOG2类的()运算符更新背景，找到前景
		mog->apply(ff, foreground, 0.001);


		/******************************************************/


		// 腐蚀
		cv::erode(foreground, foreground, cv::Mat());

		// 膨胀
		cv::dilate(foreground, foreground, cv::Mat());

		mog->getBackgroundImage(background);   // 返回当前背景图像

		cv::imshow("video", ff);
		cv::imshow("foreground", foreground);
		cv::imshow("background", background);
		
		cvtColor(foreground, foregroundRGB , COLOR_GRAY2RGB);
		outputVideo << foregroundRGB;

		if (cv::waitKey(25) > 0)
		{
			break;
		}
	}



	return 0;
}


