#pragma once
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
using namespace Eigen;

#include "parameters.h"

class TrackFitError
{
public:
	TrackFitError(int index, cv::Mat& depth, cv::Mat& normal, MatrixXd& exp_models, double fx, double fy, double cx, double cy, double xmin, double ymin, double xmax, double ymax)
		:index(index), depth(depth), normal(normal), exp_models(exp_models), fx(fx), fy(fy), cx(cx), cy(cy), xmin(xmin), ymin(ymin), xmax(xmax), ymax(ymax)
	{
	}

	template <class T>
	bool operator()(const T* const R, const T* const tr, const T* const exp_coeff, T* residuals) const {
		residuals[0] = residuals[1] = residuals[2] = residuals[3] = (T)0;

		T p1[3];
		T coeff0 = T(1.0);
		for (int i = 1; i < exp_size; i++)
			coeff0 = coeff0 - exp_coeff[i];
		for (int i = 0; i < 3; ++i){
			p1[i] = (T)exp_models(3 * index + i, 0) * coeff0;
			for (int j = 1; j < exp_size; ++j)
				p1[i] += (T)exp_models(3 * index + i, j) * exp_coeff[j];
		}
		T p2[3];
		ceres::AngleAxisRotatePoint(R, p1, p2);
		for (int i = 0; i < 3; ++i) {
			p1[i] = p2[i] + tr[i];
		}
		T x = -1.0 * p1[0] / p1[2] * fx + cx;
		T y = -1.0 * p1[1] / p1[2] * fy + cy;
		int px = *(double*)&x;
		int py = *(double*)&y;
		T wx = x - (T)px;
		T wy = y - (T)py;
		int rx = px + 1, ry = py + 1;
		if (!(px >= xmin && py >= ymin && rx <= xmax && ry <= ymax)){
			return true;
		}
		double d1 = depth.at<unsigned short>(py, px) * 1e-3f;
		double d2 = depth.at<unsigned short>(py, rx) * 1e-3f;
		double d3 = depth.at<unsigned short>(ry, px) * 1e-3f;
		double d4 = depth.at<unsigned short>(ry, rx) * 1e-3f;
		const cv::Vec3f& n1 = normal.at<cv::Vec3f>(py, px);
		const cv::Vec3f& n2 = normal.at<cv::Vec3f>(py, rx);
		const cv::Vec3f& n3 = normal.at<cv::Vec3f>(ry, px);
		const cv::Vec3f& n4 = normal.at<cv::Vec3f>(ry, rx);
		T w1 = ((T)1. - wx) * ((T)1. - wy);
		T w2 = wx * ((T)1. - wy);
		T w3 = ((T)1. - wx) * wy;
		T w4 = wx * wy;
		T d = (T)0, n = (T)0;
		if (d1)
		{
			T tmp[3];
			double dx = (px - cx) / fx;
			double dy = (py - cy) / fy;
			tmp[2] = (T)d1 /** sqrt(1 / (dx * dx + dy * dy + 1.0))*/;
			tmp[0] = -1.0 * (T)dx * tmp[2] - p1[0];
			tmp[1] = -1.0 * (T)dy * tmp[2] - p1[1];
			tmp[2] = tmp[2] - p1[2];
			residuals[0] += 5.0 * w1 * (tmp[0] * (T)n1[0] + tmp[1] * (T)n1[1] + tmp[2] * (T)n1[2]);
			residuals[1] += w1 * tmp[0];
			residuals[2] += w1 * tmp[1];
			residuals[3] += w1 * tmp[2];
		}
		if (d2)
		{
			T tmp[3];
			double dx = (rx - cx) / fx;
			double dy = (py - cy) / fy;
			tmp[2] = (T)d2 /** sqrt(1 / (dx * dx + dy * dy + 1.0))*/;
			tmp[0] = -1.0 * (T)dx * tmp[2] - p1[0];
			tmp[1] = -1.0 * (T)dy * tmp[2] - p1[1];
			tmp[2] = tmp[2] - p1[2];
			residuals[0] += 5.0 * w2 * (tmp[0] * (T)n2[0] + tmp[1] * (T)n2[1] + tmp[2] * (T)n2[2]);
			residuals[1] += w2 * tmp[0];
			residuals[2] += w2 * tmp[1];
			residuals[3] += w2 * tmp[2];
		}
		if (d3)
		{
			T tmp[3];
			double dx = (px - cx) / fx;
			double dy = (ry - cy) / fy;
			tmp[2] = (T)d3 /** sqrt(1 / (dx * dx + dy * dy + 1.0))*/;
			tmp[0] = -1.0 * (T)dx * tmp[2] - p1[0];
			tmp[1] = -1.0 * (T)dy * tmp[2] - p1[1];
			tmp[2] = tmp[2] - p1[2];
			residuals[0] += 5.0 * w3 * (tmp[0] * (T)n3[0] + tmp[1] * (T)n3[1] + tmp[2] * (T)n3[2]);
			residuals[1] += w3 * tmp[0];
			residuals[2] += w3 * tmp[1];
			residuals[3] += w3 * tmp[2];
		}
		if (d4)
		{
			T tmp[3];
			double dx = (rx - cx) / fx;
			double dy = (ry - cy) / fy;
			tmp[2] = (T)d4 /** sqrt(1 / (dx * dx + dy * dy + 1.0))*/;
			tmp[0] = -1.0 * (T)dx * tmp[2] - p1[0];
			tmp[1] = -1.0 * (T)dy * tmp[2] - p1[1];
			tmp[2] = tmp[2] - p1[2];
			residuals[0] += 5.0 * w4 * (tmp[0] * (T)n4[0] + tmp[1] * (T)n4[1] + tmp[2] * (T)n4[2]);
			residuals[1] += w4 * tmp[0];
			residuals[2] += w4 * tmp[1];
			residuals[3] += w4 * tmp[2];
		}
		//residuals[0] = (T)0;
		if (landmark.find(index) != landmark.end())
		{
			for (int i = 0; i < 4; i++)
				residuals[i] = residuals[i] * 10.0;
			residuals[0] = residuals[1] = T(0.0);
		}
		//residuals[1] = residuals[2] = residuals[3] = (T)0;
		//std::cout << *(double*)&(residuals[1]) << " " << *(double*)&(residuals[2]) << " " << *(double*)&(residuals[3]) << "\n";
		/*if (abs(*(double*)&(residuals[3])) > 1.5)
		residuals[1] = residuals[2] = residuals[3] = (T)0;*/
		//residuals[1] = residuals[2] = residuals[3] = (T)0;
		return true;
	}
	static ceres::CostFunction* Create(int index, cv::Mat& depth, cv::Mat& normal, MatrixXd& exp_models, double fx, double fy, double cx, double cy, double xmin, double ymin, double xmax, double ymax) {
		// first residual dimension, followed with parameters' dimensions
		return (new ceres::AutoDiffCostFunction<TrackFitError, 4, 3, 3, exp_size>(
			new TrackFitError(index, depth, normal, exp_models, fx, fy, cx, cy, xmin, ymin, xmax, ymax)));
	}

private:
	int index;
	cv::Mat& depth, normal;
	MatrixXd& exp_models;
	double fx, fy, cx, cy;
	double xmin, ymin, xmax, ymax;

};

#define Lamda1 0.15
#define Lamda2 0.3

class TrackSmoothError
{
public:
	TrackSmoothError(VectorXd& prev_exp_coeff_, VectorXd& prev_prev_exp_coeff_)
		:prev_exp_coeff_(prev_exp_coeff_), prev_prev_exp_coeff_(prev_prev_exp_coeff_)
	{}

	template <class T>
	bool operator()(const T* const exp_coeff, T* residuals) const
	{
		for (int i = 1; i < exp_size; i++){
			residuals[i - 1] = (exp_coeff[i] + (T)prev_prev_exp_coeff_(i) - T(2.0) * T(prev_exp_coeff_(i))) * (T)(Lamda1);
		}
		return true;
	}

	static ceres::CostFunction* Create(VectorXd& prev_exp_coeff_, VectorXd& prev_prev_exp_coeff_) {
		// first residual dimension, followed with parameters' dimensions
		return (new ceres::AutoDiffCostFunction<TrackSmoothError, exp_size - 1, exp_size>(
			new TrackSmoothError(prev_exp_coeff_, prev_prev_exp_coeff_)));
	}
private:
	VectorXd& prev_exp_coeff_, prev_prev_exp_coeff_;
};

class TrackSparseError
{
public:
	TrackSparseError()
	{}

	template <class T>
	bool operator()(const T* const exp_coeff, T* residuals) const
	{
		for (int i = 1; i < exp_size; i++)
		{
			residuals[i - 1] = ((T)exp_coeff[i]) * (T)(Lamda2);
		}
		return true;
	}

	static ceres::CostFunction* Create() {
		// first residual dimension, followed with parameters' dimensions
		return (new ceres::AutoDiffCostFunction<TrackSparseError, exp_size - 1, exp_size>(
			new TrackSparseError()));
	}
private:
};