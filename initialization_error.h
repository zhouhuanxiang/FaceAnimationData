#pragma once
#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "parameters.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
using namespace std;
using namespace Eigen;
//typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Vertices

class ProjectICPErrorPoseAndCoeff {
public:
	ProjectICPErrorPoseAndCoeff(int index, cv::Mat& depth, cv::Mat& normal, MatrixXd& pca_models, 
		double fx, double fy, double cx, double cy, double xmin, double ymin, double xmax, double ymax, bool print = false)
		:index(index), depth(depth), normal(normal), pca_models(pca_models), 
		fx(fx), fy(fy), cx(cx), cy(cy), xmin(xmin), ymin(ymin), xmax(xmax), ymax(ymax), print(print)
	{
	}

	template <class T>
	bool operator()(const T* const R, const T* const tr, const T* const pca_coeff, T* residuals) const {
		residuals[0] = residuals[1] = residuals[2] = residuals[3] = (T)0;

		T p1[3];
		for (int i = 0; i < 3; ++i){
			p1[i] = (T)pca_models(3 * index + i, 0);
			for (int j = 1; j < pca_size; ++j)
				p1[i] += (T)pca_models(3 * index + i, j) * pca_coeff[j];
		}
		T p2[3];
		ceres::AngleAxisRotatePoint(R, p1, p2);
		/*for (int i = 0; i < 3; ++i) {
			p1[i] = p2[i] + tr[i];
		}*/
		p1[0] = p2[0];
		p1[1] = p2[1];
		p1[2] = p1[2] + T(0.7);

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
		int xs[4], ys[4];
		xs[0] = px; ys[0] = py;
		xs[1] = rx; ys[1] = py;
		xs[2] = px; ys[2] = ry;
		xs[3] = rx; ys[3] = ry;
		//double d1 = depth.at<unsigned short>(py, px) * 1e-3f;
		//double d2 = depth.at<unsigned short>(py, rx) * 1e-3f;
		//double d3 = depth.at<unsigned short>(ry, px) * 1e-3f;
		//double d4 = depth.at<unsigned short>(ry, rx) * 1e-3f;
		//const cv::Vec3f& n1 = normal.at<cv::Vec3f>(py, px);
		//const cv::Vec3f& n2 = normal.at<cv::Vec3f>(py, rx);
		//const cv::Vec3f& n3 = normal.at<cv::Vec3f>(ry, px);
		//const cv::Vec3f& n4 = normal.at<cv::Vec3f>(ry, rx);
		T ws[4];
		ws[0] = ((T)1. - wx) * ((T)1. - wy);
		ws[1] = wx * ((T)1. - wy);
		ws[2] = ((T)1. - wx) * wy;
		ws[3] = wx * wy;
		T dv = T(0);
		T weight = T(0);
		T nv[3];
		nv[0] = nv[1] = nv[2] = T(0);
		for (int i = 0; i < 4; i++){
			double di = depth.at<unsigned short>(ys[i], xs[i]) * 1e-3f /** sqrt(1 / (dx * dx + dy * dy + 1.0))*/;
			cv::Vec3f ni = normal.at<cv::Vec3f>(ys[i], xs[i]);
			if (di > 0 && ni != cv::Vec3f(0, 0, 0)){
				dv += T(di) * ws[i];
				weight += ws[i];
				nv[0] += ws[i] * (T)ni[0];
				nv[1] += ws[i] * (T)ni[1];
				nv[2] += ws[i] * (T)ni[2];
			}
			//else{
			//	weight = T(0);
			//	break;
			//}
		}
		if (weight != T(0) && std::abs(*(double*)&(dv)-*(double*)&(p1[2])) < 0.008){
			residuals[0] = nv[2] * (dv - p1[2])
				+ nv[0] * (dv - p1[2]) * p1[0] / p1[2]
				+ nv[1] * (dv - p1[2]) * p1[1] / p1[2];

		}
		if (print){
			std::cout << *(double*)&(nv[2]) << " "
				<< *(double*)&(dv) << " "
				<< *(double*)&(p1[2]) << " "
				<< *(double*)&(residuals[0]) << "\n";
		}

		//if (d1)
		//{
		//	T tmp[3];
		//	double dx = (px - cx) / fx;
		//	double dy = (py - cy) / fy;
		//	tmp[2] = (T)d1 /** sqrt(1 / (dx * dx + dy * dy + 1.0))*/;
		//	tmp[0] = -1.0 * (T)dx * tmp[2] - p1[0];
		//	tmp[1] = -1.0 * (T)dy * tmp[2] - p1[1];
		//	tmp[2] = tmp[2] - p1[2];
		//	residuals[0] += 2.111 * w1 * (tmp[0] * (T)n1[0] + tmp[1] * (T)n1[1] + tmp[2] * (T)n1[2]);
		//	residuals[1] += w1 * tmp[0];
		//	residuals[2] += w1 * tmp[1];
		//	residuals[3] += w1 * tmp[2];
		//}
		//if (d2)
		//{
		//	T tmp[3];
		//	double dx = (rx - cx) / fx;
		//	double dy = (py - cy) / fy;
		//	tmp[2] = (T)d2 /** sqrt(1 / (dx * dx + dy * dy + 1.0))*/;
		//	tmp[0] = -1.0 * (T)dx * tmp[2] - p1[0];
		//	tmp[1] = -1.0 * (T)dy * tmp[2] - p1[1];
		//	tmp[2] = tmp[2] - p1[2];
		//	residuals[0] += 2.111 * w2 * (tmp[0] * (T)n2[0] + tmp[1] * (T)n2[1] + tmp[2] * (T)n2[2]);
		//	residuals[1] += w2 * tmp[0];
		//	residuals[2] += w2 * tmp[1];
		//	residuals[3] += w2 * tmp[2];
		//}
		//if (d3)
		//{
		//	T tmp[3];
		//	double dx = (px - cx) / fx;
		//	double dy = (ry - cy) / fy;
		//	tmp[2] = (T)d3 /** sqrt(1 / (dx * dx + dy * dy + 1.0))*/;
		//	tmp[0] = -1.0 * (T)dx * tmp[2] - p1[0];
		//	tmp[1] = -1.0 * (T)dy * tmp[2] - p1[1];
		//	tmp[2] = tmp[2] - p1[2];
		//	residuals[0] += 2.111 * w3 * (tmp[0] * (T)n3[0] + tmp[1] * (T)n3[1] + tmp[2] * (T)n3[2]);
		//	residuals[1] += w3 * tmp[0];
		//	residuals[2] += w3 * tmp[1];
		//	residuals[3] += w3 * tmp[2];
		//}
		//if (d4)
		//{
		//	T tmp[3];
		//	double dx = (rx - cx) / fx;
		//	double dy = (ry - cy) / fy;
		//	tmp[2] = (T)d4 /** sqrt(1 / (dx * dx + dy * dy + 1.0))*/;
		//	tmp[0] = -1.0 * (T)dx * tmp[2] - p1[0];
		//	tmp[1] = -1.0 * (T)dy * tmp[2] - p1[1];
		//	tmp[2] = tmp[2] - p1[2];
		//	residuals[0] += 2.111 * w4 * (tmp[0] * (T)n4[0] + tmp[1] * (T)n4[1] + tmp[2] * (T)n4[2]);
		//	residuals[1] += w4 * tmp[0];
		//	residuals[2] += w4 * tmp[1];
		//	residuals[3] += w4 * tmp[2];
		//}
		if (landmark.find(index) != landmark.end())
		{
			for (int i = 0; i < 4; i++)
				residuals[i] = residuals[i] * 3.0;
		}
		//std::cout << *(double*)&(residuals[0]) << "\n";


		//residuals[1] = residuals[2] = residuals[3] = (T)0;
		//std::cout << *(double*)&(residuals[1]) << " " << *(double*)&(residuals[2]) << " " << *(double*)&(residuals[3]) << "\n";
		/*if (abs(*(double*)&(residuals[3])) > 1.5)
			residuals[1] = residuals[2] = residuals[3] = (T)0;*/
		return true;
	}
	static ceres::CostFunction* Create(int index, cv::Mat& depth, cv::Mat& normal, MatrixXd& pca_models, 
		double fx, double fy, double cx, double cy, double xmin, double ymin, double xmax, double ymax, bool print = false) {
		// first residual dimension, followed with parameters' dimensions
		return (new ceres::AutoDiffCostFunction<ProjectICPErrorPoseAndCoeff, 4, 3, 3, pca_size>(
			new ProjectICPErrorPoseAndCoeff(index, depth, normal, pca_models, fx, fy, cx, cy, xmin, ymin, xmax, ymax, print)));
	}

private:
	int index;
	cv::Mat& depth, normal;
	MatrixXd& pca_models;
	double fx, fy, cx, cy;
	double xmin, ymin, xmax, ymax;
	bool print;
};

class ProjectICPErrorCoeff {
public:
	ProjectICPErrorCoeff(int index, cv::Mat& depth, cv::Mat& normal, MatrixXd& pca_models, double fx, double fy, double cx, double cy, double xmin, double ymin, double xmax, double ymax, double* param)
		:index(index), depth(depth), normal(normal), pca_models(pca_models), fx(fx), fy(fy), cx(cx), cy(cy), xmin(xmin), ymin(ymin), xmax(xmax), ymax(ymax), param(param)
	{
	}

	template <class T>
	bool operator()(const T* const pca_coeff, T* residuals) const {
		residuals[0] = residuals[1] = residuals[2] = residuals[3] = (T)0;

		T p1[3];
		for (int i = 0; i < 3; ++i){
			p1[i] = (T)pca_models(3 * index + i, 0);
			for (int j = 1; j < pca_size; ++j)
				p1[i] += (T)pca_models(3 * index + i, j) * pca_coeff[j];
		}
		T p2[3];
		T R[3], tr[3];
		R[0] = (T)param[0], R[1] = (T)param[1], R[2] = (T)param[2];
		tr[0] = (T)param[3], tr[1] = (T)param[4], tr[2] = (T)param[5];
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
			residuals[0] += 2.111 * w1 * (tmp[0] * (T)n1[0] + tmp[1] * (T)n1[1] + tmp[2] * (T)n1[2]);
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
			residuals[0] += 2.111 * w2 * (tmp[0] * (T)n2[0] + tmp[1] * (T)n2[1] + tmp[2] * (T)n2[2]);
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
			residuals[0] += 2.111 * w3 * (tmp[0] * (T)n3[0] + tmp[1] * (T)n3[1] + tmp[2] * (T)n3[2]);
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
			residuals[0] += 2.111 * w4 * (tmp[0] * (T)n4[0] + tmp[1] * (T)n4[1] + tmp[2] * (T)n4[2]);
			residuals[1] += w4 * tmp[0];
			residuals[2] += w4 * tmp[1];
			residuals[3] += w4 * tmp[2];
		}
		//residuals[0] = (T)0;
		if (landmark.find(index) != landmark.end())
		{
			for (int i = 0; i < 4; i++)
				residuals[i] = residuals[i] * 3.0;
		}
		//residuals[1] = residuals[2] = residuals[3] = (T)0;
		//std::cout << *(double*)&(residuals[1]) << " " << *(double*)&(residuals[2]) << " " << *(double*)&(residuals[3]) << "\n";
		//std::cout << (residuals[1]) << "\n" << (residuals[2]) << "\n" << (residuals[3]) << "\n\n";
		/*if (abs(*(double*)&(residuals[3])) > 1.5)
		residuals[1] = residuals[2] = residuals[3] = (T)0;*/
		//residuals[1] = residuals[2] = residuals[3] = (T)0;
		return true;
	}
	static ceres::CostFunction* Create(int index, cv::Mat& depth, cv::Mat& normal, MatrixXd& pca_models, double fx, double fy, double cx, double cy, double xmin, double ymin, double xmax, double ymax, double* param) {
		// first residual dimension, followed with parameters' dimensions
		return (new ceres::AutoDiffCostFunction<ProjectICPErrorCoeff, 4, pca_size>(
			new ProjectICPErrorCoeff(index, depth, normal, pca_models, fx, fy, cx, cy, xmin, ymin, xmax, ymax, param)));
	}

private:
	int index;
	cv::Mat& depth, normal;
	MatrixXd& pca_models;
	double fx, fy, cx, cy;
	double xmin, ymin, xmax, ymax;
	double* param;
};

class Regulation1
{
public:
	Regulation1(VectorXd& pca_weights)
		:pca_weights(pca_weights)
	{
	}

	template <class T>
	bool operator()(const T* const pca_coeff, T* residuals) const
	{
		for (int i = 1; i < pca_size; i++){
			residuals[i - 1] = ((T)pca_coeff[i]) * pca_weights(i);
			//std::cout << residuals[i] << "\n";
		}
		//std::cout << "\n";
		return true;
	}

	static ceres::CostFunction* Create(VectorXd& pca_weights) {
		// first residual dimension, followed with parameters' dimensions
		return (new ceres::AutoDiffCostFunction<Regulation1, pca_size - 1, pca_size>(
			new Regulation1(pca_weights)));
	}
private:
	VectorXd& pca_weights;
};