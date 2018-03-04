#pragma once

#include "./mlibutil/mlibutil.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>

//ml::mat4d Ceres2ML(double* param) 
//{
//	double rotation_R[9];
//	ceres::AngleAxisToRotationMatrix(param, rotation_R);
//	ml::mat4d t_extrinsic;
//	for (int i = 0; i < 3; ++i) 
//	{
//		for (int j = 0; j < 3; ++j) 
//		{
//			//Note: ceres-rotationMatrix is the transpose of opencv Matrix, or ml Matrix
//			t_extrinsic(i, j) = rotation_R[j * 3 + i];
//		}
//		t_extrinsic(i, 3) = param[3 + i];
//	}
//	t_extrinsic(3, 3) = 1;
//	return t_extrinsic;
//}
//void ML2Ceres(ml::mat4d& t_extrinsic, double* param) 
//{
//	double rotation_R[9];
//	ceres::AngleAxisToRotationMatrix(param, rotation_R);
//	for (int i = 0; i < 3; ++i) 
//	{
//		for (int j = 0; j < 3; ++j) 
//		{
//			//Note: ceres-rotationMatrix is the transpose of opencv Matrix, or ml Matrix
//			rotation_R[j * 3 + i] = t_extrinsic(i, j);
//		}
//		param[3 + i] = t_extrinsic(i, 3);
//	}
//	ceres::RotationMatrixToAngleAxis(rotation_R, param);
//}

void Ceres2ML(ml::mat3d& m, ml::vec3d& v, double* param)
{
	double rotation_R[9];
	ceres::AngleAxisToRotationMatrix(param, rotation_R);
	ml::mat4d t_extrinsic;
	for (int i = 0; i < 3; ++i) 
	{
		for (int j = 0; j < 3; ++j) 
		{
			//Note: ceres-rotationMatrix is the transpose of opencv Matrix, or ml Matrix
			m(i, j) = rotation_R[j * 3 + i];
		}
		v[i] = param[3 + i];
	}
}
void ML2Ceres(ml::mat3d& m, ml::vec3d& v, double* param)
{
	double rotation_R[9];
	ceres::AngleAxisToRotationMatrix(param, rotation_R);
	for (int i = 0; i < 3; ++i) 
	{
		for (int j = 0; j < 3; ++j) 
		{
			//Note: ceres-rotationMatrix is the transpose of opencv Matrix, or ml Matrix
			rotation_R[j * 3 + i] = m(i, j);
		}
		param[3 + i] = v[i];
	}
	ceres::RotationMatrixToAngleAxis(rotation_R, param);
}

void Ceres2Eigen(Eigen::Matrix3d& rotate, Eigen::Vector3d& tranlate, double* param)
{
	double rotation_R[9];
	ceres::AngleAxisToRotationMatrix(param, rotation_R);
	//Eigen::Matrix4d t_extrinsic;
	for (int i = 0; i < 3; ++i) 
	{
		for (int j = 0; j < 3; ++j) 
		{
			//Note: ceres-rotationMatrix is the transpose of opencv Matrix, or ml Matrix
			rotate(i, j) = rotation_R[j * 3 + i];
		}
		tranlate(i) = param[3 + i];
	}
}

void Eigen2Ceres(Eigen::Matrix3d& rotate, Eigen::Vector3d& tranlate, double* param)
{
	double rotation_R[9];
	ceres::AngleAxisToRotationMatrix(param, rotation_R);
	for (int i = 0; i < 3; ++i) 
	{
		for (int j = 0; j < 3; ++j) 
		{
			//Note: ceres-rotationMatrix is the transpose of opencv Matrix, or ml Matrix
			rotation_R[j * 3 + i] = rotate(i, j);
		}
		param[3 + i] = tranlate(i);
	}
	ceres::RotationMatrixToAngleAxis(rotation_R, param);
}