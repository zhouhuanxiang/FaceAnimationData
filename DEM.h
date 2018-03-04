#pragma once

#include "parameters.h"
#include "model_reader.h"
#include "error_func.h"
#include "initialization_error.h"
#include "track_error.h"
#include "rigid_track_error.h"
#include "face_detector.h"
#include "./imgproc/filter.h"
#include "mask_reader.h"

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>
using namespace Eigen;

struct Camera
{
	double fx = 580;
	double fy = 580;
	double cx = 320;
	double cy = 240;
};

// kinect
//struct Camera
//{
//	double fx = 365.427002;
//	double fy = 365.427002;
//	double cx = 255.713501;
//	double cy = 208.248596;
//};

class AdaptiveDynamicExpressionModel
{
public:
	AdaptiveDynamicExpressionModel()
	{
		exp_models_.resize(3 * vertex_size, exp_size);
		exp_transforms_.resize(exp_size);
		
		pca_coeff_.resize(pca_size);
		pca_coeff_.setZero();
		pca_coeff_(0) = 1;
		exp_coeff_.resize(exp_size); // -1
		exp_coeff_.setZero();
		prev_prev_exp_coeff_ = prev_exp_coeff_ = exp_coeff_;
		lap_coeff_.resize(lap_size, pca_size);
		lap_coeff_.setZero();

		//model_.resize(vertex_size);
		ml::MeshIOd::loadFromOBJ(Data_Input_Dir + "pca_obj_simplified2/pca/0.obj", mesh_);
		mesh_.m_Colors.resize(vertex_size, ml::vec4d(0.0, 0.0, 0.0, 1.0));

		//mask_reader.Test();
		mask_reader.Read(vertex_size);
		//int count = 0;
		//for (int i = 0; i < vertex_size; i++){
		//	if (mask_reader.face_mask[i])
		//		count++;
		//}
		//std::cout << count << "\n";
		rigid_rotate_.setIdentity();

		ModelReader mreader(pca_models_, exp_transforms_, pca_weights_);
		for (int i = 1; i < pca_size; i++){
			pca_models_.col(i) = pca_models_.col(i) - pca_models_.col(0);
		}
		pca_weights_.setZero();
		for (int i = 1; i < pca_size; i++){
			for (int v = 0; v < vertex_size; v++){
				if (mask_reader.face_mask[v]){
					for (int j = 0; j < 3; j++){
						pca_weights_(i) += pca_models_(3 * v + j, i) * pca_models_(3 * v + j, i);
					}
				}
			}
			pca_weights_(i) = 0.001 * sqrt(0.5 / pca_weights_(i));
		}
		//std::cout << pca_weights_;
		//system("pause");

		D_Refine.resize(pca_size - 1, pca_size - 1);
		std::vector<Tripletd> tris;
		for (int i = 0; i < pca_size - 1; i++){
			tris.push_back(Tripletd(i, i, pca_weights_(i + 1)));
		}
		D_Refine.setFromTriplets(tris.begin(), tris.end());
		M_Refine.resize(pca_size - 1, pca_size - 1);
		M_Refine.setZero();
		Y_Refine.resize(pca_size - 1);
		Y_Refine.setZero();
		S_Refine = 0;

		// pca bbox
		/*double xmin, ymin, zmin;
		double xmax, ymax, zmax;
		xmin = ymin = zmin = 10000;
		xmax = ymax = zmax = -10000;
		for (int i = 0; i < vertex_size; i++)
		{
			xmin = std::min(xmin, pca_models_[0](0, i));
			xmax = std::max(xmax, pca_models_[0](0, i));
			ymin = std::min(ymin, pca_models_[0](1, i));
			ymax = std::max(ymax, pca_models_[0](1, i));
			zmin = std::min(zmin, pca_models_[0](2, i));
			zmax = std::max(zmax, pca_models_[0](2, i));
		}
		double xw = xmax - xmin;
		double xo = (xmin + xmax) / 2;
		double yw = ymax - ymin;
		double yo = (ymin + ymax) / 2;
		double zw = zmax - zmin;
		double zo = (zmin + zmax) / 2;
		std::cout << xw << xo << yw << yo << zw << zo;*/
	}

	
	Vector3d Point2d_2_Point3d(Vector2d p2, int depth)
	{
		if (depth == 0)
			return Vector3d(0, 0, 0);
		//double x_z = (p2[0] - camera_.cx) / camera_.fx;
		//double y_z = (p2[1] - camera_.cy) / camera_.fy;
		//ml::vec3d p3;
		//p3[2] = depth * 0.001 * sqrt(1 / (x_z * x_z + y_z * y_z + 1));
		//p3[0] = -x_z * p3[2];
		//p3[1] = -y_z * p3[2];

		Vector3d p3;
		p3(2) = depth / 1000.0;
		p3(0) = -1 * (p2(0) - camera_.cx) * p3(2) / camera_.fx;
		p3(1) = -1 * (p2(1) - camera_.cy) * p3(2) / camera_.fy;

		return p3;
	}

	template <class T2, class T3 >
	T2 Point3d_2_Point2d(T3 p3)
	{
		if (p3.z == 0)
			return T2(0, 0);
		T2 p2;
		p2.x = p3.x / p3.z * (-1) * camera_.fx + camera_.cx;
		p2.y = p3.y / p3.z * (-1) * camera_.fy + camera_.cy;

		return p2;
	}

	bool UpdateFaceRegion(cv::Mat& infrared_frame)
	{
		if (!fd.Detect(infrared_frame, face_rect_))
			return false;
		return true;
	}

	bool InitRigidMotion(cv::Mat& infrared_frame, cv::Mat& depth_frame)	
	{
		if (!UpdateFaceRegion(infrared_frame))
			return false;

		Vector2d p2(face_rect_.x + face_rect_.width / 2, face_rect_.y + face_rect_.height / 2);
		double depth = depth_frame.at<ushort>(p2[1], p2[0]);// at(y, x)
		Vector3d p3 = Point2d_2_Point3d(p2, depth);
		rigid_translate_ = p3 - pca_center;

		return true;
	}

	void InitializePhase1(cv::Mat& frame, cv::Mat& normal, double* param)
	{
		ceres::Problem problem;
		ceres::Solver::Options options;
		options.linear_solver_type = ceres::SPARSE_SCHUR;
		options.minimizer_progress_to_stdout = false;
		options.max_num_iterations = 20;
		options.num_threads = 1;
		ceres::LossFunctionWrapper* loss_function_wrapper = new ceres::LossFunctionWrapper(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
		for (int residual_count = 0; residual_count < vertex_size; residual_count++)
		{
			if (!mask_reader.face_mask[residual_count])
				continue;
			//auto error = ProjectICPErrorPoseAndCoeff(residual_count, frame, normal, pca_models_, camera_.fx, camera_.fy, camera_.cx, camera_.cy, face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height);
			//double* residual = new double[4];
			//error(param, param + 3, pca_coeff_.data(), residual);
			//model_attribute_(0, residual_count) = residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3];
			////if (model_attribute_(0, residual_count) > 0)
			////	model_attribute_(0, residual_count) = 0.0255;
			//delete[] residual;

			problem.AddResidualBlock(
				ProjectICPErrorPoseAndCoeff::Create(residual_count, frame, normal, pca_models_, camera_.fx, camera_.fy, camera_.cx, camera_.cy, face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height),
				loss_function_wrapper, param, param + 3, pca_coeff_.data()
				);
		}

		problem.AddResidualBlock(
			Regulation1::Create(pca_weights_),
			0, pca_coeff_.data()
			);

		ceres::Solver::Summary summary;
		ceres::Solve(options, &problem, &summary);
	}

	void InitializePhase2(cv::Mat& frame, cv::Mat& normal, double* param)
	{
		ceres::Problem problem;
		ceres::Solver::Options options;
		options.linear_solver_type = ceres::SPARSE_SCHUR;
		options.minimizer_progress_to_stdout = false;
		options.max_num_iterations = 20;
		options.num_threads = 1;
		ceres::LossFunctionWrapper* loss_function_wrapper = new ceres::LossFunctionWrapper(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
		for (int residual_count = 0; residual_count < vertex_size; residual_count++)
		{
			if (!mask_reader.face_mask[residual_count])
				continue;

			auto error = ProjectICPErrorCoeff(residual_count, frame, normal, pca_models_, camera_.fx, camera_.fy, camera_.cx, camera_.cy, face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height, param);
			double* residual = new double[4];
			error(pca_coeff_.data(), residual);
			double value = residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3];
			delete[] residual;

			if (value > 0.000005)
				continue;

			problem.AddResidualBlock(	
				ProjectICPErrorCoeff::Create(residual_count, frame, normal, pca_models_, camera_.fx, camera_.fy, camera_.cx, camera_.cy, face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height, param),
				loss_function_wrapper, pca_coeff_.data()
				);
		}

		problem.AddResidualBlock(
			Regulation1::Create(pca_weights_),
			0, pca_coeff_.data()
			);

		ceres::Solver::Summary summary;
		ceres::Solve(options, &problem, &summary);
	}

	void TwoPhaseInitialize(cv::Mat& frame)
	{
		double param[6];
		Eigen2Ceres(rigid_rotate_, rigid_translate_, param);

		// calc depth map
		depth_map = FilterDepth(frame);
		cv::Mat points = ComputePoints(depth_map, camera_.fx, camera_.fy, camera_.cx, camera_.cy);
		normal_map = ComputeNormal(points);

		InitializePhase1(frame, normal_map, param);
		Ceres2Eigen(rigid_rotate_, rigid_translate_, param);////
		for (int residual_count = 0; residual_count < vertex_size; residual_count++)
		{
			auto error = ProjectICPErrorPoseAndCoeff(residual_count, frame, normal_map, pca_models_, 
				camera_.fx, camera_.fy, camera_.cx, camera_.cy, 
				face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height, false);
			double* residual = new double[4];
			error(param, param + 3, pca_coeff_.data(), residual);
			//model_attribute_(0, residual_count) = residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3];
			mesh_.m_Colors[residual_count][0] = std::min(1.0, 4000 * (residual[0] * residual[0] + residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3]));
			delete[] residual;
		}

		//InitializePhase2(frame, normal_map, param);
		//for (int residual_count = 0; residual_count < vertex_size; residual_count++)
		//{
		//	auto error = ProjectICPErrorCoeff(residual_count, frame, normal, pca_models_, camera_.fx, camera_.fy, camera_.cx, camera_.cy, face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height, param);
		//	double* residual = new double[4];
		//	error(pca_coeff_.data(), residual);
		//	//model_attribute_(0, residual_count) = residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3];
		//	mesh_.m_Colors[residual_count][0] = std::min(1.0, 2000 * (residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3]));
		//	delete[] residual;
		//}
		
		for (int i = 0; i < pca_size; i++)
			cout << pca_coeff_[i] << " ";
		cout << rigid_rotate_ << "\n" << rigid_translate_ << "\n";
		std::cout << "\n** frame " << frame_count << " done **\n\n";

		WriteNeturalFace();
	}

	void Track(cv::Mat& frame)
	{
		UpdateNeturalFace();
		UpdateBlendshape();

		// calc depth map
		depth_map = FilterDepth(frame);
		cv::Mat points = ComputePoints(depth_map, camera_.fx, camera_.fy, camera_.cx, camera_.cy);
		normal_map = ComputeNormal(points);

		//WriteDepthImage(frame, iframe);
		//return;

		double param[6];
		Eigen2Ceres(rigid_rotate_, rigid_translate_, param);

		ceres::Problem problem;
		ceres::Solver::Options options;
		options.linear_solver_type = ceres::SPARSE_SCHUR;
		options.minimizer_progress_to_stdout = false;
		options.max_num_iterations = 20;
		options.num_threads = 1;
		ceres::LossFunctionWrapper* loss_function_wrapper = new ceres::LossFunctionWrapper(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
		for (int residual_count = 0; residual_count < vertex_size; residual_count++)
		{
			if (!mask_reader.face_mask[residual_count])
				continue;

			problem.AddResidualBlock(
				TrackFitError::Create(residual_count, frame, normal_map, exp_models_, camera_.fx, camera_.fy, camera_.cx, camera_.cy, face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height),
				loss_function_wrapper, param, param + 3, exp_coeff_.data()
				);
		}

		problem.AddResidualBlock(
			TrackSmoothError::Create(prev_exp_coeff_, prev_prev_exp_coeff_),
			0, exp_coeff_.data()
			);

		problem.AddResidualBlock(
			TrackSparseError::Create(),
			0, exp_coeff_.data()
			);

		ceres::Solver::Summary summary;
		ceres::Solve(options, &problem, &summary);

		prev_prev_exp_coeff_ = prev_exp_coeff_;
		prev_exp_coeff_ = exp_coeff_;

		//for (int i = 0; i < exp_size; i++)
		//	cout << exp_coeff_[i] << " ";
		std::cout << "\n** frame "<< frame_count <<" done **\n\n";

		//for (int residual_count = 0; residual_count < vertex_size; residual_count++)
		//{
		//	auto error = TrackFitError(residual_count, frame, normal_map, exp_models_, camera_.fx, camera_.fy, camera_.cx, camera_.cy, face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height);
		//	double* residual = new double[4];
		//	error(param, param + 3, exp_coeff_.data(), residual);
		//	//model_attribute_(0, residual_count) = residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3];
		//	mesh_.m_Colors[residual_count][0] = std::min(1.0, 2000 * (residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3]));
		//	delete[] residual;
		//}

		UpdateExpressionFace();
		if (frame_count > 0){
			WriteExpressionFace();
		}
	}

	void Refine()
	{
		UpdateExpressionFace();
		//if (frame_count > 20){
			//WriteExpressionFace();
		//}
		//return;

		SpMat A(4 * vertex_size, 3 * vertex_size);
		std::vector<Tripletd> tris;
		tris.reserve(12 * vertex_size);
		MatrixXd C(4 * vertex_size, 1);

		for (int i = 0; i < vertex_size; i++){
			if (!mask_reader.face_mask[i]){
				for (int j = 0; j < 4; j++){
					for (int k = 0; k < 3; k++)
						tris.push_back(Tripletd(4 * i + j, 3 * i + k, 0));
					C(4 * i + j, 0) = 0;
				}
				continue;
			}
			Vector3d p, pp;
			p(0) = model_(3 * i + 0);
			p(1) = model_(3 * i + 1);
			p(2) = model_(3 * i + 2);
			//p = rigid_rotate_ * p + rigid_translate_;
			pp(0) = p(0) / p(2) * (-1) * camera_.fx + camera_.cx;
			pp(1) = p(1) / p(2) * (-1) * camera_.fy + camera_.cy;
			double wx = pp(0) - (int)pp(0);
			double wy = pp(1) - (int)pp(1);
			// depth
			double zz1 = depth_map.at<unsigned short>(pp(1), pp(0)) / 1000.0;
			double zz2 = depth_map.at<unsigned short>(pp(1) + 1, pp(0)) / 1000.0;
			double zz3 = depth_map.at<unsigned short>(pp(1), pp(0) + 1) / 1000.0;
			double zz4 = depth_map.at<unsigned short>(pp(1) + 1, pp(0) + 1) / 1000.0;
			pp(2) = (1 - wx) * (1 - wy) * zz1
				+ (1 - wx) * wy * zz2
				+ wx * (1 - wy) * zz3
				+ wx * wy * zz4;
			// normal
			cv::Vec3f n1 = normal_map.at<cv::Vec3f>(pp(1), pp(0));
			cv::Vec3f n2 = normal_map.at<cv::Vec3f>(pp(1) + 1, pp(0));
			cv::Vec3f n3 = normal_map.at<cv::Vec3f>(pp(1), pp(0) + 1);
			cv::Vec3f n4 = normal_map.at<cv::Vec3f>(pp(1) + 1, pp(0) + 1);
			cv::Vec3f n = (1 - wx) * (1 - wy) * n1
				+ (1 - wx) * wy * n2
				+ wx * (1 - wy) * n3
				+ wx * wy * n4;
			Matrix<double, 4, 3> mix;
			mix.setZero();
			mix(0, 0) = mix(1, 1) = mix(2, 2) = 1;
			mix(3, 0) = n[0] * 1.11;
			mix(3, 1) = n[1] * 1.11;
			mix(3, 2) = n[2] * 1.11;
			//Matrix<double, 1, 3> nn(n[0], n[1], n[2]);

			// get closest point to (x, y, z) as p(xx, yy, zz)
			pp(0) = -1 * (pp(0) - camera_.cx) * pp(2) / camera_.fx;
			pp(1) = -1 * (pp(1) - camera_.cy) * pp(2) / camera_.fy;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
			//
			Matrix<double, 4, 3> lhs = mix * rigid_rotate_;
			Matrix<double, 4, 1> rhs = mix * (pp - rigid_translate_);
			//std::cout << lhs << "\n" << rhs << "\n";
			for (int j = 0; j < 4; j++){
				for (int k = 0; k < 3; k++){
					tris.push_back(Tripletd(4 * i + j, 3 * i + k, lhs(j, k)));
				}
				C(4 * i + j, 0) =  rhs(j, 0);
			}
			/*auto tmp = lhs * p - rhs;
			std::cout << p(0) << " " << p(1) << " " << p(2) << "\n";
			std::cout << pp(0) << " " << pp(1) << " " << pp(2) << "\n";
			std::cout << tmp(0) << " " << tmp(1) << " " << tmp(2) << " " << tmp(3) << "\n";*/
		}
		A.setFromTriplets(tris.begin(), tris.end());

		SpMat sm(3 * vertex_size, 3 * vertex_size);
		sm.setIdentity();
		sm = sm * exp_coeff_[0];
		for (int i = 1; i < exp_size; i++){
			sm += exp_coeff_[i] * exp_transforms_[i];
		}
		MatrixXd A_ = A * sm * pca_models_.rightCols(pca_size - 1);
		MatrixXd C_ = C - A * sm * pca_models_.block(0, 0, vertex_size * 3, 1);

		double gamma = 0.9;

		double St_new = gamma * S_Refine + 1;
		M_Refine = gamma * S_Refine / St_new * M_Refine + 1 / St_new * A_.transpose() * A_;
		Y_Refine = gamma * S_Refine / St_new * Y_Refine + 1 / St_new * A_.transpose() * C_;
		S_Refine = St_new;
		//std::cout << Mt.rows() << " " << Mt.cols() << "\n";
		//std::cout << A_.rows() << " " << A_.cols() << "\n";

		ConjugateGradient<MatrixXd, Lower | Upper> cg;
		MatrixXd tmp = M_Refine + 1.0 *D_Refine;
		cg.compute(tmp);
		auto result = cg.solveWithGuess(Y_Refine, pca_coeff_.segment(1, pca_size - 1));

		for (int i = 0; i < pca_size - 1; i++)
			std::cout << result(i) << " ";
		std::cout << "\n";
		std::cout << "#iterations:     " << cg.iterations() << std::endl;
		std::cout << "estimated error: " << cg.error() << std::endl;
		//system("pause");

		//saveMarket(A, Desktop_Path + "A");
		//saveMarket(A_, Desktop_Path + "A_");
		//saveMarket(C, Desktop_Path + "C");
		//saveMarket(C_, Desktop_Path + "C_");
		//saveMarket(sm, Desktop_Path + "sm");

		for (int i = 1; i < pca_size; i++)
			pca_coeff_(i) = result(i - 1);
	}

	void UpdateNeturalFace()
	{
		exp_models_.block(0, 0, vertex_size * 3, 1) = pca_models_ * pca_coeff_;
	}

	void UpdateBlendshape()
	{
		//for (int i = 0; i < vertex_size; i++){
		//	exp_models_[0][3 * i + 2] *= -1;
		//}
		for (int i = 1; i < exp_size; i++){
			exp_models_.block(0, i, vertex_size * 3, 1) = exp_transforms_[i] * exp_models_.block(0, 0, vertex_size * 3, 1);
		}
		//for (int j = 0; j < exp_size; j++){
		//	for (int i = 0; i < vertex_size; i++){
		//		exp_models_[j][3 * i + 2] *= -1;
		//	}
		//}
	}

	void UpdateExpressionFace()
	 {
		 exp_coeff_(0) = 1;
		 for (int i = 1; i < exp_size; i++)
			 exp_coeff_(0) -= exp_coeff_(i);
		 model_ = exp_models_ * exp_coeff_;

		 for (int i = 0; i < vertex_size; i++){
			 Vector3d v(model_(3 * i + 0), model_(3 * i + 1), model_(3 * i + 2));
			 v = rigid_rotate_ * v + rigid_translate_;
			 model_(3 * i + 0) = v(0);
			 model_(3 * i + 1) = v(1);
			 model_(3 * i + 2) = v(2);
		 }
	 }

	void WriteExpressionFace()
	{
		UpdateMeshVertex(model_, mesh_);
		char str[100];  
		sprintf(str, "e%d.obj", frame_count);
		ml::MeshIOd::saveToOBJ(Test_Output_Dir + str, mesh_);
	}

	void WriteNeturalFace()
	{
		UpdateNeturalFace();
		UpdateMeshVertex(exp_models_.col(0), mesh_);
		char str[100];
		sprintf(str, "netrural%d.obj", frame_count);
		ml::MeshIOd::saveToOBJ(Desktop_Path + str, mesh_);
	}

	void WriteDepthImage(cv::Mat& dm, cv::Mat& im)
	{
		ml::MeshDatad tm;
		tm.m_Vertices.resize(dm.cols * dm.rows);
		tm.m_Colors.resize(dm.cols * dm.rows);
		for (int i = 0; i < dm.cols; i++){
			for (int j = 0; j < dm.rows; j++){
				Vector3d tmp = Point2d_2_Point3d(Vector2d(i, j), dm.at<ushort>(j, i));
				tm.m_Vertices[i * dm.rows + j] = ml::vec3d(tmp(0), tmp(1), tmp(2));
				double color = std::min(((int)im.at<uchar>(j, i)) * 3 / 256.0, 1.0);
				tm.m_Colors[i * dm.rows + j] = ml::vec4d(color, color, color, 1.0);
			}
		}
		char str[200];
		sprintf(str, "dem_pc%d.obj", frame_count);
		ml::MeshIOd::saveToOBJ(Test_Output_Dir + str, tm);
	}

	void WriteDepthImage(cv::Mat& dm)
	{
		depth_map = FilterDepth(dm);
		cv::Mat points = ComputePoints(depth_map, camera_.fx, camera_.fy, camera_.cx, camera_.cy);
		normal_map = ComputeNormal(points);

		ml::MeshDatad tm;
		tm.m_Vertices.resize(dm.cols * dm.rows);
		//tm.m_Normals.resize(dm.cols * dm.rows);
		tm.m_Colors.resize(dm.cols * dm.rows);
		for (int i = 0; i < dm.cols; i++){
			for (int j = 0; j < dm.rows; j++){
				Vector3d tmp = Point2d_2_Point3d(Vector2d(i, j), dm.at<ushort>(j, i));
				tm.m_Vertices[i * dm.rows + j] = ml::vec3d(tmp(0), tmp(1), tmp(2));
				cv::Vec3f n = normal_map.at<cv::Vec3f>(j, i);
				//if (n[0] > 0.5)
				tm.m_Colors[i * dm.rows + j] = ml::vec4d((n[2] + 1) / 2, 0, 0, 1);
			}
		}
		char str[200];
		sprintf(str, "dem_pc%d.obj", frame_count);
		ml::MeshIOd::saveToOBJ(Test_Output_Dir + str, tm);
	}

public:
	FaceDetector fd;
	Camera camera_;
	cv::Rect face_rect_;
	MaskReader mask_reader;
	int frame_count;
	cv::Mat normal_map;
	cv::Mat depth_map;

public:
	// pca basic info
	double xw = 0.180881292;
	double yw = 0.200453609;
	double zw = 0.154334202;
	Vector3d pca_center = Vector3d(-0.000269647688, -0.00578920171, -0.0547919050);

	//ml::mat3d rigid_rotate_;
	//ml::vec3d rigid_tranlate_;
	Matrix3d rigid_rotate_;
	Vector3d rigid_translate_;

	VectorXd pca_weights_;
	SpMat D_Refine;
	MatrixXd M_Refine;
	VectorXd Y_Refine;
	double S_Refine;
	/*MatrixXd lap_weights_;*/

	ml::MeshDatad mesh_;
	VectorXd model_;
	MatrixXd pca_models_;
	MatrixXd lap_models_;
	MatrixXd exp_models_;
	std::vector<SpMat> exp_transforms_;

	VectorXd exp_coeff_, prev_exp_coeff_, prev_prev_exp_coeff_;//x
	VectorXd pca_coeff_;//y
	MatrixXd lap_coeff_;//z
};