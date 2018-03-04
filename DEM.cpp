#include "DEM.h"

AdaptiveDynamicExpressionModel::AdaptiveDynamicExpressionModel()
{
	//mask_reader.Test();

	rigid_rotate_.setIdentity();
	rigid_tranlate_.setZero();

	pca_coeff_.resize(pca_size, 0);
	pca_coeff_[0] = 1;
	lap_coeff_.resize(pca_size, std::vector<double>(lap_size, 0));

	model_.resize(Eigen::NoChange, vertex_size);
	model_attribute_.resize(Eigen::NoChange, vertex_size);
	model_attribute_.setZero();
	ModelReader mreader(pca_models_, lap_models_, exp_models_, pca_weights_);
	/*float xmin, ymin, zmin;
	float xmax, ymax, zmax;
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
	float xw = xmax - xmin;
	float xo = (xmin + xmax) / 2;
	float yw = ymax - ymin;
	float yo = (ymin + ymax) / 2;
	float zw = zmax - zmin;
	float zo = (zmin + zmax) / 2;
	std::cout << xw << xo << yw << yo << zw << zo;*/
}

Eigen::Vector3f AdaptiveDynamicExpressionModel::Point2d_2_Point3f(Eigen::Vector2f p2, int depth)
{
	if (depth == 0)
		return Eigen::Vector3f(0, 0, 0);
	float x_z = (p2(0) - camera_.cx) / camera_.fx;
	float y_z = (p2(1) - camera_.cy) / camera_.fy;
	Eigen::Vector3f p3;
	p3(2) = depth * 0.001 * sqrt(1 / (x_z * x_z + y_z * y_z + 1));
	p3(0) = -x_z * p3(2);
	p3(1) = -y_z * p3(2);
	return p3;
}

cv::Point2f AdaptiveDynamicExpressionModel::Point3f_2_Point2d(cv::Point3f p3)
{
	if (p3.z == 0)
		return cv::Point2f(0, 0);
	cv::Point2f p2;
	p2.x = p3.x / p3.z * (-1) * camera_.fx + camera_.cx;
	p2.y = p3.y / p3.z * (-1) * camera_.fy + camera_.cy;

	return p2;
}

bool AdaptiveDynamicExpressionModel::UpdateFaceRegion(cv::Mat& gray_frame)
{
	if (!fd.Detect(gray_frame, face_rect_))
		return false;
	return true;
}

bool AdaptiveDynamicExpressionModel::InitRigidMotion(cv::Mat& gray_frame, cv::Mat& depth_frame)
{
	if (!UpdateFaceRegion(gray_frame))
		return false;
	//rectangle(depth_frame, cv::Point(face_rect_.x, face_rect_.y), cv::Point(face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height), cv::Scalar(0));
	//showImage(depth_frame);

	Eigen::Vector2f p2(face_rect_.x + face_rect_.width / 2, face_rect_.y + face_rect_.height / 2);
	float depth = depth_frame.at<ushort>(p2(1), p2(0));
	Eigen::Vector3f p3 = Point2d_2_Point3f(p2, depth);
	rigid_tranlate_ = p3 - pca_center;

	//int delta = 2;
	//face_rect_.x -= delta;
	//face_rect_.y -= delta;
	//face_rect_.width += delta * 2;
	//face_rect_.height += delta * 2;

	return true;
}

void AdaptiveDynamicExpressionModel::InitializePhase1(cv::Mat& frame, cv::Mat& normal, double* param)
{
	ceres::Problem problem;
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::SPARSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	options.max_num_iterations = 20;
	options.num_threads = 1;
	ceres::LossFunctionWrapper* loss_function_wrapper = new ceres::LossFunctionWrapper(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
	for (int residual_count = 0; residual_count < vertex_size; residual_count++)
	{
		if (!mask_reader.face_mask[residual_count * vertex_step])
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

void AdaptiveDynamicExpressionModel::InitializePhase2(cv::Mat& frame, cv::Mat& normal, double* param)
{
	ceres::Problem problem;
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::SPARSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	options.max_num_iterations = 20;
	options.num_threads = 1;
	ceres::LossFunctionWrapper* loss_function_wrapper = new ceres::LossFunctionWrapper(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
	for (int residual_count = 0; residual_count < vertex_size; residual_count++)
	{
		if (!mask_reader.face_mask[residual_count * vertex_step])
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

void AdaptiveDynamicExpressionModel::TwoPhaseInitialize(cv::Mat& frame)
{
	double param[6];
	Eigen2Ceres(rigid_rotate_, rigid_tranlate_, param);

	// calc depth map
	cv::Mat filter_depth = FilterDepth(frame);
	cv::Mat points = ComputePoints(filter_depth, camera_.fx, camera_.fy, camera_.cx, camera_.cy);
	cv::Mat normal = ComputeNormal(points);

	InitializePhase1(frame, normal, param);
	Ceres2Eigen(rigid_rotate_, rigid_tranlate_, param);////
	model_attribute_.setZero();
	for (int residual_count = 0; residual_count < vertex_size; residual_count++)
	{
		auto error = ProjectICPErrorPoseAndCoeff(residual_count, frame, normal, pca_models_, camera_.fx, camera_.fy, camera_.cx, camera_.cy, face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height);
		double* residual = new double[4];
		error(param, param + 3, pca_coeff_.data(), residual);
		model_attribute_(0, residual_count) = residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3];

		delete[] residual;
	}
	WritePointCloud("C:\\Users\\zhx\\Desktop\\phase1.ply");

	InitializePhase2(frame, normal, param);
	model_attribute_.setZero();
	for (int residual_count = 0; residual_count < vertex_size; residual_count++)
	{
		auto error = ProjectICPErrorCoeff(residual_count, frame, normal, pca_models_, camera_.fx, camera_.fy, camera_.cx, camera_.cy, face_rect_.x, face_rect_.y, face_rect_.x + face_rect_.width, face_rect_.y + face_rect_.height, param);
		double* residual = new double[4];
		error(pca_coeff_.data(), residual);
		model_attribute_(0, residual_count) = residual[1] * residual[1] + residual[2] * residual[2] + residual[3] * residual[3];

		delete[] residual;
	}
	WritePointCloud("C:\\Users\\zhx\\Desktop\\phase2.ply");
}

void AdaptiveDynamicExpressionModel::WritePointCloud(std::string path)
{
	Vertices output_model = pca_models_[0] * pca_coeff_[0];
	for (int i = 1; i < pca_size; i++)
	{
		output_model += pca_models_[i] * pca_coeff_[i];
	}
	PrintPointCloud(output_model, model_attribute_, rigid_rotate_, rigid_tranlate_, path);
}

void AdaptiveDynamicExpressionModel::WriteMesh()
{
	std::string head;
	std::string end;
	std::fstream ifs;
	ifs.open("F:\\FaceAnimation\\pca\\0.ply");
	std::string line;
	for (int i = 0; i < 14; i++)
	{
		std::getline(ifs, line);
		head += line + "\n";
	}
	for (int i = 0; i < vertex_total_size; i++)
		std::getline(ifs, line);
	while (std::getline(ifs, line))
	{
		end += line + "\n";
	}
	ifs.close();

	Vertices output_model = pca_models_[0] * pca_coeff_[0];
	for (int i = 1; i < pca_size; i++)
	{
		output_model += pca_models_[i] * pca_coeff_[i];
	}
	output_model = rigid_rotate_ * output_model;
	std::ofstream ofs;
	ofs.open("C:\\Users\\zhx\\Desktop\\mesh.ply");
	ofs << head;
	for (int i = 0; i < vertex_size; i++)
	{
		ofs << output_model(0, i) + rigid_tranlate_(0) << " "
			<< output_model(1, i) + rigid_tranlate_(1) << " "
			<< output_model(2, i) + rigid_tranlate_(2) << " 220 220 220 255\n";
	}
	ofs << end;
	ofs.close();
}

void AdaptiveDynamicExpressionModel::WriteInitResult()
{
	std::ofstream ofs;
	ofs.open("C:\\Users\\zhx\\Desktop\\init.txt");
	ofs << rigid_rotate_ << "\n" << rigid_tranlate_;
	ofs << "\n";
	for (int i = 0; i < pca_coeff_.size(); i++)
		ofs << pca_coeff_[i] << " ";
	ofs.close();
}

void AdaptiveDynamicExpressionModel::ReadInitResult()
{
	std::ifstream ifs;
	ifs.open("C:\\Users\\zhx\\Desktop\\init.txt");
	ifs >> rigid_rotate_(0, 0) >> rigid_rotate_(0, 1) >> rigid_rotate_(0, 2) >>
		rigid_rotate_(1, 0) >> rigid_rotate_(1, 1) >> rigid_rotate_(1, 2) >>
		rigid_rotate_(2, 0) >> rigid_rotate_(2, 1) >> rigid_rotate_(2, 2);
	ifs >> rigid_tranlate_(0) >> rigid_tranlate_(1) >> rigid_tranlate_(2);
	for (int i = 0; i < pca_coeff_.size(); i++)
		ifs >> pca_coeff_[i];
	ifs.close();
}