#include "image_reader_kinect.h"

#include "expression_generator.h"
#include "pca_generator_3d.h"
#include "expression_transfer.h"
#include "model_preprocess.h"
#include "eigen_binary_io.h"

int main()
{
	//ml::MeshDatad mesh1, mesh2;
	//ml::MeshIOd::loadFromOBJ("F:/FaceAnimation/pca_obj/0.obj", mesh1);
	//ml::MeshIOd::loadFromOBJ("F:/FaceAnimation/pca_obj/tmp/00sim.obj", mesh2);
	//ml::BoundingBox3d bb1 = mesh1.computeBoundingBox();
	//ml::BoundingBox3d bb2 = mesh2.computeBoundingBox();
	//ml::vec3d v =  bb1.getCenter() - bb2.getCenter() / 10.0;
	////std::cout << bb1.getExtentX() << " " << bb2.getExtentX() << " " << bb1.getExtentX() / bb2.getExtentY();
	////std::cout << bb1.getExtentY() << " " << bb2.getExtentY() << " " << bb1.getExtentY() / bb2.getExtentY();
	////std::cout << bb1.getExtentZ() << " " << bb2.getExtentZ() << " " << bb1.getExtentZ() / bb2.getExtentZ();
	//for (int i = 0; i < mesh2.m_Vertices.size(); i++)
	//{
	//	mesh2.m_Vertices[i] = mesh2.m_Vertices[i] / 10.0 + v;
	//	mesh2.m_Vertices[i] = mesh2.m_Vertices[i] * 1000;
	//	mesh2.m_Vertices[i].z *= -1;
	//}
	//ml::MeshIOd::saveToOBJ("F:/FaceAnimation/pca_obj/02sim.obj", mesh2);

	//ml::MeshDatad mesh;
	//ml::MeshIOd::loadFromOBJ(Desktop_Path + "ss.obj", mesh);
	//mesh.makeTriMesh();
	//ml::MeshIOd::saveToOBJ(Desktop_Path + "2.obj", mesh);
	//return 0;


	//Eigen::VectorXd v(9);
	//v << 1, 1, 1, 2, 2, 2, 3, 3, 3;
	//Eigen::Vector3d t;
	//t << 1, 1, 1;
	//Eigen::Matrix3d m;
	//m << 1, 0, 0, 0, 1, 0, 0, 0, 1;

	//Eigen::Ref<Eigen::MatrixXd> block = m.block(1, 1, 1, 2);
	//block(0, 0) = -1;
	//std::cout << m;

	//system("pause");
	//return 0;


#if 0
	ModelPreprocess mp;
	//mp.RewriteMaskFromDae();
	mp.AdjustSimplifiedPcaModel();
	return 0;
#endif

#if 0
	PcaGenerator3D pg;
	pg.Run();
	return 0;
#endif

#if 1
	ExpressionTransfer et(Data_Input_Dir + "pca_obj_simplified1/original.obj", Data_Input_Dir + "pca_obj_simplified1/expression");
	et.Solve();
	return 0;
	for (et.e_index = 0; et.e_index <= 28; et.e_index++){
		et.TestCompress();
		break;
	}
	system("pause");
	return 0;
#endif

#if 0
	ExpressionGenerator eg;
	eg.Run();
	return 0;
#endif

#if 0
	//SpMat v, t;
	//std::vector<T> tripletList;
	//t.resize(3 * vertex_size, 3 * vertex_size);
	//tripletList.resize(3 * vertex_size);
	//for (int i = 0; i < vertex_size; i++){
	//	tripletList.push_back(T(3 * i + 0, 3 * i + 0, 1));
	//	tripletList.push_back(T(3 * i + 1, 3 * i + 1, 1));
	//	tripletList.push_back(T(3 * i + 2, 3 * i + 2, -1));
	//}
	//t.setFromTriplets(tripletList.begin(), tripletList.end());
	//for (int i = 0; i < exp_size - 1; i++){
	//	char str[100];
	//	sprintf(str, "transformation%d_9", i);
	//	Deserialize<double>(v, Data_Input_Dir + "exp_transform/" + str);
	//	if (v.cols() != vertex_size * 3){
	//		printf("bad exp transform %d\n", i);
	//		system("pause");
	//	}
	//	v = v * t;
	//	Serialize<double>(v, Data_Input_Dir + "exp_transform/" + str);
	//}

	//cv::Mat infrared_frame;
	cv::Mat depth_frame;
	ImageReaderKinect image_reader("D:/project2013/render_kinect-master/result/");
	//
	//image_reader.GetFrame(0, infrared_frame, depth_frame);
	//FaceDetector fd;
	//cv::Rect rect;
	//fd.Detect(infrared_frame, rect);
	//return 0;

	AdaptiveDynamicExpressionModel dem;
	//
	//dem.ReadInitResult();
	//dem.WriteMesh();
	//return 0;

	int frame_count = 0;
	image_reader.GetFrame(frame_count, depth_frame);
	// infrared 1
	/*if (!dem.InitRigidMotion(infrared_frame, depth_frame)){
		std::cout << "no face detected";
	}*/

	// infrared 1
	//dem.face_rect_ = cv::Rect(247, 100, 393 - 247, 265 - 100);
	dem.face_rect_ = cv::Rect(cv::Point2d(258, 144), cv::Point2d(382, 325));
	dem.rigid_translate_ = Eigen::Vector3d(0, 0, 0.7);

	for (;frame_count <= 3;){
		frame_count += 1;
		// infrared 3
		/*image_reader.GetFrame(frame_count, infrared_frame, depth_frame);
		dem.UpdateFaceRegion(infrared_frame);*/
		image_reader.GetFrame(frame_count, depth_frame);

		dem.frame_count = frame_count;
		//dem.WriteDepthImage(depth_frame);
		//continue;
		dem.TwoPhaseInitialize(depth_frame);
	}
	system("pause");
	return 0;

	for (; frame_count < 150;)
	{
		frame_count += 3;
		// infrared 4
		/*image_reader.GetFrame(frame_count, infrared_frame, depth_frame);
		dem.UpdateFaceRegion(infrared_frame);*/
		image_reader.GetFrame(frame_count, depth_frame);

		dem.frame_count = frame_count;
		dem.Track(depth_frame);
		dem.Refine();
	}
	system("pause");
	return 0;
#endif
}