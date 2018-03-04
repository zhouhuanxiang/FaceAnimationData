#pragma once

#include <assimp/Importer.hpp>      // 导入器在该头文件中定义
#include <assimp/scene.h>           // 读取到的模型数据都放在scene中
#include <assimp/postprocess.h>     // 该头文件中包含后处理的标志位定义

#include "defines.h"
#include "file_func.h"
#include "./SICP/nanoflann.hpp"
#include "./mlibutil/mlibutil.h"
#include "eigen_binary_io.h"

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace nanoflann;

template <typename T>
struct MeshPointCloud
{
	MeshPointCloud()
	{}

	MeshPointCloud(ml::MeshData<T>& mesh)
		:mesh(mesh)
	{}

	ml::MeshData<T>& mesh;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return mesh.m_Vertices.size(); }

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		return mesh.m_Vertices[idx][dim];
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

};

class ModelPreprocess
{
public:
	//ModelPreprocess()
	//{}

	ModelPreprocess()
	{}

	template <class T>
	void WriteMesh(ml::MeshData<T>& mesh)
	{
		ml::MeshIO<T>::saveToOBJ(Test_Output_Dir + "sha.obj", mesh);
	}

	void RewriteMaskFromDae()
	{
		// load mesh1
		ml::MeshDatad mesh1;
		ml::MeshIOd::loadFromOBJ(Test_Output_Dir + "Normal_tri.obj", mesh1);

		// load mesh2
		ml::MeshDatad mesh2;
		Assimp::Importer import;
		const aiScene *scene = import.ReadFile(Test_Output_Dir + "sha.dae", aiProcess_Triangulate | aiProcess_FlipUVs);
		if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
		{
			cout << "ERROR::ASSIMP::" << import.GetErrorString() << endl;
		}
		aiMesh *mesh = scene->mMeshes[0];
		mesh2.m_Colors.resize(mesh->mNumVertices);
		mesh2.m_Vertices.resize(mesh->mNumVertices);
		mesh2.m_FaceIndicesVertices.clear();
		for (unsigned int i = 0; i < mesh->mNumVertices; i++){
			mesh2.m_Vertices[i][0] = mesh->mVertices[i].x;
			mesh2.m_Vertices[i][1] = mesh->mVertices[i].y;
			mesh2.m_Vertices[i][2] = mesh->mVertices[i].z;

			mesh2.m_Colors[i][0] = mesh->mColors[0][i].r;
			mesh2.m_Colors[i][1] = mesh->mColors[0][i].g;
			mesh2.m_Colors[i][2] = mesh->mColors[0][i].b;
		}

		// build KD-tree
		MeshPointCloud<double> cloud(mesh2);
		typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, MeshPointCloud<double> >, MeshPointCloud<double>, 3 /* dim */> my_kd_tree_t;
		my_kd_tree_t index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
		index.buildIndex();
		size_t num_results = 3;
		std::vector<size_t> ret_index(num_results);
		std::vector<double> out_dist_sqr(num_results);

		// query KD-tree and acquire mask
		int vnum = mesh1.m_Vertices.size();
		int fnum = mesh1.m_FaceIndicesVertices.size();
		std::vector<bool> vmask1(vnum, false);
		std::vector<bool> vmask2(vnum, false);
		for (int j = 0; j < mesh1.m_Vertices.size(); j++){
			index.knnSearch((double*)&mesh1.m_Vertices[j], num_results, &ret_index[0], &out_dist_sqr[0]);
			// In case of less points in the tree than requested:
			ret_index.resize(num_results);
			out_dist_sqr.resize(num_results);

			/*for (size_t i = 0; i < num_results; i++)
				cout << out_dist_sqr[i] << " ";
			cout << "\n" << mesh1.m_Vertices[j][0] << " " << mesh1.m_Vertices[j][1] << " " << mesh1.m_Vertices[j][2] << "\n" <<
				mesh2.m_Vertices[ret_index[0]][0] << " " << mesh2.m_Vertices[ret_index[0]][1] << " " << mesh2.m_Vertices[ret_index[0]][2] << "\n#####\n";*/
			if (mesh2.m_Colors[ret_index[0]][1] == 0)
				vmask1[j] = true;
		}
		vmask2 = vmask1;
		
		// spread mask
		for (int iter = 0; iter < 3; iter++){
			vmask2 = vmask1;
			for (int i = 0; i < fnum; i++){
				int tsize = mesh1.m_FaceIndicesVertices[i].size();
				for (int j = 0; j < tsize; j++){
					if (vmask1[mesh1.m_FaceIndicesVertices[i][j]]){
						for (int k = 0; k < tsize; k++){
							vmask2[mesh1.m_FaceIndicesVertices[i][k]] = true;
						}
						break;
					}
				}
			}
			vmask1 = vmask2;
		}

		// write color as mask
		mesh1.m_Colors.clear();
		mesh1.m_Colors.resize(vnum, ml::vec4f(1.0, 1.0, 1.0, 1.0));
		for (int i = 0; i < vnum; i++){
			if (vmask2[i])
				mesh1.m_Colors[i] = ml::vec4f(1.0, 0.0, 0.0, 0.0);
			// optional
			//ml::vec3f v = mesh1.m_Vertices[i];
			//mesh1.m_Vertices[i][1] = v[2];
			//mesh1.m_Vertices[i][2] = v[1] * -1;
			//mesh1.m_Vertices[i] = mesh1.m_Vertices[i] * 1000.0;
		}
		WriteMesh(mesh1);
	}

	// only use first 50 pca face
	void AdjustSimplifiedPcaModel()
	{
		Eigen::MatrixXd pcas(3 * vertex_size, pca_size);
		for (int i = 0; i <= 50; i++)
		{
			ml::MeshDatad mesh;
			char str[200];
			sprintf(str, "pca_obj_simplified1/pca/%d.obj", i);
			ml::MeshIOd::loadFromOBJ(Data_Input_Dir + str, mesh);

			for (int j = 0; j < mesh.m_Vertices.size(); j++){
				pcas(3 * j + 0, i) = mesh.m_Vertices[j].x * -0.001;
				pcas(3 * j + 1, i) = mesh.m_Vertices[j].y * 0.001;
				pcas(3 * j + 2, i) = mesh.m_Vertices[j].z * -0.001;
				mesh.m_Vertices[j].x *= -0.001;
				mesh.m_Vertices[j].y *= 0.001;
				mesh.m_Vertices[j].z *= -0.001;
			}

			ml::MeshIOd::saveToOBJ(Data_Input_Dir + str, mesh);
		}
		write_binary<Eigen::MatrixXd>(Data_Input_Dir + "pca_obj_simplified1/pca/pca", pcas);
	}

	// only use first 50 pca face
	void CalcPcaWeight()
	{
		ml::MeshDataf mesh;
		ml::MeshIOf::loadFromPLY(Data_Input_Dir + "pca_ply/0.ply", mesh);
		std::vector<double> w;
		for (int i = 0; i < 50; i++){
			ml::MeshDataf tm;
			char str[200];
			sprintf(str, "pca_ply/%d.ply", i + 1);
			ml::MeshIOf::loadFromPLY(Data_Input_Dir + str, tm);
			float error = 0;
			for (int j = 0; j < mesh.m_Vertices.size(); j++){
				for (int k = 0; k < 3; k++){
					error += pow(tm.m_Vertices[j][k] - mesh.m_Vertices[j][k], 2);
				}
			}
			std::cout << error << "\n";
			w.push_back(std::sqrt(error));
		}
		std::ofstream ofs;
		ofs.open(Test_Output_Dir + "pca_weight.txt");
		ofs << "0.0 ";
		for (int i = 0; i < 50; i++)
		{
			ofs << w[i] / w[0] << " ";
		}
		ofs.close();
	}

private:

};