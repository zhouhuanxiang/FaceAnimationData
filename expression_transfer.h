#pragma once

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>
#include <unsupported/Eigen/SparseExtra>

//#include "./mlibutil/mlibutil.h"
#include "./file_func.h"
#include "defines.h"
#include "deformation_transfer.h"
#include "mask_reader.h"
#include "eigen_binary_io.h"

#include <map>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


class ExpressionTransfer
{
public:
	ExpressionTransfer() {}

	ExpressionTransfer(string nfilename, string dirname)
		:nfilename(nfilename)
	{
		GetAllFormatFiles(dirname, filenames, ".obj");
		ml::MeshIOd::loadFromOBJ(nfilename, nface);
		for (int i = 0; i < nface.m_Vertices.size(); i++){
			nface.m_Vertices[i].x *= -1;
			nface.m_Vertices[i].z *= -1;
		}
		transformation.resize(3 * nface.m_Vertices.size(), 3 * nface.m_Vertices.size());
		transformation.setZero();

		mask_reader.Read(nface.m_Vertices.size());
	}

	void ReadDeformation()
	{
		ifstream ifs(Test_Output_Dir + "deformation.txt");
		deformation.resize(nface.m_FaceIndicesVertices.size());
		for (int i = 0; i < deformation.size(); i++){
			for (int j = 0; j < 9; j++){
				ifs >> deformation[i][j];
			}
		}
		ifs.close();
		//for (int i = 0; i < 10; i++){
		//	for (int j = 0; j < 9; j++){
		//		cout << deformation[i][j] << " ";
		//	}
		//}
		//std::cout << "\n";
	}

	void CalcDeformation()
	{
		DeformationTransfer dt;
		dt.Solve(nface, eface);
		//dt.WriteDeformation();
		//dt.Test();

		deformation = dt.tar_deformation;
	}

	void Test()
	{
		int vnum = nface.m_Vertices.size();
		Eigen::VectorXd pts0(3 * vnum);
		for (int i = 0; i < vnum; i++)
		{
			pts0[i * 3 + 0] = nface.m_Vertices[i].x;
			pts0[i * 3 + 1] = nface.m_Vertices[i].y;
			pts0[i * 3 + 2] = nface.m_Vertices[i].z;
		}
		Eigen::VectorXd pts1 = transformation * pts0;

		ml::MeshDatad tmp;
		tmp.m_Vertices.resize(vnum);
		tmp.m_FaceIndicesVertices = nface.m_FaceIndicesVertices;
		double* ptr = &pts1[0];
		for (int i = 0; i < vnum; i++)
		{
			tmp.m_Vertices[i].x = *(double*)(ptr + i * 3 + 0);
			tmp.m_Vertices[i].y = *(double*)(ptr + i * 3 + 1);
			tmp.m_Vertices[i].z = *(double*)(ptr + i * 3 + 2);
		}
		char str[100];
		sprintf(str, "test_expression_transfer%d_%d.obj", e_index, compress_level);
		ml::MeshIOd::saveToOBJ(Test_Output_Dir + str, tmp);
	}

	void WriteTransformation()
	{
		//Eigen::SparseMatrix<float> tm = compressed_transformation.cast<float>();
		char str[200];
		if (binary_format)
		{
			sprintf(str, "transformation%d_%d", e_index, compress_level);
			Serialize<double>(transformation, Test_Output_Dir + str);
		}
		else
		{
			sprintf(str, "transformation%d_%d.mtx", e_index, compress_level);
			Eigen::saveMarket(transformation, Test_Output_Dir + str);
		}
	}

	void ReadTransformation(int e, int c)
	{
		char str[200];
		if (binary_format)
		{
			sprintf(str, "transformation%d_%d", e_index, compress_level);
			Deserialize<double>(transformation, Test_Output_Dir + str);
		}
		else
		{
			sprintf(str, "transformation%d_%d.mtx", e_index, compress_level);
			Eigen::loadMarket(transformation, Test_Output_Dir + str);
		}
	}

	//
	void TestCompress()
	{
		binary_format = true;
		//e_index = 28;
		compress_level = 11;
		ReadTransformation(e_index, compress_level);

		//
		/*ml::MeshDatad mesh1, mesh2;
		char str[100];
		sprintf(str, "F:/FaceAnimation/pca_obj_simplified2/original.obj");
		ml::MeshIOd::loadFromOBJ(str, mesh1);
		sprintf(str, "F:/FaceAnimation/pca_obj_simplified2/expression/%d.obj", e_index);
		ml::MeshIOd::loadFromOBJ(str, mesh2);
		std::vector<double> distance(vertex_size, 0);
		for (int i = 0; i < vertex_size; i++){
			
			ml::vec3d& v1 = mesh1.m_Vertices[i];
			ml::vec3d& v2 = mesh2.m_Vertices[i];
			double dist = (v1 - v2).lengthSq();
			distance[i] = dist;
		}*/
		//

		int vnum = nface.m_Vertices.size();
		SpMat prev = transformation;
		for (compress_level = 9; compress_level >= 9; compress_level--){
			std::vector<T> tripletList;
			double threshold = pow(0.5, compress_level);
			for (int i = 0; i < vnum; i++){
				if (!mask_reader.face_mask[i])
					continue;
				/*if (distance[i] < threshold){
					for (int k = 0; k < 3; k++){
						tripletList.push_back(T(3 * i + k, 3 * i + k, 1));
					}
					continue;
				}*/
				for (int j = 0; j < vnum; j++){
					ml::mat3d m;
					ml::vec3d v;
					for (int k = 0; k < 3; k++){
						for (int l = 0; l < 3; l++){
							m[3 * k + l] = prev.coeff(3 * i + k, 3 * j + l);
						}
						v[k] = nface.m_Vertices[j][k];
					}
					v = m * v;
					double weight = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
					if (weight > threshold){
						for (int k = 0; k < 3; k++)
							for (int l = 0; l < 3; l++)
								tripletList.push_back(T(3 * i + k, 3 * j + l, m[3 * k + l]));
					}
				}
			}

			transformation.setZero();
			transformation.setFromTriplets(tripletList.begin(), tripletList.end());
			Test();
			WriteTransformation();
		}
	}

	void Compress()
	{
		//
		//binary_format = true;
		//compress_level = 9;
		//for (e_index = 0; e_index <= 22; e_index++)
		//{
		//	ReadTransformation(e_index, compress_level);
		//	Test();
		//}
		//return;
		//
		int vnum = nface.m_Vertices.size();
		std::vector<T> tripletList;
		double threshold = pow(0.5, compress_level);
		for (int i = 0; i < vnum; i++){
			for (int j = 0; j < vnum; j++){
				ml::mat3d m;
				ml::vec3d v;
				for (int k = 0; k < 3; k++){
					for (int l = 0; l < 3; l++){
						m[3 * k + l] = transformation.coeff(3 * i + k, 3 * j + l);
					}
					v[k] = nface.m_Vertices[j][k];
				}
				v = m * v;
				double weight = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
				if (weight > threshold){
					for (int k = 0; k < 3; k++)
						for (int l = 0; l < 3; l++)
							tripletList.push_back(T(3 * i + k, 3 * j + l, m[3 * k + l]));
				}
			}
		}
		transformation.setZero();
		transformation.setFromTriplets(tripletList.begin(), tripletList.end());
	}

	void TestThreshold()
	{
		e_index = 2;
		compress_level = 0;
		binary_format = true;
		ReadTransformation(e_index, compress_level);

		compress_level = 25;
		double threshold = pow(0.5, compress_level);
		vector<int> indices = { 24, 579, 236 };
		int vnum = nface.m_Vertices.size();
		for (int ii = 0; ii < indices.size(); ii++){
			int i = indices[ii];
			vector<ml::vec3d> vs(vnum);
			vector<ml::vec4d> cs(vnum);
			for (int j = 0; j < vnum; j++){
				ml::mat3d m;
				ml::vec3d v;
				for (int k = 0; k < 3; k++){
					for (int l = 0; l < 3; l++){
						m[3 * k + l] = transformation.coeff(3 * i + k, 3 * j + l);
					}
					v[k] = nface.m_Vertices[j][k];
				}
				v = m * v;
				vs[j] = v;
				double weight = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
				if (weight > threshold){
					cs[j] = ml::vec4d(0.8, 0, 0, 1);
				}
				else{
					cs[j] = ml::vec4d(0.0, 0, 0, 1);
				}
			}
			ml::MeshDatad tm;
			tm.m_Vertices = vs;
			tm.m_Colors = cs;
			char str[200];
			sprintf(str, "compress_%d_vertex_%d.obj", compress_level, i);
			ml::MeshIOd::saveToOBJ(Test_Output_Dir + str, tm);
		}
	}

	void SolveOne()
	{
		int vnum = nface.m_Vertices.size();
		int fnum = nface.m_FaceIndicesVertices.size();
		int ednum = fnum * 2;

		SpMat g(3 * ednum, 3 * vnum);
		SpMat h(3 * ednum, 3 * ednum);
		SpMat f(3 * vnum, 3 * vnum);

		std::vector<Tripletd> tripletList1, tripletList2, tripletList3;
		tripletList1.reserve(12 * fnum);
		for (int i = 0; i < fnum; i++){
			// group by 6 * 2 -- g
			ml::MeshDatad::Indices::Face& ind = nface.m_FaceIndicesVertices[i];
			tripletList1.push_back(T(6 * i + 0, ind[0] * 3 + 0, -1));
			tripletList1.push_back(T(6 * i + 1, ind[0] * 3 + 1, -1));
			tripletList1.push_back(T(6 * i + 2, ind[0] * 3 + 2, -1));
			tripletList1.push_back(T(6 * i + 3, ind[0] * 3 + 0, -1));
			tripletList1.push_back(T(6 * i + 4, ind[0] * 3 + 1, -1));
			tripletList1.push_back(T(6 * i + 5, ind[0] * 3 + 2, -1));

			tripletList1.push_back(T(6 * i + 0, ind[1] * 3 + 0, 1));
			tripletList1.push_back(T(6 * i + 1, ind[1] * 3 + 1, 1));
			tripletList1.push_back(T(6 * i + 2, ind[1] * 3 + 2, 1));
			tripletList1.push_back(T(6 * i + 3, ind[2] * 3 + 0, 1));
			tripletList1.push_back(T(6 * i + 4, ind[2] * 3 + 1, 1));
			tripletList1.push_back(T(6 * i + 5, ind[2] * 3 + 2, 1));
		}
		g.setZero();
		g.setFromTriplets(tripletList1.begin(), tripletList1.end());

		tripletList2.reserve(18 * fnum);
		for (int i = 0; i < fnum; i++){
			ml::mat3d& s = deformation[i];
			int x = 6 * i;
			for (int iter = 0; iter < 2; iter++){
				for (int j = 0; j < 3; j++){
					for (int k = 0; k < 3; k++){
						tripletList2.push_back(T(x + j, x + k, s.at(j, k)));
					}
				}
				x += 3;
			}
		}
		h.setZero();
		h.setFromTriplets(tripletList2.begin(), tripletList2.end());

		tripletList3.reserve(3 * vnum);
		for (int i = 0; i < vnum; i++){
			// group by 3 * 1 -- f
			if (!mask_reader.face_mask[i]){
				tripletList3.push_back(T(i * 3 + 0, i * 3 + 0, 100));
				tripletList3.push_back(T(i * 3 + 1, i * 3 + 1, 100));
				tripletList3.push_back(T(i * 3 + 2, i * 3 + 2, 100));
			}
		}
		f.setZero();
		f.setFromTriplets(tripletList3.begin(), tripletList3.end());

		Eigen::CholmodSupernodalLLT<SpMat, Eigen::Lower> chol;
		SpMat A = g.transpose() * g + f;
		SpMat rhs(3 * vnum, 3 * vnum);
		rhs = g.transpose() * h * g + f;
		chol.compute(A);
		if (chol.info() != Eigen::Success) {
			std::cout << "Oh: Very bad" << endl;
		}
		for (int i = 0; i < rhs.cols(); i++)
		{
			transformation.col(i) = chol.solve(rhs.col(i));

			if (i % 100 == 0)
				std::cout << i << " cols\n";
		}
	}

	void Solve()
	{
		//for (int i = 0; i < filenames.size(); i++)
		//	cout << i << " " << filenames[i] << "\n";

		//for (e_index = filenames.size() - 1; e_index >= 0; e_index--){
		int start = 0;
		std::cout << "e_index: ";
		std::cin >> start;
		int end = 0;
		std::cin >> end;

		for (e_index = start; e_index >= end; e_index--){
			cout << filenames.size() << "/" << e_index << ": " << filenames[e_index] <<"\n";
			//e_index = 27;
			ml::MeshIOd::loadFromOBJ(filenames[e_index], eface);
			for (int i = 0; i < eface.m_Vertices.size(); i++){
				eface.m_Vertices[i].x *= -1;
				eface.m_Vertices[i].z *= -1;
			}
			//ReadDeformation();
			CalcDeformation();
			SolveOne();

			binary_format = true;
			compress_level = 11;
			Compress();
			Test();
			WriteTransformation();
			return;
		}
	}

public:
	string nfilename;
	vector<string> filenames;
	ml::MeshDatad nface;
	ml::MeshDatad eface;
	int e_index;
	int compress_level;
	bool binary_format;

	vector<ml::mat3d> deformation;
	SpMat transformation;

	MaskReader mask_reader;
};