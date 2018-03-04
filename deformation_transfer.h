#pragma once

#include <Eigen/Eigen>
#include <Eigen/Sparse>

//#include "./mlibutil/mlibutil.h"
#include "defines.h"

//#include<CGAL/Simple_cartesian.h>
//#include<CGAL/Polyhedron_incremental_builder_3.h>
//#include<CGAL/Polyhedron_3.h>
//#include<CGAL/IO/Polyhedron_iostream.h>
//typedef CGAL::Simple_cartesian<double>     Kernel;
//typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
//typedef Polyhedron::HalfedgeDS             HalfedgeDS;

#include <map>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
	std::vector<double> &coords;
	std::vector<int>    &tris;
	polyhedron_builder(std::vector<double> &_coords, std::vector<int> &_tris) : coords(_coords), tris(_tris) {}
	void operator()(HDS& hds) {
		typedef typename HDS::Vertex   Vertex;
		typedef typename Vertex::Point Point;

		// create a cgal incremental builder
		CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
		B.begin_surface(coords.size() / 3, tris.size() / 3);

		// add the polyhedron vertices
		for (int i = 0; i<(int)coords.size(); i += 3){
			B.add_vertex(Point(coords[i + 0], coords[i + 1], coords[i + 2]));
		}

		// add the polyhedron triangles
		for (int i = 0; i<(int)tris.size(); i += 3){
			B.begin_facet();
			B.add_vertex_to_facet(tris[i + 0]);
			B.add_vertex_to_facet(tris[i + 1]);
			B.add_vertex_to_facet(tris[i + 2]);
			B.end_facet();
		}

		// finish up the surface
		B.end_surface();
	}
};

class DeformationTransfer
{
public:
	// ceres-error class
	//class MatrixError
	//{
	//public:
	//	MatrixError(int index, vector<ml::mat3d>& src_deformation)
	//		:index(index), src_deformation(src_deformation)
	//	{}

	//	template <class T>
	//	bool operator()(const T* const deformation, T* residuals) const {
	//		for (int i = 0; i < 9; i++){
	//			residuals[i] = deformation[i] - (double)src_deformation[index][i];
	//		}
	//		return true;
	//	}

	//	static ceres::CostFunction* Create(int index, vector<ml::mat3d>& src_deformation) {
	//		// first residual dimension, followed with parameters' dimensions
	//		return (new ceres::AutoDiffCostFunction<MatrixError, 9, 9>(
	//			new MatrixError(index, src_deformation)));
	//	}
	//private:
	//	int index;
	//	vector<ml::mat3d>& src_deformation;
	//};

	//class ConstraintError
	//{
	//public:
	//	ConstraintError(int i1, int i2, int i3, ml::MeshDatad& src_mesh)
	//		:vertex_index(i1), face_index1(i2), face_index2(i3), src_mesh(src_mesh)
	//	{}

	//	template <class T>
	//	bool operator()(const T* const d1, const T* const t1, const T* const d2, const T* const t2, T* residuals) const {
	//		T v[3], v1[3], v2[3];
	//		for (int i = 0; i < 3; i++)
	//			v[i] = (T)src_mesh.m_Vertices[vertex_index][i];
	//		for (int i = 0; i < 3; i++){
	//			v1[i] = t1[i];
	//			for (int j = 0; j < 3; j++)
	//				v1[i] += d1[i * 3 + j] * v[j];
	//		}
	//		for (int i = 0; i < 3; i++){
	//			v2[i] = t2[i];
	//			for (int j = 0; j < 3; j++)
	//				v2[i] += d2[i * 3 + j] * v[j];
	//		}
	//		for (int i = 0; i < 3; i++){
	//			residuals[i] = 1000.0 * (v1[i] - v2[i]);
	//			//cout << *(double*)(&residuals[i]) << " ";
	//		}
	//		//cout << "\n";
	//		return true;
	//	}

	//	static ceres::CostFunction* Create(int i1, int i2, int i3, ml::MeshDatad& src_mesh) {
	//		// first residual dimension, followed with parameters' dimensions
	//		return (new ceres::AutoDiffCostFunction<ConstraintError, 3, 9, 3, 9, 3>(
	//			new ConstraintError(i1, i2, i3, src_mesh)));
	//	}
	//private:
	//	int vertex_index, face_index1, face_index2;
	//	ml::MeshDatad& src_mesh;
	//};

	DeformationTransfer() {}

	DeformationTransfer(string src_file, string tar_file)
	{
		// surport OBJ OFF format
		if (src_file.find(".obj") != -1)
		{
			ml::MeshIOd::loadFromOBJ(src_file, src_mesh);
			ml::MeshIOd::loadFromOBJ(tar_file, tar_mesh);
		}
		else if (src_file.find(".off") != -1)
		{
			ml::MeshIOd::loadFromOFF(src_file, src_mesh);
			ml::MeshIOd::loadFromOFF(tar_file, tar_mesh);
		}

		//// coord and tris 
		//vector<double> coords;
		//vector<int> tris;
		//coords.resize(3 * src_mlmesh.m_Vertices.size());
		//tris.resize(3 * src_mlmesh.m_FaceIndicesVertices.size());
		//// build source mesh
		//for (int i = 0; i < src_mlmesh.m_Vertices.size(); i++)
		//{
		//	coords[3 * i + 0] = src_mlmesh.m_Vertices[i].x;
		//	coords[3 * i + 1] = src_mlmesh.m_Vertices[i].y;
		//	coords[3 * i + 2] = src_mlmesh.m_Vertices[i].z;
		//}
		//for (int i = 0; i < src_mlmesh.m_FaceIndicesVertices.size(); i++)
		//{
		//	tris[3 * i + 0] = src_mlmesh.m_FaceIndicesVertices[i][0];
		//	tris[3 * i + 1] = src_mlmesh.m_FaceIndicesVertices[i][1];
		//	tris[3 * i + 2] = src_mlmesh.m_FaceIndicesVertices[i][2];
		//}
		//polyhedron_builder<HalfedgeDS> b1(coords, tris);
		//src_mesh.delegate(b1);
		//// build target mesh
		//for (int i = 0; i < src_mlmesh.m_Vertices.size(); i++)
		//{
		//	coords[3 * i + 0] = src_mlmesh.m_Vertices[i].x;
		//	coords[3 * i + 1] = src_mlmesh.m_Vertices[i].y;
		//	coords[3 * i + 2] = src_mlmesh.m_Vertices[i].z;
		//}
		//for (int i = 0; i < src_mlmesh.m_FaceIndicesVertices.size(); i++)
		//{
		//	tris[3 * i + 0] = src_mlmesh.m_FaceIndicesVertices[i][0];
		//	tris[3 * i + 1] = src_mlmesh.m_FaceIndicesVertices[i][1];
		//	tris[3 * i + 2] = src_mlmesh.m_FaceIndicesVertices[i][2];
		//}
		//polyhedron_builder<HalfedgeDS> b2(coords, tris);
		//tar_mesh.delegate(b1);

		for (int i = 0; i < src_mesh.m_FaceIndicesVertices.size(); i++)
		{
			vertex2face.insert(pair<int, int>(src_mesh.m_FaceIndicesVertices[i][0], i));
			vertex2face.insert(pair<int, int>(src_mesh.m_FaceIndicesVertices[i][1], i));
			vertex2face.insert(pair<int, int>(src_mesh.m_FaceIndicesVertices[i][2], i));
		}

		//
		translation.resize(src_mesh.m_FaceIndicesVertices.size());
		src_deformation.resize(src_mesh.m_FaceIndicesVertices.size());
		tar_deformation.resize(src_mesh.m_FaceIndicesVertices.size());
	}

	ml::mat3d CalcTetrahedon(ml::MeshDatad& mesh, ml::MeshDatad& mesh1, int ii)
	{
		ml::MeshDatad::Indices::Face& ind = mesh1.m_FaceIndicesVertices[ii];
		ml::vec3d v1 = mesh.m_Vertices[ind[0]];
		ml::vec3d v2 = mesh.m_Vertices[ind[1]];
		ml::vec3d v3 = mesh.m_Vertices[ind[2]];
		ml::vec3d vtmp = ml::vec3d::cross(v2 - v1, v3 - v1);
		ml::vec3d v4 = v1 + vtmp / vtmp.length();

		ml::mat3d s(v2 - v1, v3 - v1, v4 - v1);
		s.transpose();

		return s;
	}

	// build problem and solve for the translation and deformation from source mesh to the target
	void Solve(ml::MeshDatad& src_mesh1, ml::MeshDatad& tar_mesh1)
	{
		////
		translation.resize(src_mesh1.m_FaceIndicesVertices.size());
		src_deformation.resize(src_mesh1.m_FaceIndicesVertices.size());
		tar_deformation.resize(src_mesh1.m_FaceIndicesVertices.size());

		// initialize
		for (int i = 0; i < src_deformation.size(); i++)
		{
			////
			ml::MeshDatad::Indices::Face& ind = src_mesh1.m_FaceIndicesVertices[i];
			ml::vec3d& v11 = src_mesh1.m_Vertices[ind[0]];
			ml::vec3d& v12 = src_mesh1.m_Vertices[ind[1]];
			ml::vec3d& v13 = src_mesh1.m_Vertices[ind[2]];
			ml::vec3d& v21 = tar_mesh1.m_Vertices[ind[0]];
			ml::vec3d& v22 = tar_mesh1.m_Vertices[ind[1]];
			ml::vec3d& v23 = tar_mesh1.m_Vertices[ind[2]];
			double dist = (v11 - v21).lengthSq() + (v12 - v22).lengthSq() + (v13 - v23).lengthSq();
			if (dist < 0.0001){
				src_deformation[i].setIdentity();
			}
			ml::mat3d v1 = CalcTetrahedon(src_mesh1, src_mesh1, i);
			ml::mat3d v2 = CalcTetrahedon(tar_mesh1, src_mesh1, i);
			src_deformation[i] = v2 * v1.getInverse();

			translation[i] = tar_mesh1.m_Vertices[src_mesh1.m_FaceIndicesVertices[i][0]]
				- src_mesh1.m_Vertices[src_mesh1.m_FaceIndicesVertices[i][0]];
		}
		tar_deformation = src_deformation;
		return;

		/*ceres::Problem problem;
		ceres::Solver::Options options;
		options.linear_solver_type = ceres::SPARSE_SCHUR;
		options.minimizer_progress_to_stdout = true;
		options.max_num_iterations = 20;
		options.num_threads = 4;
		for (int i = 0; i < tar_deformation.size(); i++){
			problem.AddResidualBlock(
				MatrixError::Create(i, src_deformation),
				0, tar_deformation[i].getData()
				);
		}
		for (int i = 0; i < src_mesh.m_Vertices.size(); i++){
			auto ret = vertex2face.equal_range(i);
			vector<int> vret;
			for (auto it = ret.first; it != ret.second; ++it){
				vret.push_back(it->second);
			}
			for (int j = 0; j < vret.size(); j++){
				for (int k = j + 1; k < vret.size(); k++){
					problem.AddResidualBlock(
						ConstraintError::Create(i, vret[j], vret[k], src_mesh),
						0, 
						tar_deformation[vret[j]].getData(), 
						translation[vret[j]].ptr(),
						tar_deformation[vret[k]].getData(), 
						translation[vret[k]].ptr()
						);
				}
			}
		}
		ceres::Solver::Summary summary;
		ceres::Solve(options, &problem, &summary);*/
	}

	// write out the result polyhedron
	// e.g. tar_deformation * face@source + translation
	void Test()
	{
		int fnum = src_mesh.m_FaceIndicesVertices.size();
		int vnum = 3 * fnum;
		tar_mesh.m_Vertices.resize(vnum);
		//tar_mesh.m_FaceIndicesVertices.resize(fnum);
		for (int i = 0; i < fnum; i++){
			std::vector<unsigned int> face(3);
			for (int j = 0; j < 3; j++){
				tar_mesh.m_Vertices[3 * i + j] = tar_deformation[i] * src_mesh.m_Vertices[src_mesh.m_FaceIndicesVertices[i][j]] + translation[i];
				face[j] = 3 * i + j;
			}
			tar_mesh.m_FaceIndicesVertices.push_back(face);
		}
		ml::MeshIOd::saveToOBJ(Test_Output_Dir + "test_deformation.obj", tar_mesh);
	}

	void WriteDeformation()
	{
		int fnum = src_mesh.m_FaceIndicesVertices.size();
		ofstream ofs(Test_Output_Dir + "deformation.txt");
		for (int i = 0; i < fnum; i++){
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					ofs << tar_deformation[i][j * 3 + k] << " ";
				}
			}
			ofs << "\n";
		}
		ofs.close();
	}

public:
	ml::MeshDatad src_mesh, tar_mesh;
	multimap<int, int> vertex2face;
	//Polyhedron src_mesh, tar_mesh;
	vector<ml::vec3d> translation;
	vector<ml::mat3d> src_deformation;
	vector<ml::mat3d> tar_deformation;
};