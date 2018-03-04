#pragma once

#include "mask_reader.h"
#include "file_func.h"

#include "deformation_transfer.h"

#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Eigen>
#include <Eigen/Sparse>
//#include <Eigen/PaStiXSupport>
#include <Eigen/CholmodSupport>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

//Eigen
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;
//CGAL
typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Vector_3 Vector;
typedef K::Triangle_3 Triangle;
typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
//AABB tree
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

class PcaGenerator3D
{
public:
	PcaGenerator3D()
	{
		ml::MeshIOd::loadFromOBJ(Data_Input_Dir + "pca_obj/0.obj", nface1);
		ml::MeshIOd::loadFromOBJ(Data_Input_Dir + "pca_obj_simplified1/original.obj", sface1);

		n_vnum = nface1.m_Vertices.size();
		s_vnum = sface1.m_Vertices.size();
		n_fnum = nface1.m_FaceIndicesVertices.size();
		s_fnum = sface1.m_FaceIndicesVertices.size();

		//
		mask_reader.Read(s_vnum);
		//mask_reader.Test(Data_Input_Dir + "pca_obj_simplified1/original.obj");
		//system("pause");

		// collect vertices of mesh for easy access
		std::vector<Point> npts1(n_vnum);
		std::vector<Point> spts1(s_vnum);
		for (int i = 0; i < n_vnum; i++){
			ml::vec3d v = nface1.m_Vertices[i];
			npts1[i] = Point(v.x, v.y, v.z);
		}
		for (int i = 0; i < s_vnum; i++){
			ml::vec3d v = sface1.m_Vertices[i];
			spts1[i] = Point(v.x, v.y, v.z);
		}
		// build polyhedron and AABB tree
		Polyhedron P1;
		std::map<Halfedge_handle, int> hh2index1;
		for (int i = 0; i < nface1.m_FaceIndicesVertices.size(); i++){
			ml::MeshDatad::Indices::Face ind = nface1.m_FaceIndicesVertices[i];
			Halfedge_handle hh1 = P1.make_triangle(npts1[ind[0]], npts1[ind[1]], npts1[ind[2]]);
			hh2index1[hh1] = i * 2;
		}
		Tree e_tree(faces(P1).begin(), faces(P1).end(), P1);

		indices.resize(s_vnum);
		for (int i = 0; i < s_vnum; i++){

			int index;
			Polyhedron::Face_handle ff;
			Point_and_primitive_id pp = e_tree.closest_point_and_primitive(spts1[i]);
			ff = pp.second;
			index = hh2index1[ff->halfedge()];
			//ml::MeshDatad::Indices::Face ind = nface1.m_FaceIndicesVertices[index / 2];

			// remain the same if not masked
			if (!mask_reader.face_mask[i])
				continue;

			indices[i] = index / 2;
		}
	}

	void RunOne(string input_filename, string output_filename)
	{
		ml::MeshDatad nface2, sface2;
		ml::MeshIOd::loadFromOBJ(input_filename, nface2);
		sface2.m_Vertices = sface1.m_Vertices;
		sface2.m_FaceIndicesVertices = sface1.m_FaceIndicesVertices;

		DeformationTransfer dt;
		dt.Solve(nface1, nface2);

		for (int i = 0; i < s_vnum; i++){
			// remain the same if not masked
			if (!mask_reader.face_mask[i])
				continue;

			ml::MeshDatad::Indices::Face ind = nface2.m_FaceIndicesVertices[indices[i]];
			ml::vec3d p0 = nface1.m_Vertices[ind[0]];
			ml::vec3d p00 = nface2.m_Vertices[ind[0]];
			ml::vec3d p1 = sface1.m_Vertices[i];
			sface2.m_Vertices[i] = dt.tar_deformation[indices[i]] * (p1 - p0) + p00;
		}

		ml::MeshIOd::saveToOBJ(output_filename, sface2);
	}

	void Run()
	{
		for (int i = 1; i <= 50; i++)
		{
			char str1[100], str2[100];
			sprintf(str1, "pca_obj/%d.obj", i);
			sprintf(str2, "pca_obj_simplified1/pca/%d.obj", i);

			RunOne(Data_Input_Dir + str1, Data_Input_Dir + str2);
			std::cout << str2 << "\n";
			//break;
		}
		//system("pause");
	}


private:
	MaskReader mask_reader;

	int n_vnum, s_vnum, n_fnum, s_fnum;
	ml::MeshDatad nface1, sface1;
	std::vector<int> indices;
};