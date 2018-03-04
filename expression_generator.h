#pragma once

#include "mask_reader.h"
#include "file_func.h"

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

class ExpressionGenerator
{
public:
	double VectorLength(const Vector& v)
	{
		return std::sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
	}

	void BarycenticCoord(Point pt, Point pt1, Point pt2, Point pt3, double& a, double& b)
	{
		Vector v = pt - pt1;
		Plane plane(pt1, pt2, pt3);
		Vector orth = plane.orthogonal_vector();
		orth = orth / VectorLength(orth);
		v = v - (Vector(v.x(), v.y(), v.z()) * orth) * orth;
		double l = VectorLength(v);
		if (l == 0){
			a = 1; b = 0;
			return;
		}

		Vector v2 = pt2 - pt1;
		Vector v3 = pt3 - pt1;
		double l2 = VectorLength(v2);
		double l3 = VectorLength(v3);
		double cos = std::min(1.0, v2 * v3 / l2 / l3);
		double sin = std::sqrt(std::max(0.0, 1 - cos * cos));
		double cosv = std::min(1.0, v2 * v / l2 / l);
		double sinv = std::sqrt(std::max(0.0, 1 - cosv * cosv));
		double x1, x2, x3, y1, y2, y3, x, y;
		x1 = 0;
		y1 = 0;
		x2 = l2;
		y2 = 0;
		x3 = l3 * cos;
		y3 = l3 * sin;
		x = l * cosv;
		y = l * sinv;

		a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3))
			/ ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
		b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3))
			/ ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
	}

	void RunOne(string input_filename, string output_filename)
	{
		//int left_corner = 45;
		//int right_corner = 138;
		static set<int> upper_lip = { 74, 70, 69, 71, 4, 499, 497, 498, 502, /**/};
		static set<int> lower_lip = { 425, 84, 513 };
		upper_lip.clear();
		lower_lip.clear();

		ml::MeshDatad eface1, eface2, nface1, nface2;
		ml::MeshIOd::loadFromOBJ(Data_Input_Dir + "expression/Normal.obj", eface1);
		ml::MeshIOd::loadFromOBJ(Data_Input_Dir + "pca_obj_simplified1/original.obj", nface1);

		ml::MeshIOd::loadFromOBJ(input_filename, eface2);

		int e_vnum = eface1.m_Vertices.size();
		int n_vnum = nface1.m_Vertices.size();
		int e_fnum = eface1.m_FaceIndicesVertices.size();
		int n_fnum = nface1.m_FaceIndicesVertices.size();

		MaskReaderRaw mask_reader;
		std::vector<bool> face_mask(e_vnum, false);
		std::vector<bool> vts_mask[5];
		for (int i = 0; i < 5; i++)
			vts_mask[i] = mask_reader.face_mask;
		for (int iter = 0; iter < 3; iter++){
			for (int i = 0; i < e_fnum; i++){
				ml::MeshDatad::Indices::Face ind = eface1.m_FaceIndicesVertices[i];
				for (int j = 0; j < 4; j++){
					if (vts_mask[iter][ind[j]] == true){
						for (int k = 0; k < 4; k++){
							vts_mask[iter + 1][ind[k]] = true;
						}
						face_mask[i] = true;
						break;
					}
				}
			}
		}

		//ml::MeshDataf mesh;
		//ml::MeshIOf::loadFromOBJ("F:/FaceAnimation/BS_2017_11_29/mask3.obj", mesh);
		//mesh.m_Colors.resize(8957);
		//for (int i = 0; i < 8957; i++)
		//{
		//	if (vts_mask[1][i] == true)
		//		mesh.m_Colors[i] = ml::vec4f(1, 1, 1, 1);
		//	else
		//		mesh.m_Colors[i] = ml::vec4f(0, 1, 0, 1);
		//}
		//ml::MeshIOf::saveToOFF("C:\\Users\\zhx\\Desktop\\zhx.off", mesh);
		//return;

		// collect vertices of mesh for easy access
		std::vector<Point> epts1(e_vnum);
		std::vector<Point> epts2(e_vnum);
		std::vector<Point> npts1(n_vnum);
		for (int i = 0; i < e_vnum; i++){
			ml::vec3d v = eface1.m_Vertices[i];
			epts1[i] = Point(v.x, v.y, v.z);
		}
		for (int i = 0; i < e_vnum; i++){
			ml::vec3d v = eface2.m_Vertices[i];
			epts2[i] = Point(v.x, v.y, v.z);
		}
		for (int i = 0; i < n_vnum; i++){
			ml::vec3d v = nface1.m_Vertices[i];
			npts1[i] = Point(v.x, v.y, v.z);
		}
		// build polyhedron and AABB tree
		Polyhedron P1;
		std::map<Halfedge_handle, int> hh2index1;
		for (int i = 0; i < eface1.m_FaceIndicesVertices.size(); i++){
			ml::MeshDatad::Indices::Face ind = eface1.m_FaceIndicesVertices[i];
			Halfedge_handle hh1 = P1.make_triangle(epts1[ind[0]], epts1[ind[1]], epts1[ind[2]]);
			Halfedge_handle hh2 = P1.make_triangle(epts1[ind[0]], epts1[ind[2]], epts1[ind[3]]);
			hh2index1[hh1] = i * 2;
			hh2index1[hh2] = i * 2 + 1;
		}
		Tree e_tree(faces(P1).begin(), faces(P1).end(), P1); 
		//
		Polyhedron P2;
		ml::MeshDatad lip_mesh;
		std::map<Halfedge_handle, int> hh2index2;
		ml::MeshIOd::loadFromOBJ(Data_Input_Dir + "txt_info/sha.obj", lip_mesh);
		for (int i = 0; i < eface1.m_FaceIndicesVertices.size(); i++){
			ml::MeshDatad::Indices::Face ind = eface1.m_FaceIndicesVertices[i];
			if (lip_mesh.m_Colors[ind[0]][1] == 0
				|| lip_mesh.m_Colors[ind[1]][1] == 0
				|| lip_mesh.m_Colors[ind[2]][1] == 0
				|| lip_mesh.m_Colors[ind[3]][1] == 0){
				Halfedge_handle hh1 = P2.make_triangle(epts1[ind[0]], epts1[ind[1]], epts1[ind[2]]);
				Halfedge_handle hh2 = P2.make_triangle(epts1[ind[0]], epts1[ind[2]], epts1[ind[3]]);
				hh2index2[hh1] = i * 2;
				hh2index2[hh2] = i * 2 + 1;
				break;
			}
		}
		Tree upper_tree(faces(P2).begin(), faces(P2).end(), P2);
		//
		Polyhedron P3;
		std::map<Halfedge_handle, int> hh2index3;
		lip_mesh.clear();
		ml::MeshIOd::loadFromOBJ(Data_Input_Dir + "txt_info/xia.obj", lip_mesh);
		for (int i = 0; i < eface1.m_FaceIndicesVertices.size(); i++){
			ml::MeshDatad::Indices::Face ind = eface1.m_FaceIndicesVertices[i];
			if (lip_mesh.m_Colors[ind[0]][1] == 0
				|| lip_mesh.m_Colors[ind[1]][1] == 0
				|| lip_mesh.m_Colors[ind[2]][1] == 0
				|| lip_mesh.m_Colors[ind[3]][1] == 0){
				Halfedge_handle hh1 = P3.make_triangle(epts1[ind[0]], epts1[ind[1]], epts1[ind[2]]);
				Halfedge_handle hh2 = P3.make_triangle(epts1[ind[0]], epts1[ind[2]], epts1[ind[3]]);
				hh2index3[hh1] = i * 2;
				hh2index3[hh2] = i * 2 + 1;
				break;
			}
		}
		Tree lower_tree(faces(P3).begin(), faces(P3).end(), P3);

		//
		std::vector<bool> fm(n_vnum, true);
		std::vector<bool> hm(n_vnum, false);
		//
		nface2.m_Vertices.resize(n_vnum, ml::vec3d(0, 0, 0));
		nface2.m_FaceIndicesVertices = nface1.m_FaceIndicesVertices;
		for (int i = 0; i < n_vnum; i++){

			int index;
			Polyhedron::Face_handle ff;
			if (upper_lip.find(i) != upper_lip.end()){
				Point_and_primitive_id pp = upper_tree.closest_point_and_primitive(npts1[i]);
				ff = pp.second;
				index = hh2index2[ff->halfedge()];
			}
			else if (lower_lip.find(i) != lower_lip.end()){
				Point_and_primitive_id pp = lower_tree.closest_point_and_primitive(npts1[i]);
				ff = pp.second;
				index = hh2index3[ff->halfedge()];
			}
			else{
				Point_and_primitive_id pp = e_tree.closest_point_and_primitive(npts1[i]);
				ff = pp.second;
				index = hh2index1[ff->halfedge()];
			}
			ml::MeshDatad::Indices::Face ind = eface1.m_FaceIndicesVertices[index / 2];

			////
			if (face_mask[index / 2] == false)
				fm[i] = false;
			for (int j = 0; j < 4; j++)
				if (mask_reader.half_mask[ind[j]])
					hm[i] = true;

			// remain the same if not masked
			if (!face_mask[index / 2]){
				nface2.m_Vertices[i] = ml::vec3d(npts1[i].x(), npts1[i].y(), npts1[i].z());
				continue;
			}

			// calculate barycentic coord
			double a, b, c;
			BarycenticCoord(npts1[i],
				ff->halfedge()->vertex()->point(),
				ff->halfedge()->next()->vertex()->point(),
				ff->halfedge()->next()->next()->vertex()->point(),
				a, b);
			c = 1 - a - b;

			// calculate final coord
			if (index % 2 == 0){
				double x = a * epts2[ind[0]].x() + b * epts2[ind[1]].x() + c * epts2[ind[2]].x();
				double y = a * epts2[ind[0]].y() + b * epts2[ind[1]].y() + c * epts2[ind[2]].y();
				double z = a * epts2[ind[0]].z() + b * epts2[ind[1]].z() + c * epts2[ind[2]].z();
				nface2.m_Vertices[i] = ml::vec3d(x, y, z);
			}
			else{
				double x = a * epts2[ind[0]].x() + b * epts2[ind[2]].x() + c * epts2[ind[3]].x();
				double y = a * epts2[ind[0]].y() + b * epts2[ind[2]].y() + c * epts2[ind[3]].y();
				double z = a * epts2[ind[0]].z() + b * epts2[ind[2]].z() + c * epts2[ind[3]].z();
				nface2.m_Vertices[i] = ml::vec3d(x, y, z);
			}
			
			//if (i == 45)
			//	nface2.m_Vertices[i] = ml::vec3d(epts2[2514][0], epts2[2514][1], epts2[2514][2]);
			//else if (i == 138)
			//	nface2.m_Vertices[i] = ml::vec3d(epts2[5917][0], epts2[5917][1], epts2[5917][2]);
		}

		MaskReader tmp;
		tmp.face_mask = fm;
		tmp.half_mask = hm;
		tmp.Write(fm, hm);
		tmp.Test(Data_Input_Dir + "pca_obj_simplified1/original.obj");

		//for (int i = 0; i < n_vnum; i++){
		//	ml::vec3d& v = nface2.m_Vertices[i];
		//	v = v / 1000.0;
		//	v.z = v.z * -1;
		//}
		ml::MeshIOd::saveToOBJ(output_filename, nface2);

	}

	void Run()
	{
		vector<string> filenames;
		GetAllFormatFiles(Data_Input_Dir + "expression", filenames, "obj");

		int count = 0;
		for (int i = 0; i < filenames.size(); i++)
		{
			//i = 5;
			if (filenames[i] == Data_Input_Dir + "expression/Normal.obj")
				continue;
			char str[100];
			sprintf(str, "pca_obj_simplified1/expression/%d.obj", count++);
			RunOne(filenames[i], Data_Input_Dir + str);

			std::cout << count - 1 << ".obj ~ " << filenames[i] << "\n";
			//break;
		}
		//system("pause");
	}
private:

};