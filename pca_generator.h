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

class PcaGenerator
{
public:
	struct Coeff
	{
		double a, b, c, index;

		Coeff()
		{}

		Coeff(double a, double b, double c, double index)
			:a(a), b(b), c(c), index(index)
		{}
	};

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

	void Init()
	{
		ml::MeshDatad nface1;
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

		coeff.resize(s_vnum, Coeff(-1, -1, -1, -1));
		for (int i = 0; i < s_vnum; i++){

			int index;
			Polyhedron::Face_handle ff;
			Point_and_primitive_id pp = e_tree.closest_point_and_primitive(spts1[i]);
			ff = pp.second;
			index = hh2index1[ff->halfedge()];
			ml::MeshDatad::Indices::Face ind = nface1.m_FaceIndicesVertices[index / 2];

			// remain the same if not masked
			if (!mask_reader.face_mask[i])
				continue;

			// calculate barycentic coord
			double a, b, c;
			BarycenticCoord(spts1[i],
				ff->halfedge()->vertex()->point(),
				ff->halfedge()->next()->vertex()->point(),
				ff->halfedge()->next()->next()->vertex()->point(),
				a, b);
			c = 1 - a - b;

			// calculate final coord
			if (index % 2 == 0){
				coeff[i] = Coeff(a, b, c, index / 2);
			}
		}
	}

	void RunOne(string input_filename, string output_filename)
	{
		ml::MeshDatad nface2, sface2;
		ml::MeshIOd::loadFromOBJ(input_filename, nface2);
		sface2.m_Vertices = sface1.m_Vertices;
		sface2.m_FaceIndicesVertices = sface1.m_FaceIndicesVertices;

		for (int i = 0; i < s_vnum; i++){
			// remain the same if not masked
			if (!mask_reader.face_mask[i])
				continue;

			ml::MeshDatad::Indices::Face ind = nface2.m_FaceIndicesVertices[coeff[i].index];
			sface2.m_Vertices[i] = coeff[i].a * nface2.m_Vertices[ind[0]]
				+ coeff[i].b * nface2.m_Vertices[ind[1]]
				+ coeff[i].c * nface2.m_Vertices[ind[2]];
		}

		ml::MeshIOd::saveToOBJ(output_filename, sface2);
	}

	void Run()
	{
		Init();
		for (int i = 1; i <= 50; i++)
		{
			char str1[100], str2[100];
			sprintf(str1, "pca_obj/%d.obj", i);
			sprintf(str2, "pca_obj_simplified1/pca/%d.obj", i);

			RunOne(Data_Input_Dir + str1, Data_Input_Dir + str2);
			std::cout << str2 << "\n";
			break;
		}
		//system("pause");
	}
private:
	MaskReader mask_reader;
	int n_vnum, s_vnum, n_fnum, s_fnum;
	ml::MeshDatad sface1;
	std::vector<Coeff> coeff;
};