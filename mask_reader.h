#pragma once

#include "./mlibutil/mlibutil.h"
#include "defines.h"
#include "parameters.h"

#include <string>
#include <fstream>
#include <vector>
using namespace std;

class MaskReaderRaw
{
public:
	MaskReaderRaw()
	{
		int size = 8957;
		face_mask.resize(size, false);
		half_mask.resize(size, false);

		fstream ifs1(Data_Input_Dir + "txt_info/mask1.txt");
		for (int i = 0; i < size; i++)
		{	
			string line;
			getline(ifs1, line);
			istringstream iss(line);
			double tmp, r, g, b, a;
			iss >> tmp >> tmp >> tmp >> tmp >> r >> g >> b >> a;
			if (r < 0.9 || g < 0.9 || b < 0.9)
				face_mask[i] = true;
		}


		fstream ifs2(Data_Input_Dir + "txt_info/mask2.txt");
		//ml::MeshDataf mesh;
		//ml::MeshIOf::loadFromOBJ("F:/FaceAnimation/BS_2017_11_29/mask3.obj", mesh);
		//mesh.m_Colors.resize(size);
		for (int i = 0; i < size; i++)
		{
			string line;
			getline(ifs2, line);
			istringstream iss(line);
			double tmp, r, g, b, a;
			iss >> tmp >> tmp >> tmp >> tmp >> r >> g >> b >> a;
			if (r < 0.9 || g < 0.9 || b < 0.9)
				half_mask[i] = true;
			//mesh.m_Colors[i] = ml::vec4f(r, g, b, a);
		}
		//ml::MeshIOf::saveToOFF("C:\\Users\\zhx\\Desktop\\mask3.off", mesh);
	}

public:
	vector<bool> face_mask;
	vector<bool> half_mask;
};

class MaskReader
{
public:
	void Read(int vnum)
	{
		face_mask.resize(vertex_size);
		half_mask.resize(vertex_size);

		ifstream ifs;
		ifs.open(Data_Input_Dir + "pca_obj_simplified1/mask");

		int tmp;
		for (int i = 0; i < vnum; i++)
		{
			ifs >> tmp;
			face_mask[i] = tmp > 0;
		}
		for (int i = 0; i < vnum; i++)
		{
			ifs >> tmp;
			half_mask[i] = tmp > 0;
		}
		ifs.close();
	}

	void Write(vector<bool>& fm, vector<bool>& hm)
	{
		ofstream ofs;
		ofs.open(Data_Input_Dir + "pca_obj_simplified1/mask");
		for (int i = 0; i < fm.size(); i++)
			ofs << fm[i] << " ";
		for (int i = 0; i < hm.size(); i++)
			ofs << hm[i] << " ";
		ofs.close();
	}

	void Test(string filename)
	{
		ml::MeshDataf mesh;
		ml::MeshIOf::loadFromOBJ(filename, mesh);
		mesh.m_Colors.resize(mesh.m_Vertices.size());
		for (int i = 0; i < mesh.m_Vertices.size(); i++)
		{
			if (face_mask[i] == true)
				mesh.m_Colors[i] = ml::vec4f(1, 1, 1, 1);
			else
				mesh.m_Colors[i] = ml::vec4f(0, 1, 0, 1);
		}
		ml::MeshIOf::saveToOFF(Test_Output_Dir + "mask_reader_test1.off", mesh);

		ml::MeshIOf::loadFromOBJ(filename, mesh);
		mesh.m_Colors.resize(mesh.m_Vertices.size());
		for (int i = 0; i < mesh.m_Vertices.size(); i++)
		{
			if (half_mask[i] == true)
				mesh.m_Colors[i] = ml::vec4f(1, 1, 1, 1);
			else
				mesh.m_Colors[i] = ml::vec4f(0, 1, 0, 1);
		}
		ml::MeshIOf::saveToOFF(Test_Output_Dir + "mask_reader_test2.off", mesh);

		return;
	}

	vector<bool> face_mask;
	vector<bool> half_mask;
};