#pragma once

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include "./mlibutil/mlibutil.h"

#include <string>
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Tripletd;


template<class Matrix>
void write_binary(string filename, const Matrix& matrix){
	std::ofstream out(filename, ios::out | ios::binary | ios::trunc);
	typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
	out.write((char*)(&rows), sizeof(typename Matrix::Index));
	out.write((char*)(&cols), sizeof(typename Matrix::Index));
	std::cout << rows << " " << cols << "\n";
	out.write((char*)matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
	out.close();
}
template<class Matrix>
void read_binary(string filename, Matrix& matrix){
	std::ifstream in(filename, ios::in | std::ios::binary);
	typename Matrix::Index rows = 0, cols = 0;
	in.read((char*)(&rows), sizeof(typename Matrix::Index));
	in.read((char*)(&cols), sizeof(typename Matrix::Index));
	matrix.resize(rows, cols);
	//std::cout << rows << " " << cols << "\n";
	in.read((char *)matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
	in.close();
}

//
template <typename TT>
void Serialize(SpMat& m, string filename) {
	std::vector<T> res;
	int sz = m.nonZeros();
	m.makeCompressed();
	fstream writeFile;
	writeFile.open(filename, ios::binary | ios::out);
	if (writeFile.is_open())
	{
		int rows, cols, nnzs, outS, innS;
		rows = m.rows();
		cols = m.cols();
		nnzs = m.nonZeros();
		outS = m.outerSize();
		innS = m.innerSize();
		writeFile.write((const char *)&(rows), sizeof(int));
		writeFile.write((const char *)&(cols), sizeof(int));
		writeFile.write((const char *)&(nnzs), sizeof(int));
		writeFile.write((const char *)&(outS), sizeof(int));
		writeFile.write((const char *)&(innS), sizeof(int));
		writeFile.write((const char *)(m.valuePtr()), sizeof(TT) * m.nonZeros());
		writeFile.write((const char *)(m.outerIndexPtr()), sizeof(int) * m.outerSize());
		writeFile.write((const char *)(m.innerIndexPtr()), sizeof(int) * m.nonZeros());
		writeFile.close();
	}
}

template <typename TT>
void Deserialize(SpMat& m, string filename) {
	fstream readFile;
	readFile.open(filename, ios::binary | ios::in);
	if (readFile.is_open())
	{
		int rows, cols, nnz, inSz, outSz;
		readFile.read((char*)&rows, sizeof(int));
		readFile.read((char*)&cols, sizeof(int));
		readFile.read((char*)&nnz, sizeof(int));
		readFile.read((char*)&inSz, sizeof(int));
		readFile.read((char*)&outSz, sizeof(int));
		m.resize(rows, cols);
		m.makeCompressed();
		m.resizeNonZeros(nnz);
		readFile.read((char*)(m.valuePtr()), sizeof(TT) * nnz);
		readFile.read((char*)(m.outerIndexPtr()), sizeof(int) * outSz);
		readFile.read((char*)(m.innerIndexPtr()), sizeof(int) * nnz);
		m.finalize();
		readFile.close();
	} // file is open
}

template <class T>
void UpdateMeshVertex(T& v, ml::MeshDatad& mesh)
//void UpdateMeshVertex(Eigen::Matrix<double, Eigen::Dynamic, 1>& v, ml::MeshDatad& mesh)
{
	int vnum = v.rows() / 3;
	mesh.m_Vertices.resize(vnum);
	for (int i = 0; i < vnum; i++){
		for (int j = 0; j < 3; j++){
			mesh.m_Vertices[i][j] = v(i * 3 + j);

			if (j == 2)
				mesh.m_Vertices[i][j] += 0.7;
		}
	}
}