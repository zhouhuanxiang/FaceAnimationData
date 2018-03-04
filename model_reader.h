#pragma once

#include "defines.h"
#include "parameters.h"
#include "eigen_binary_io.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
using namespace Eigen;
typedef SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Triplet<double> Trupletd;

#include <vector>
#include <string>
#include <fstream>
using namespace std;

class ModelReader
{
public:
	ModelReader(MatrixXd& pca_models_,
		vector<SpMat>& exp_tranforms_,
		VectorXd& pca_weights_)
	{
		printf("pca model read begin...\n");
		read_binary<MatrixXd>(Data_Input_Dir + "pca_obj_simplified1/pca/pca", pca_models_);
		if (pca_models_.rows() != vertex_size * 3){
			printf("bad pca model\n");
			system("pause");
		}
		printf("pca model read done!\n\n");

		printf("pca weights read begin...\n");
		ifstream ifs1;
		ifs1.open(Data_Input_Dir + "pca_weight.txt");
		pca_weights_.resize(pca_size);
		for (int i = 0; i < pca_size; i++){
			double tmp;
			ifs1 >> tmp;
			//pca_weights_(i) = sqrt(53149.0) / pow(tmp, 0.1) / 1000 * 0.2;
			pca_weights_(i) = 50 / pow(tmp, 2) / 1000;
			std::cout << pca_weights_(i) << " ";
		}
		std::cout << "\n";
		//system("pause");
		ifs1.close();
		printf("pca weights read done!\n\n");

		printf("exp transform read begin...\n");
		for (int i = 1; i < exp_size; i++){
			char str[100];
			sprintf(str, "transformation%d_9", i - 1);
			Deserialize<double>(exp_tranforms_[i], Data_Input_Dir + "pca_obj_simplified1/expression/" + str);
			if (exp_tranforms_[i].cols() != vertex_size * 3){
				printf("bad exp transform %d\n", i);
				system("pause");
			}
		}
		printf("exp transform read done!\n\n");
	}

private:
};