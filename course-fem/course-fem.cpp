#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <list>
#include <math.h>  

using namespace std;

// Матрица в строчном формате для ЛОС
struct CRSMatrix {

	vector<int> ig, jg;
	vector<double> ggl, ggu, di;
	vector<double> l, u, d;

	vector<double> F, x, temp, temp0;

	vector<double> r, z, p;

	int n, maxiter;
	double epsilon;

};

// Структура для хранения локальных матриц в плотном формате
struct DenseMatrix {

	vector<vector<double>> G;
	vector<vector<double>> M;
	vector<vector<double>> A;
	vector<vector<double>> Acond3;
	vector<int> globNumVert;
	vector<double> b;

	DenseMatrix() {
		G.resize(3);
		M.resize(3);
		A.resize(3);
		Acond3.resize(2);
		globNumVert.resize(2);

		for (int i = 0; i < 3; i++) {
			G[i].resize(3);
			M[i].resize(3);
			A[i].resize(3);
		}

		for (int i = 0; i < 2; i++) {
			Acond3[i].resize(2);
		}

		b.resize(3);
	}
};

struct Node {

	double r, phi, t;
	int globalNumber;

	Node(double rr, double phii, double tt) {
		r = rr;
		phi = phii;
		t = tt;
	}
	
	Node() {}
};

// Ввод упорядоченный по локальным базисным векторам
struct FinitElement {

	vector<Node> nodes;
	int formulaNumber;

};

struct BoundCond1 {

	int globNum1, globNum2, formulaNumber;

};

struct BoundCond2 {

	int globNum1, globNum2, formulaNumber;

};

struct BoundCond3 {

	int globNum1, globNum2, formulaNumber;

};

struct Grid {

	vector<Node> nodes;
	vector<FinitElement> finitElements;
	vector<BoundCond1> boundConds1;
	vector<BoundCond2> boundConds2;
	vector<BoundCond3> boundConds3;

};

double F(Node& node, int formulaNumber) {

	//return 3 * (node.r) ;
	return 2 * node.r + 2 * node.t * node.r;
	//return -6 * node.r * node.phi + node.r * node.r * node.phi;
	//return (formulaNumber == 0) ? -20 : 0;

}

double lambda(Node& node, int formulaNumber) {

	return 0;
	//return (formulaNumber == 0) ? 10 : 1;

}

double gamma(int formulaNumber) {

	//return (formulaNumber == 0) ? 5 : 0;
	return 0;

}

double xi(int formulaNumber) {

	//return (formulaNumber == 0) ? 5 : 0;
	return 1;

}

double sigma(int formulaNumber) {

	//return (formulaNumber == 0) ? 5 : 0;
	return 1;

}

double tetta(Node& node, int formulaNumber) {

	//return (formulaNumber == 0) ? -6 : ((formulaNumber == 1) ? -1 : 6);
	return 0;

}

double betaF() {

	//return 10;
	return 0;
}

double uBeta(Node& node) {

	//return 6 * node.phi + 2.1;
	return 0;
}

double u1(Node& node, int formulaNumber) {

	return node.t * node.t * node.r;

}

// Ввод данных
void Input(Grid& grid) {
	ifstream fInNodes("nodes.txt");
	ifstream fInFinitElems("finit_elems.txt");
	ifstream fInBoundConds1("boundary_conds1.txt");
	ifstream fInBoundConds2("boundary_conds2.txt");
	ifstream fInBoundConds3("boundary_conds3.txt");

	int count, temp;

	fInNodes >> count;

	grid.nodes.resize(count);

	for (int i = 0; i < count; i++) {
		fInNodes >> grid.nodes[i].r;
		fInNodes >> grid.nodes[i].phi;
		fInNodes >> grid.nodes[i].globalNumber;
	}

	fInNodes.close();

	fInFinitElems >> count;

	grid.finitElements.resize(count);

	for (int i = 0; i < count; i++) {
		grid.finitElements[i].nodes.resize(3);

		fInFinitElems >> temp;

		grid.finitElements[i].nodes[0] = grid.nodes[temp];

		fInFinitElems >> temp;

		grid.finitElements[i].nodes[1] = grid.nodes[temp];

		fInFinitElems >> temp;

		grid.finitElements[i].nodes[2] = grid.nodes[temp];

		fInFinitElems >> temp;

		grid.finitElements[i].formulaNumber = temp;
	}

	fInFinitElems.close();

	fInBoundConds1 >> count;

	grid.boundConds1.resize(count);

	for (int i = 0; i < count; i++) {
		fInBoundConds1 >> grid.boundConds1[i].globNum1;
		fInBoundConds1 >> grid.boundConds1[i].globNum2;
		fInBoundConds1 >> grid.boundConds1[i].formulaNumber;
	}

	fInBoundConds1.close();

	fInBoundConds2 >> count;

	grid.boundConds2.resize(count);

	for (int i = 0; i < count; i++) {
		fInBoundConds2 >> grid.boundConds2[i].globNum1;
		fInBoundConds2 >> grid.boundConds2[i].globNum2;
		fInBoundConds2 >> grid.boundConds2[i].formulaNumber;
	}

	fInBoundConds2.close();

	fInBoundConds3 >> count;

	grid.boundConds3.resize(count);

	for (int i = 0; i < count; i++) {
		fInBoundConds3 >> grid.boundConds3[i].globNum1;
		fInBoundConds3 >> grid.boundConds3[i].globNum2;
		fInBoundConds3 >> grid.boundConds3[i].formulaNumber;
	}

	fInBoundConds3.close();

}

void Portrait(Grid& grid, CRSMatrix& crsMatrix) {

	int elemAmount = grid.finitElements.size();
	int funcAmount = grid.nodes.size();

	crsMatrix.d.resize(funcAmount);
	crsMatrix.r.resize(funcAmount);
	crsMatrix.F.resize(funcAmount);
	crsMatrix.z.resize(funcAmount);
	crsMatrix.x.resize(funcAmount);
	crsMatrix.p.resize(funcAmount);
	crsMatrix.di.resize(funcAmount);
	crsMatrix.temp.resize(funcAmount);
	crsMatrix.temp0.resize(funcAmount);
	crsMatrix.n = funcAmount;
	crsMatrix.maxiter = 1E9;
	crsMatrix.epsilon = 1E-15;

	vector<int> listbeg;
	vector<int> list1;
	vector<int> list2;

	listbeg.resize(funcAmount);
	list1.resize(8 * funcAmount);
	list2.resize(8 * funcAmount);
	int listsize = -1;
	int iaddr;

	for (int ielem = 0; ielem < elemAmount; ielem++) {
		for (int i = 0; i < 3; i++) {
			int k = grid.finitElements[ielem].nodes[i].globalNumber;

			for (int j = i + 1; j < 3; j++) {
				int ind1 = k;
				int ind2 = grid.finitElements[ielem].nodes[j].globalNumber;

				if (ind2 < ind1) {
					ind1 = ind2;
					ind2 = k;
				}

				iaddr = listbeg[ind2];

				if (iaddr == 0) {
					listsize++;
					listbeg[ind2] = listsize;
					list1[listsize] = ind1;
					list2[listsize] = 0;
				}
				else {
					while (list1[iaddr] < ind1 && list2[iaddr] > 0) {
						iaddr = list2[iaddr];
					}

					if (list1[iaddr] > ind1) {
						listsize++;
						list1[listsize] = list1[iaddr];
						list2[listsize] = list2[iaddr];
						list1[iaddr] = ind1;
						list2[iaddr] = listsize;
					}
					else {
						if (list1[iaddr] < ind1) {
							listsize++;
							list2[iaddr] = listsize;
							list1[listsize] = ind1;
							list2[listsize] = 0;
						}
					}
				}

			}

		}
	}

	crsMatrix.ig.resize(funcAmount + 1);

	crsMatrix.ggl.resize(list1.size());
	crsMatrix.ggu.resize(list1.size());
	crsMatrix.u.resize(list1.size());
	crsMatrix.l.resize(list1.size());
	crsMatrix.jg.resize(list1.size());

	crsMatrix.ig[0] = 0;

	for (int i = 0; i < funcAmount; i++) {
		crsMatrix.ig[i + 1] = crsMatrix.ig[i];

		iaddr = listbeg[i];

		while (iaddr != 0) {
			crsMatrix.jg[crsMatrix.ig[i + 1] + 1] = list1[iaddr];
			crsMatrix.ig[i + 1]++;
			iaddr = list2[iaddr];
		}

	}

	for (int i = 2; i < funcAmount + 1; i++) {
		crsMatrix.ig[i]++;
	}

}

// LU Факторизация
void CalcLU(CRSMatrix& matrix) {
	double sumU, sumL, sumD;
	int n = matrix.n;


	for (int i = 0; i < n; i++) {

		sumD = 0;

		int begI = matrix.ig[i];
		int endI = matrix.ig[i + 1];
		for (int igi = begI; igi < endI; igi++) {

			sumU = 0;
			sumL = 0;

			int Jindex = matrix.jg[igi];

			for (int igj = begI; igj < igi; igj++) {

				int begJ = matrix.ig[Jindex];
				int endJ = matrix.ig[Jindex + 1];

				for (int jgi = begJ; jgi < endJ; jgi++) {

					if (matrix.jg[igj] == matrix.jg[jgi]) {

						sumL += matrix.l[igj] * matrix.u[jgi];

						sumU += matrix.l[jgi] * matrix.u[igj];

					}

				}
			}
			matrix.l[igi] -= sumL;
			matrix.u[igi] -= sumU;
			matrix.u[igi] /= matrix.d[Jindex];
			sumD += matrix.l[igi] * matrix.u[igi];
		}

		matrix.d[i] -= sumD;
	}
}

// Прямой ход Ly = F
void CalcDir(CRSMatrix& matrix, vector<double>& y, vector<double>& F) {
	double sum, buf;
	int n = matrix.n;


	for (int i = 0; i < n; i++) {
		y[i] = F[i];
	}

	for (int i = 0; i < n; i++) {

		sum = 0;

		int begI = matrix.ig[i];
		int endI = matrix.ig[i + 1];

		for (int igi = begI; igi < endI; igi++) {

			sum += y[matrix.jg[igi]] * matrix.l[igi];

		}

		buf = y[i] - sum;

		y[i] = buf / matrix.d[i];

	}

}

// Обратный ход Ux = y
void CalcRev(CRSMatrix& matrix, vector<double>& x, vector<double>& y) {
	int n = matrix.n;

	for (int i = 0; i < n; i++) {
		x[i] = y[i];
	}

	for (int i = matrix.n - 1; i >= 0; i--) {

		int begI = matrix.ig[i];
		int endI = matrix.ig[i + 1];

		for (int igi = begI; igi < endI; igi++) {

			x[matrix.jg[igi]] -= x[i] * matrix.u[igi];

		}

	}

}

// Умножение матрицы на вектор Ax = res
void MultMV(CRSMatrix& matrix, vector<double>& x, vector<double>& res) {
	int n = matrix.n;

	for (int i = 0; i < n; i++) {

		res[i] = matrix.di[i] * x[i];

		int begI = matrix.ig[i];
		int endI = matrix.ig[i + 1];

		for (int igi = begI; igi < endI; igi++) {

			int Jindex = matrix.jg[igi];

			res[i] += matrix.ggl[igi] * x[Jindex];
			res[Jindex] += matrix.ggu[igi] * x[i];

		}
	}
}

// Скалярное произведение двух векторов
double ScalarProd(vector<double>& x, vector<double>& y) {
	int n = x.size();

	double result = 0;

	for (int i = 0; i < n; i++) {
		result += x[i] * y[i];
	}

	return result;
}

// Локально-оптимальная схема c факторизацией LU
void LOS_LU(CRSMatrix& matrix) {

	double alpha, beta, norm, temp_nev = 0;

	int n = matrix.n, maxiter = matrix.maxiter;
	double epsilon = matrix.epsilon;

	CalcLU(matrix);
	// A * x0
	MultMV(matrix, matrix.x, matrix.temp);

	// f - A * x0
	for (int i = 0; i < n; i++) {
		matrix.temp[i] = matrix.F[i] - matrix.temp[i];
	}

	// L * r0 = f - A * x0
	CalcDir(matrix, matrix.r, matrix.temp);

	// U * z0 = r0
	CalcRev(matrix, matrix.z, matrix.r);

	// A * z0
	MultMV(matrix, matrix.z, matrix.temp);

	// L * p0 = A * z0
	CalcDir(matrix, matrix.p, matrix.temp);

	norm = ScalarProd(matrix.r, matrix.r);

	int k;

	for (k = 0; k < maxiter && norm > epsilon && norm != temp_nev; k++) {

		// если невязка не изменилась, то выходим из итерационного процесса
		temp_nev = norm;

		alpha = ScalarProd(matrix.p, matrix.r) / ScalarProd(matrix.p, matrix.p);

		for (int i = 0; i < n; i++) {

			matrix.x[i] = matrix.x[i] + alpha * matrix.z[i];
			matrix.r[i] = matrix.r[i] - alpha * matrix.p[i];

		}

		// U * temp = r
		CalcRev(matrix, matrix.temp, matrix.r);

		// A * U-1 * r = temp0
		MultMV(matrix, matrix.temp, matrix.temp0);

		// L * temp = A * U-1 * r 
		CalcDir(matrix, matrix.temp, matrix.temp0);

		beta = -1 * ScalarProd(matrix.p, matrix.temp) / ScalarProd(matrix.p, matrix.p);

		// U * temp0 = r
		CalcRev(matrix, matrix.temp0, matrix.r);

		norm = norm - alpha * alpha * ScalarProd(matrix.p, matrix.p);

		for (int i = 0; i < n; i++) {

			matrix.z[i] = matrix.temp0[i] + beta * matrix.z[i];
			matrix.p[i] = matrix.temp[i] + beta * matrix.p[i];

		}

	}

}

double det(FinitElement& finitElement) {

	return (finitElement.nodes[1].r - finitElement.nodes[0].r) * (finitElement.nodes[2].phi - finitElement.nodes[0].phi) -
		(finitElement.nodes[2].r - finitElement.nodes[0].r) * (finitElement.nodes[1].phi - finitElement.nodes[0].phi);

}

// Построение локальной матрицы массы и правой части
void MMatrix(FinitElement& finitElement, vector<vector<double>>& M, vector<double>& b) {

	double detD = det(finitElement);

	double coefB = abs(detD) / 120, sumR;

	double coef = coefB;

	vector<double> Fvalue;

	Fvalue.resize(3);

	for (int i = 0; i < 3; i++) {
		Fvalue[i] = F(finitElement.nodes[i], finitElement.formulaNumber);
		b[i] = 0;
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {

			sumR = (i == j) ? 6 * finitElement.nodes[i].r + 2 * finitElement.nodes[(i + 1) % 3].r + 2 * finitElement.nodes[(i + 2) % 3].r :
				2 * finitElement.nodes[i].r + 2 * finitElement.nodes[j].r + finitElement.nodes[3 - i - j].r;

			M[i][j] = sumR * coef;

			b[i] += Fvalue[j] * sumR * coefB;
		}
	}

}

void CalcGlobalB(Grid& grid, CRSMatrix& crsMatrix, DenseMatrix& denseMatrix) {
	for (int i = 0; i < grid.finitElements.size(); i++) {

		MMatrix(grid.finitElements[i], denseMatrix.M, denseMatrix.b);

		for (int k = 0; k < 3; k++) {
			crsMatrix.F[grid.finitElements[i].nodes[k].globalNumber] += denseMatrix.b[k];
		}
	}
}

struct TimeLayer {
	double deltaT, deltaT1, deltaT0;
	int layerCount;
	vector<vector<double>> q;
	vector<double> b, fLocal;
	vector<vector<double>> f;
	vector<double> t;

	TimeLayer(Grid& grid, CRSMatrix& matrix, DenseMatrix& denseMatrix, double startTime, double stopTime, int _layerCount) {

		int size = grid.nodes.size();
		q.resize(3);
		b.resize(3);
		f.resize(3);
		fLocal.resize(3);

		layerCount = _layerCount;

		t.resize(_layerCount);

		deltaT = (stopTime - startTime) / layerCount;

		for (int i = 0; i < layerCount; i++) {
			t[i] = deltaT * (i + 1);
		}

		for (int i = 0; i < 2; i++) {
			q[i].resize(size);
			f[i].resize(size);

			for (int j = 0; j < size; j++) {
				grid.nodes[j].t = t[i];
				q[i][j] = u1(grid.nodes[j], 0);
			}

			for (int j = 0; j < grid.finitElements.size(); j++) {
				for (int k = 0; k < 3; k++) {
					grid.finitElements[j].nodes[k].t = t[i];
				}
			}

			CalcGlobalB(grid, matrix, denseMatrix);

			f[i] = matrix.F;

			matrix.ggl.assign(matrix.ggl.size(), 0);
			matrix.ggu.assign(matrix.ggl.size(), 0);
			matrix.x.assign(matrix.di.size(), 0);
			matrix.temp.assign(matrix.di.size(), 0);
			matrix.temp0.assign(matrix.di.size(), 0);
			matrix.r.assign(matrix.di.size(), 0);
			matrix.p.assign(matrix.di.size(), 0);
			matrix.z.assign(matrix.di.size(), 0);
			matrix.F.assign(matrix.di.size(), 0);
			matrix.di.assign(matrix.di.size(), 0);

		}
		q[2].resize(size);
		f[2].resize(size);
	}
	TimeLayer() = default;
};

// Построение локальной матрицы жесткости 
void GMatrix(FinitElement& finitElement, vector<vector<double>>& G) {

	double detD = det(finitElement);

	double coef1 = (2 * lambda(finitElement.nodes[0], finitElement.formulaNumber) * finitElement.nodes[0].r
						+ 2 * lambda(finitElement.nodes[1], finitElement.formulaNumber) * finitElement.nodes[1].r
						+ 2 * lambda(finitElement.nodes[2], finitElement.formulaNumber) * finitElement.nodes[2].r

						+ lambda(finitElement.nodes[0], finitElement.formulaNumber) * finitElement.nodes[1].r
						+ lambda(finitElement.nodes[0], finitElement.formulaNumber) * finitElement.nodes[2].r

						+ lambda(finitElement.nodes[1], finitElement.formulaNumber) * finitElement.nodes[0].r
						+ lambda(finitElement.nodes[1], finitElement.formulaNumber) * finitElement.nodes[2].r

						+ lambda(finitElement.nodes[2], finitElement.formulaNumber) * finitElement.nodes[0].r
						+ lambda(finitElement.nodes[2], finitElement.formulaNumber) * finitElement.nodes[1].r) * abs(detD) / 24;

	double coef2 = (2 * lambda(finitElement.nodes[0], finitElement.formulaNumber) / finitElement.nodes[0].r
						+ 2 * lambda(finitElement.nodes[1], finitElement.formulaNumber) / finitElement.nodes[1].r
						+ 2 * lambda(finitElement.nodes[2], finitElement.formulaNumber) / finitElement.nodes[2].r

						+ lambda(finitElement.nodes[0], finitElement.formulaNumber) / finitElement.nodes[1].r
						+ lambda(finitElement.nodes[0], finitElement.formulaNumber) / finitElement.nodes[2].r

						+ lambda(finitElement.nodes[1], finitElement.formulaNumber) / finitElement.nodes[0].r
						+ lambda(finitElement.nodes[1], finitElement.formulaNumber) / finitElement.nodes[2].r

						+ lambda(finitElement.nodes[2], finitElement.formulaNumber) / finitElement.nodes[0].r
						+ lambda(finitElement.nodes[2], finitElement.formulaNumber) / finitElement.nodes[1].r) * abs(detD) / 24;

	// Первая строка

	G[0][0] = (coef1 * (finitElement.nodes[1].phi - finitElement.nodes[2].phi) * (finitElement.nodes[1].phi - finitElement.nodes[2].phi) +
				 coef2 * (finitElement.nodes[2].r - finitElement.nodes[1].r) * (finitElement.nodes[2].r - finitElement.nodes[1].r)) / (detD * detD);


	G[0][1] = (coef1 * (finitElement.nodes[1].phi - finitElement.nodes[2].phi) * (finitElement.nodes[2].phi - finitElement.nodes[0].phi) +
					coef2 * (finitElement.nodes[2].r - finitElement.nodes[1].r) * (finitElement.nodes[0].r - finitElement.nodes[2].r)) / (detD * detD);


	G[0][2] = (coef1 * (finitElement.nodes[1].phi - finitElement.nodes[2].phi) * (finitElement.nodes[0].phi - finitElement.nodes[1].phi) +
					coef2 * (finitElement.nodes[2].r - finitElement.nodes[1].r) * (finitElement.nodes[1].r - finitElement.nodes[0].r)) / (detD * detD);

	// Вторая строка

	G[1][0] = (coef1 * (finitElement.nodes[2].phi - finitElement.nodes[0].phi) * (finitElement.nodes[1].phi - finitElement.nodes[2].phi) +
					coef2 * (finitElement.nodes[0].r - finitElement.nodes[2].r) * (finitElement.nodes[2].r - finitElement.nodes[1].r)) / (detD * detD);


	G[1][1] = (coef1 * (finitElement.nodes[2].phi - finitElement.nodes[0].phi) * (finitElement.nodes[2].phi - finitElement.nodes[0].phi) +
					coef2 * (finitElement.nodes[0].r - finitElement.nodes[2].r) * (finitElement.nodes[0].r - finitElement.nodes[2].r)) / (detD * detD);


	G[1][2] = (coef1 * (finitElement.nodes[2].phi - finitElement.nodes[0].phi) * (finitElement.nodes[0].phi - finitElement.nodes[1].phi) +
					coef2 * (finitElement.nodes[0].r - finitElement.nodes[2].r) * (finitElement.nodes[1].r - finitElement.nodes[0].r)) / (detD * detD);

	// Третья строка

	G[2][0] = (coef1 * (finitElement.nodes[0].phi - finitElement.nodes[1].phi) * (finitElement.nodes[1].phi - finitElement.nodes[2].phi) +
					coef2 * (finitElement.nodes[1].r - finitElement.nodes[0].r) * (finitElement.nodes[2].r - finitElement.nodes[1].r)) / (detD * detD);


	G[2][1] = (coef1 * (finitElement.nodes[0].phi - finitElement.nodes[1].phi) * (finitElement.nodes[2].phi - finitElement.nodes[0].phi) +
					coef2 * (finitElement.nodes[1].r - finitElement.nodes[0].r) * (finitElement.nodes[0].r - finitElement.nodes[2].r)) / (detD * detD);


	G[2][2] = (coef1 * (finitElement.nodes[0].phi - finitElement.nodes[1].phi) * (finitElement.nodes[0].phi - finitElement.nodes[1].phi) +
					coef2 * (finitElement.nodes[1].r - finitElement.nodes[0].r) * (finitElement.nodes[1].r - finitElement.nodes[0].r)) / (detD * detD);

}

// Построение локальных матриц
void LocalMatrix(FinitElement& finitElement, DenseMatrix& denseMatrix, TimeLayer& timeLayer) {

	GMatrix(finitElement, denseMatrix.G);
	MMatrix(finitElement, denseMatrix.M, denseMatrix.b);

	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			denseMatrix.A[j][k] = denseMatrix.M[j][k] + denseMatrix.G[j][k];
		}
	}
}

void CalcGlobalG(Grid& grid, CRSMatrix& G, DenseMatrix& denseMatrix) {
	int temp;

	for (int i = 0; i < grid.finitElements.size(); i++) {

		GMatrix(grid.finitElements[i], denseMatrix.G);

		for (int k = 0; k < 3; k++) {

			int begI = grid.finitElements[i].nodes[k].globalNumber;

			for (int j = k + 1; j < 3; j++) {

				int endI = grid.finitElements[i].nodes[j].globalNumber;

				if (begI < endI) {

					temp = G.ig[endI];
					while (G.jg[temp++] - begI);
					temp--;
					G.ggl[temp] += denseMatrix.G[k][j];
					G.ggu[temp] += denseMatrix.G[j][k];

				}

				else {

					temp = G.ig[begI];
					while (G.jg[temp++] - endI);
					temp--;
					G.ggl[temp] += denseMatrix.G[k][j];
					G.ggu[temp] += denseMatrix.G[j][k];

				}

			}

			G.di[begI] += denseMatrix.G[k][k];

		}
	}
}

void CalcGlobalM(Grid& grid, CRSMatrix& M, DenseMatrix& denseMatrix) {
	int temp;

	for (int i = 0; i < grid.finitElements.size(); i++) {

		MMatrix(grid.finitElements[i], denseMatrix.M, denseMatrix.b);

		for (int k = 0; k < 3; k++) {

			int begI = grid.finitElements[i].nodes[k].globalNumber;

			for (int j = k + 1; j < 3; j++) {

				int endI = grid.finitElements[i].nodes[j].globalNumber;

				if (begI < endI) {

					temp = M.ig[endI];
					while (M.jg[temp++] - begI);
					temp--;
					M.ggl[temp] += denseMatrix.M[k][j];
					M.ggu[temp] += denseMatrix.M[j][k];

				}

				else {

					temp = M.ig[begI];
					while (M.jg[temp++] - endI);
					temp--;
					M.ggl[temp] += denseMatrix.M[k][j];
					M.ggu[temp] += denseMatrix.M[j][k];

				}

			}

			M.di[begI] += denseMatrix.M[k][k];

		}
	}
}

void GlobalMatrix(Grid& grid, CRSMatrix& crsMatrix, DenseMatrix& denseMatrix, TimeLayer& timeLayer) {

	int temp;

	for (int i = 0; i < grid.finitElements.size(); i++) {

		LocalMatrix(grid.finitElements[i], denseMatrix, timeLayer);

		for (int k = 0; k < 3; k++) {

			int begI = grid.finitElements[i].nodes[k].globalNumber;

			for (int j = k + 1; j < 3; j++) {

				int endI = grid.finitElements[i].nodes[j].globalNumber;

				if (begI < endI) {

					temp = crsMatrix.ig[endI];
					while (crsMatrix.jg[temp++] - begI);
					temp--;
					crsMatrix.ggl[temp] += denseMatrix.A[k][j];
					crsMatrix.ggu[temp] += denseMatrix.A[j][k];

				}

				else {

					temp = crsMatrix.ig[begI];
					while (crsMatrix.jg[temp++] - endI);
					temp--;
					crsMatrix.ggl[temp] += denseMatrix.A[k][j];
					crsMatrix.ggu[temp] += denseMatrix.A[j][k];

				}

			}

			crsMatrix.di[begI] += denseMatrix.A[k][k];

		}

		for (int k = 0; k < 3; k++) {

			crsMatrix.F[grid.finitElements[i].nodes[k].globalNumber] += denseMatrix.b[k];

		}
	}

}

double mesG(Node& node1, Node& node2) {
	return sqrt((node1.r - node2.r) * (node1.r - node2.r) + (node1.phi - node2.phi) * (node1.phi - node2.phi));
}

void CalcboundCond2(Grid& grid, CRSMatrix& crsMatrix, DenseMatrix& denseMatrix) {

	double hm;
	for (int i = 0; i < grid.boundConds2.size(); i++) {

		denseMatrix.globNumVert[0] = grid.boundConds2[i].globNum1;
		denseMatrix.globNumVert[1] = grid.boundConds2[i].globNum2;

		hm = mesG(grid.nodes[denseMatrix.globNumVert[0]], grid.nodes[denseMatrix.globNumVert[1]]);

		crsMatrix.F[denseMatrix.globNumVert[0]] += hm * (3 * tetta(grid.nodes[denseMatrix.globNumVert[0]], grid.boundConds2[i].formulaNumber) * grid.nodes[denseMatrix.globNumVert[0]].r
			+ tetta(grid.nodes[denseMatrix.globNumVert[0]], grid.boundConds2[i].formulaNumber) * grid.nodes[denseMatrix.globNumVert[1]].r
			+ tetta(grid.nodes[denseMatrix.globNumVert[1]], grid.boundConds2[i].formulaNumber) * grid.nodes[denseMatrix.globNumVert[0]].r
			+ tetta(grid.nodes[denseMatrix.globNumVert[1]], grid.boundConds2[i].formulaNumber) * grid.nodes[denseMatrix.globNumVert[1]].r) / 12;




		crsMatrix.F[denseMatrix.globNumVert[1]] += hm * (tetta(grid.nodes[denseMatrix.globNumVert[0]], grid.boundConds2[i].formulaNumber) * grid.nodes[denseMatrix.globNumVert[0]].r
			+ tetta(grid.nodes[denseMatrix.globNumVert[0]], grid.boundConds2[i].formulaNumber) * grid.nodes[denseMatrix.globNumVert[1]].r
			+ tetta(grid.nodes[denseMatrix.globNumVert[1]], grid.boundConds2[i].formulaNumber) * grid.nodes[denseMatrix.globNumVert[0]].r
			+ 3 * tetta(grid.nodes[denseMatrix.globNumVert[1]], grid.boundConds2[i].formulaNumber) * grid.nodes[denseMatrix.globNumVert[1]].r) / 12;
	}

}

void CalcboundCond3(Grid& grid, CRSMatrix& crsMatrix, DenseMatrix& denseMatrix) {

	double hm;

	int temp; double coef;

	for (int i = 0; i < grid.boundConds3.size(); i++) {

		denseMatrix.globNumVert[0] = grid.boundConds3[i].globNum1;
		denseMatrix.globNumVert[1] = grid.boundConds3[i].globNum2;

		hm = mesG(grid.nodes[denseMatrix.globNumVert[0]], grid.nodes[denseMatrix.globNumVert[1]]);

		coef = betaF() * hm / 12;

		denseMatrix.Acond3[0][0] = coef * (3 * grid.nodes[denseMatrix.globNumVert[0]].r + grid.nodes[denseMatrix.globNumVert[1]].r);
		denseMatrix.Acond3[0][1] = coef * (grid.nodes[denseMatrix.globNumVert[0]].r + grid.nodes[denseMatrix.globNumVert[1]].r);
		denseMatrix.Acond3[1][0] = coef * (grid.nodes[denseMatrix.globNumVert[0]].r + grid.nodes[denseMatrix.globNumVert[1]].r);
		denseMatrix.Acond3[1][1] = coef * (grid.nodes[denseMatrix.globNumVert[0]].r + 3 * grid.nodes[denseMatrix.globNumVert[1]].r);

		for (int k = 0; k < 2; k++) {

			int begI = denseMatrix.globNumVert[k];

			for (int j = k + 1; j < 2; j++) {

				int endI = denseMatrix.globNumVert[j];

				if (begI < endI) {

					temp = crsMatrix.ig[endI];
					while (crsMatrix.jg[temp++] - begI);
					temp--;
					crsMatrix.ggl[temp] += denseMatrix.Acond3[k][j];
					crsMatrix.ggu[temp] += denseMatrix.Acond3[j][k];

				}

				else {

					temp = crsMatrix.ig[begI];
					while (crsMatrix.jg[temp++] - endI);
					temp--;
					crsMatrix.ggl[temp] += denseMatrix.Acond3[k][j];
					crsMatrix.ggu[temp] += denseMatrix.Acond3[j][k];

				}

			}

			crsMatrix.di[begI] += denseMatrix.Acond3[k][k];

		}

		crsMatrix.F[denseMatrix.globNumVert[0]] += coef * (3 * uBeta(grid.nodes[denseMatrix.globNumVert[0]]) * grid.nodes[denseMatrix.globNumVert[0]].r
			+ uBeta(grid.nodes[denseMatrix.globNumVert[0]]) * grid.nodes[denseMatrix.globNumVert[1]].r
			+ uBeta(grid.nodes[denseMatrix.globNumVert[1]]) * grid.nodes[denseMatrix.globNumVert[0]].r
			+ uBeta(grid.nodes[denseMatrix.globNumVert[1]]) * grid.nodes[denseMatrix.globNumVert[1]].r);


		crsMatrix.F[denseMatrix.globNumVert[1]] += coef * (uBeta(grid.nodes[denseMatrix.globNumVert[0]]) * grid.nodes[denseMatrix.globNumVert[0]].r
			+ uBeta(grid.nodes[denseMatrix.globNumVert[0]]) * grid.nodes[denseMatrix.globNumVert[1]].r
			+ uBeta(grid.nodes[denseMatrix.globNumVert[1]]) * grid.nodes[denseMatrix.globNumVert[0]].r
			+ 3 * uBeta(grid.nodes[denseMatrix.globNumVert[1]]) * grid.nodes[denseMatrix.globNumVert[1]].r);

	}

}

void CalcboundCond1(Grid& grid, CRSMatrix& crsMatrix) {

	int temp, begI, endI;

	for (int i = 0; i < grid.boundConds1.size(); i++) {

		crsMatrix.di[grid.boundConds1[i].globNum1] = 1;
		crsMatrix.di[grid.boundConds1[i].globNum2] = 1;

		crsMatrix.F[grid.boundConds1[i].globNum1] = u1(grid.nodes[grid.boundConds1[i].globNum1], grid.boundConds1[i].formulaNumber);
		crsMatrix.F[grid.boundConds1[i].globNum2] = u1(grid.nodes[grid.boundConds1[i].globNum2], grid.boundConds1[i].formulaNumber);

		temp = crsMatrix.ig[grid.boundConds1[i].globNum1 + 1] - crsMatrix.ig[grid.boundConds1[i].globNum1];

		for (int j = 0; j < temp; j++) {
			crsMatrix.ggl[crsMatrix.ig[grid.boundConds1[i].globNum1] + j] = 0;
		}

		temp = crsMatrix.ig[grid.boundConds1[i].globNum2 + 1] - crsMatrix.ig[grid.boundConds1[i].globNum2];

		for (int j = 0; j < temp; j++) {
			crsMatrix.ggl[crsMatrix.ig[grid.boundConds1[i].globNum2] + j] = 0;
		}

		for (int j = grid.boundConds1[i].globNum1 + 1; j < grid.nodes.size(); j++) {

			begI = crsMatrix.ig[j];
			endI = crsMatrix.ig[j + 1];

			for (int k = begI; k < endI; k++) {

				if (crsMatrix.jg[k] == grid.boundConds1[i].globNum1) {
					crsMatrix.ggu[k] = 0;
					continue;
				}

			}

		}

		for (int j = grid.boundConds1[i].globNum2 + 1; j < grid.nodes.size(); j++)
		{

			begI = crsMatrix.ig[j];
			endI = crsMatrix.ig[j + 1];

			for (int k = begI; k < endI; k++) {

				if (crsMatrix.jg[k] == grid.boundConds1[i].globNum2) {
					crsMatrix.ggu[k] = 0;
					continue;
				}

			}

		}

	}
}

void Output(CRSMatrix& crsmatrix) {
	ofstream fOut("output.txt");

	for (int i = 0; i < crsmatrix.x.size(); i++) {
		fOut << fixed << scientific << setprecision(6) << crsmatrix.x[i] << endl;
	}
		
	fOut.close();
}

void procLU(CRSMatrix& crsMatrix) {

	for (int i = 0; i < crsMatrix.ggl.size(); i++) {
		crsMatrix.u[i] = crsMatrix.ggu[i];
		crsMatrix.l[i] = crsMatrix.ggl[i];
	}

	for (int i = 0; i < crsMatrix.n; i++) {
		crsMatrix.d[i] = crsMatrix.di[i];
	}

}

double getValue(FinitElement& finitElement, CRSMatrix& crsMatrix, double r, double phi) {
 
	double detD, s01, s12, s20, result = 0;
	vector<double> L(3);

	detD = det(finitElement);

	s01 = abs((finitElement.nodes[1].r - finitElement.nodes[0].r) * (phi - finitElement.nodes[0].phi) -
			(r - finitElement.nodes[0].r) * (finitElement.nodes[1].phi - finitElement.nodes[0].phi));

	s12 = abs((finitElement.nodes[2].r - finitElement.nodes[1].r) * (phi - finitElement.nodes[1].phi) -
			(r - finitElement.nodes[1].r) * (finitElement.nodes[2].phi - finitElement.nodes[1].phi));

	s20 = abs((finitElement.nodes[0].r - finitElement.nodes[2].r) * (phi - finitElement.nodes[2].phi) -
			(r - finitElement.nodes[2].r) * (finitElement.nodes[0].phi - finitElement.nodes[2].phi));

	L[0] = s12 / abs(detD);
	L[1] = s20 / abs(detD);
	L[2] = s01 / abs(detD);

	for (int j = 0; j < 3; j++) {
			result += crsMatrix.x[finitElement.nodes[j].globalNumber] * L[j];
	}
	
	return result;
}

void CalcA(CRSMatrix& crsMatrix, CRSMatrix& Mx, CRSMatrix& MSigma, CRSMatrix& G, TimeLayer& timeLayer) {

	for (int i = 0; i < crsMatrix.d.size(); i++) {
		crsMatrix.di[i] = Mx.di[i] / (timeLayer.deltaT0 * timeLayer.deltaT0) + MSigma.di[i] / (timeLayer.deltaT0 * 2) + G.di[i] / 2;
	}

	for (int i = 0; i < crsMatrix.ggl.size(); i++) {
		crsMatrix.ggl[i] = Mx.ggl[i] / (timeLayer.deltaT0 * timeLayer.deltaT0) + MSigma.ggl[i] / (timeLayer.deltaT0 * 2) + G.ggl[i] / 2;
		crsMatrix.ggu[i] = Mx.ggu[i] / (timeLayer.deltaT0 * timeLayer.deltaT0) + MSigma.ggu[i] / (timeLayer.deltaT0 * 2) + G.ggu[i] / 2;
	}

}

void CalcD(CRSMatrix& crsMatrix, CRSMatrix& Mx, CRSMatrix& MSigma, CRSMatrix& G, TimeLayer& timeLayer) {

	double coefXi = xi(0);
	double coefSigma = sigma(0);
	int n = crsMatrix.n;

	vector<double> MXqJ_1, MXqJ_2, MSqJ_1, MSqJ_2, GqJ_2;

	MXqJ_1.resize(n);
	MXqJ_2.resize(n);
	MSqJ_1.resize(n);
	MSqJ_2.resize(n);
	GqJ_2.resize(n);

	MultMV(Mx, timeLayer.q[1], MXqJ_1);
	MultMV(Mx, timeLayer.q[0], MXqJ_2);
	MultMV(MSigma, timeLayer.q[1], MSqJ_1);
	MultMV(MSigma, timeLayer.q[0], MSqJ_2);
	MultMV(G, timeLayer.q[0], GqJ_2);
	
	for (int i = 0; i < crsMatrix.F.size(); i++) {
		crsMatrix.F[i] = crsMatrix.F[i] / 2 
							+ timeLayer.f[0][i] / 2 
							+ MXqJ_1[i] * 2 / (timeLayer.deltaT0 * timeLayer.deltaT0) 
							- MXqJ_2[i] * 1 / (timeLayer.deltaT0 * timeLayer.deltaT0)
							+ MSqJ_2[i] / (timeLayer.deltaT0 * 2)
							//- MSqJ_1[i] * (timeLayer.deltaT0) / (timeLayer.deltaT1 * timeLayer.deltaT0)
							- GqJ_2[i] / 2;
	}

}

void SimpleIteration(Grid& grid, CRSMatrix& crsMatrix, CRSMatrix& M, CRSMatrix& Mx, CRSMatrix& MSigma, CRSMatrix& G, DenseMatrix& denseMatrix, TimeLayer& timeLayer) {

	double coefXi = xi(0);
	double coefSigma = sigma(0);

	CalcGlobalG(grid, G, denseMatrix);
	CalcGlobalM(grid, M, denseMatrix);

	for (int i = 0; i < crsMatrix.d.size(); i++) {
		Mx.di[i] = coefXi * M.di[i];
		MSigma.di[i] = coefSigma * M.di[i];
	}

	for (int i = 0; i < crsMatrix.ggl.size(); i++) {
		Mx.ggl[i] = coefXi * M.ggl[i];
		Mx.ggu[i] = coefXi * M.ggu[i];
		MSigma.ggl[i] = coefSigma * M.ggl[i];
		MSigma.ggu[i] = coefSigma * M.ggu[i];
	}

	for (int i = 2; i < timeLayer.layerCount; i++) {
		for (int j = 0; j < grid.nodes.size(); j++) {
			grid.nodes[j].t = timeLayer.t[i];
		}
		for (int j = 0; j < grid.finitElements.size(); j++) {
			for (int k = 0; k < 3; k++) {
				grid.finitElements[j].nodes[k].t = timeLayer.t[i];
			}
		}

		timeLayer.deltaT = timeLayer.t[i] - timeLayer.t[i - 2];
		timeLayer.deltaT1 = timeLayer.t[i - 1] - timeLayer.t[i - 2];
		timeLayer.deltaT0 = timeLayer.t[i] - timeLayer.t[i - 1];

		CalcA(crsMatrix, Mx, MSigma, G, timeLayer);
		CalcGlobalB(grid, crsMatrix, denseMatrix);

		timeLayer.f[2].assign(crsMatrix.F.begin(), crsMatrix.F.end());

 		CalcD(crsMatrix, Mx, MSigma, G, timeLayer);

		CalcboundCond1(grid, crsMatrix);

		procLU(crsMatrix);

		LOS_LU(crsMatrix);

		cout << i << "\t" << "временой слой" << endl;
		Output(crsMatrix);

		timeLayer.q[2].swap(crsMatrix.x);
 		timeLayer.q[0].swap(timeLayer.q[1]);
		timeLayer.q[1].swap(timeLayer.q[2]);
		
		timeLayer.f[0].swap(timeLayer.f[1]);
		timeLayer.f[1].swap(timeLayer.f[2]);

		crsMatrix.ggl.assign(crsMatrix.ggl.size(), 0);
		crsMatrix.ggu.assign(crsMatrix.ggl.size(), 0);
		crsMatrix.x.assign(crsMatrix.di.size(), 0);
		crsMatrix.temp.assign(crsMatrix.di.size(), 0);
		crsMatrix.temp0.assign(crsMatrix.di.size(), 0);
		crsMatrix.r.assign(crsMatrix.di.size(), 0);
		crsMatrix.p.assign(crsMatrix.di.size(), 0);
		crsMatrix.z.assign(crsMatrix.di.size(), 0);
		crsMatrix.F.assign(crsMatrix.di.size(), 0);
		crsMatrix.di.assign(crsMatrix.di.size(), 0);


	}

}

int main()
{
	Grid grid;
	DenseMatrix denseMatrix = DenseMatrix();
	CRSMatrix crsMatrix, layerMatrix, M, Mx, Msigma, G;
	TimeLayer timeLayer;
	
	Input(grid);

	Portrait(grid, crsMatrix);
	Portrait(grid, layerMatrix);
	Portrait(grid, Mx);
	Portrait(grid, M);
	Portrait(grid, Msigma);
	Portrait(grid, G);

	timeLayer = TimeLayer(grid, layerMatrix, denseMatrix, 0.0, 10.0, 10);

	SimpleIteration(grid, crsMatrix, M, Mx, Msigma, G, denseMatrix, timeLayer);
}