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

	double r, phi;
	double x, y;
	int globalNumber;

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
	
	return (formulaNumber == 0) ? -20 : 0;

}

double lambda(Node& node, int formulaNumber) {

	return (formulaNumber == 0) ? 10 : 1;

}

double gamma(int formulaNumber) {

	return 0;

}

double tetta(int formulaNumber) {
	return (formulaNumber == 0) ? 20 : 0;
}

double beta() {
	return 2;
}

double uBeta(Node& node) {
	return 20 * node.y - 27;
}

double u1(Node& node) {
	return node.y * node.y;
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
		fInNodes >> grid.nodes[i].x;
		fInNodes >> grid.nodes[i].y;
		grid.nodes[i].r = sqrt(grid.nodes[i].x * grid.nodes[i].x + grid.nodes[i].y * grid.nodes[i].y);
		grid.nodes[i].phi = acos(grid.nodes[i].x / grid.nodes[i].r) * 180.0 / M_PI;
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
	list1.resize(2 * (funcAmount - 1));
	list2.resize(2 * (funcAmount - 1));
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
void CalcDir(CRSMatrix& matrix, vector<double>&y, vector<double>& F) {
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

// Обратный ход Ux = phi
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

// Построение локальной матрицы жесткости 
void GMatrix(FinitElement& finitElement, vector<vector<double>>& G) {
	
	double detD = (finitElement.nodes[1].x - finitElement.nodes[0].x) * (finitElement.nodes[2].y - finitElement.nodes[0].y) -
		(finitElement.nodes[2].x - finitElement.nodes[0].x) * (finitElement.nodes[1].y - finitElement.nodes[0].y);

	double coef = (lambda(finitElement.nodes[0], finitElement.formulaNumber) * finitElement.nodes[0].r 
					+ lambda(finitElement.nodes[1], finitElement.formulaNumber) * finitElement.nodes[1].r
					+ lambda(finitElement.nodes[2], finitElement.formulaNumber) * finitElement.nodes[2].r

					+ lambda(finitElement.nodes[0], finitElement.formulaNumber) * finitElement.nodes[1].r / 2
					+ lambda(finitElement.nodes[0], finitElement.formulaNumber) * finitElement.nodes[2].r / 2

					+ lambda(finitElement.nodes[1], finitElement.formulaNumber) * finitElement.nodes[0].r / 2
					+ lambda(finitElement.nodes[1], finitElement.formulaNumber) * finitElement.nodes[2].r / 2

					+ lambda(finitElement.nodes[2], finitElement.formulaNumber) * finitElement.nodes[0].r / 2
					+ lambda(finitElement.nodes[2], finitElement.formulaNumber) * finitElement.nodes[1].r / 2 ) * abs(detD) / 12;

	// Первая строка

	G[0][0] = coef * ((finitElement.nodes[1].y - finitElement.nodes[2].y) * (finitElement.nodes[1].y - finitElement.nodes[2].y) +
		(finitElement.nodes[2].x - finitElement.nodes[1].x) * (finitElement.nodes[2].x - finitElement.nodes[1].x)) / (detD * detD);


	G[0][1] = coef * ((finitElement.nodes[1].y - finitElement.nodes[2].y) * (finitElement.nodes[2].y - finitElement.nodes[0].y) +
		(finitElement.nodes[2].x - finitElement.nodes[1].x) * (finitElement.nodes[0].x - finitElement.nodes[2].x)) / (detD * detD);


	G[0][2] = coef * ((finitElement.nodes[1].y - finitElement.nodes[2].y) * (finitElement.nodes[0].y - finitElement.nodes[1].y) +
		(finitElement.nodes[2].x - finitElement.nodes[1].x) * (finitElement.nodes[1].x - finitElement.nodes[0].x)) / (detD * detD);

	// Вторая строка

	G[1][0] = coef * ((finitElement.nodes[2].y - finitElement.nodes[0].y) * (finitElement.nodes[1].y - finitElement.nodes[2].y) +
		(finitElement.nodes[0].x - finitElement.nodes[2].x) * (finitElement.nodes[2].x - finitElement.nodes[1].x)) / (detD * detD);


	G[1][1] = coef * ((finitElement.nodes[2].y - finitElement.nodes[0].y) * (finitElement.nodes[2].y - finitElement.nodes[0].y) +
		(finitElement.nodes[0].x - finitElement.nodes[2].x) * (finitElement.nodes[0].x - finitElement.nodes[2].x)) / (detD * detD);


	G[1][2] = coef * ((finitElement.nodes[2].y - finitElement.nodes[0].y) * (finitElement.nodes[0].y - finitElement.nodes[1].y) +
		(finitElement.nodes[0].x - finitElement.nodes[2].x) * (finitElement.nodes[1].x - finitElement.nodes[0].x)) / (detD * detD);

	// Третья строка

	G[2][0] = coef * ((finitElement.nodes[0].y - finitElement.nodes[1].y) * (finitElement.nodes[1].y - finitElement.nodes[2].y) +
		(finitElement.nodes[1].x - finitElement.nodes[0].x) * (finitElement.nodes[2].x - finitElement.nodes[1].x)) / (detD * detD);


	G[2][1] = coef * ((finitElement.nodes[0].y - finitElement.nodes[1].y) * (finitElement.nodes[2].y - finitElement.nodes[0].y) +
		(finitElement.nodes[1].x - finitElement.nodes[0].x) * (finitElement.nodes[0].x - finitElement.nodes[2].x)) / (detD * detD);


	G[2][2] = coef * ((finitElement.nodes[0].y - finitElement.nodes[1].y) * (finitElement.nodes[0].y - finitElement.nodes[1].y) +
		(finitElement.nodes[1].x - finitElement.nodes[0].x) * (finitElement.nodes[1].x - finitElement.nodes[0].x)) / (detD * detD);
	
}

// Построение локальной матрицы массы и правой части
void MMatrix(FinitElement& finitElement, vector<vector<double>>& M, vector<double>& b) {

	double detD = (finitElement.nodes[1].x - finitElement.nodes[0].x) * (finitElement.nodes[2].y - finitElement.nodes[0].y) -
		(finitElement.nodes[2].x - finitElement.nodes[0].x) * (finitElement.nodes[1].y - finitElement.nodes[0].y);

	double coefB = abs(detD) / 120, sumR;

	double coef = gamma(finitElement.formulaNumber) * coefB;

	vector<double> Fvalue;

	Fvalue.resize(3);

	for (int i = 0; i < 3; i++) {
		Fvalue[i] = F(finitElement.nodes[i], finitElement.formulaNumber);
		b[i] = 0;
	}
	

	/*b[0] = coefB * (6 * Fvalue[0] * finitElement.nodes[0].r + 2 * Fvalue[0] * finitElement.nodes[1].r + 2 * Fvalue[0] * finitElement.nodes[2].r +
						 2 * Fvalue[1] * finitElement.nodes[0].r + 2 * Fvalue[1] * finitElement.nodes[1].r + Fvalue[1] * finitElement.nodes[2].r + 
						 2 * Fvalue[2] * finitElement.nodes[0].r + Fvalue[2] * finitElement.nodes[1].r + 2 * Fvalue[2] * finitElement.nodes[2].r);

	b[1] = coefB * (2 * Fvalue[0] * finitElement.nodes[0].r + 2 * Fvalue[0] * finitElement.nodes[1].r + Fvalue[0] * finitElement.nodes[2].r +
						 2 * Fvalue[1] * finitElement.nodes[0].r + 6 * Fvalue[1] * finitElement.nodes[1].r + 2 * Fvalue[1] * finitElement.nodes[2].r +
						 Fvalue[2] * finitElement.nodes[0].r + 2* Fvalue[2] * finitElement.nodes[1].r + 2 * Fvalue[2] * finitElement.nodes[2].r);

	b[2] = coefB * (2 * Fvalue[0] * finitElement.nodes[0].r + Fvalue[0] * finitElement.nodes[1].r + 2 * Fvalue[0] * finitElement.nodes[2].r +
						 Fvalue[1] * finitElement.nodes[0].r + 2 * Fvalue[1] * finitElement.nodes[1].r + 2 * Fvalue[1] * finitElement.nodes[2].r +
						 2 * Fvalue[2] * finitElement.nodes[0].r + 2 * Fvalue[2] * finitElement.nodes[1].r + 6 * Fvalue[2] * finitElement.nodes[2].r);*/

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			
			sumR = (i == j) ? 6 * finitElement.nodes[i].r + 2 * finitElement.nodes[(i + 1) % 3].r + 2 * finitElement.nodes[(i + 2) % 3].r :
									2 * finitElement.nodes[i].r + 2 * finitElement.nodes[j].r + finitElement.nodes[3 - i - j].r;
			
			M[i][j] = sumR * coef;

			b[i] += Fvalue[j] * sumR * coefB;
		}
	}

}
 
// Построение локальных матриц
void LocalMatrix(FinitElement& finitElement, DenseMatrix& denseMatrix) {

	GMatrix(finitElement, denseMatrix.G);
	MMatrix(finitElement, denseMatrix.M, denseMatrix.b);

	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			denseMatrix.A[j][k] = denseMatrix.G[j][k] + denseMatrix.M[j][k];
		}
	}
}

void GlobalMatrix(Grid& grid, CRSMatrix& crsMatrix, DenseMatrix& denseMatrix) {

	int temp;

	for (int i = 0; i < grid.finitElements.size(); i++) {

		LocalMatrix(grid.finitElements[i], denseMatrix);

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

double mesG(Node node1, Node node2) {
	return sqrt((node1.x - node2.x) * (node1.x - node2.x) + (node1.y - node2.y) * (node1.y - node2.y));
}

void CalcboundCond2(Grid& grid, CRSMatrix& crsMatrix, DenseMatrix& denseMatrix) {

	double hm;
	for (int i = 0; i < grid.boundConds2.size(); i++) {
	
		denseMatrix.globNumVert[0] = grid.boundConds2[i].globNum1;
		denseMatrix.globNumVert[1] = grid.boundConds2[i].globNum2;

		hm = mesG(grid.nodes[denseMatrix.globNumVert[0]], grid.nodes[denseMatrix.globNumVert[1]]);
	
		crsMatrix.F[denseMatrix.globNumVert[0]] += hm * (2 * tetta(grid.boundConds2[i].formulaNumber) + tetta(grid.boundConds2[i].formulaNumber)) / 6;
		crsMatrix.F[denseMatrix.globNumVert[1]] += hm * (tetta(grid.boundConds2[i].formulaNumber) + 2 * tetta(grid.boundConds2[i].formulaNumber)) / 6;
	}

}

void CalcboundCond3(Grid& grid, CRSMatrix& crsMatrix, DenseMatrix& denseMatrix) {

	double hm;
	//int globNumVert1, globNumVert2;
	int temp; double coef;

	for (int i = 0; i < grid.boundConds3.size(); i++) {

		denseMatrix.globNumVert[0] = grid.boundConds3[i].globNum1;
		denseMatrix.globNumVert[1] = grid.boundConds3[i].globNum2;

		hm = mesG(grid.nodes[denseMatrix.globNumVert[0]], grid.nodes[denseMatrix.globNumVert[1]]);

		coef = beta() * hm / 6;
	
		denseMatrix.Acond3[0][0] = 2 * coef;
		denseMatrix.Acond3[0][1] = coef;
		denseMatrix.Acond3[1][0] = coef;
		denseMatrix.Acond3[1][1] = 2 * coef;

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

		crsMatrix.F[denseMatrix.globNumVert[0]] += coef * (2 * uBeta(grid.nodes[denseMatrix.globNumVert[0]]) + uBeta(grid.nodes[denseMatrix.globNumVert[1]]));
		crsMatrix.F[denseMatrix.globNumVert[1]] += coef * (uBeta(grid.nodes[denseMatrix.globNumVert[0]]) + 2 * uBeta(grid.nodes[denseMatrix.globNumVert[1]]));

	}

}

void CalcboundCond1(Grid& grid, CRSMatrix& crsMatrix, DenseMatrix& denseMatrix) {

	int temp, begI, endI;

	for (int i = 0; i < grid.boundConds1.size(); i++) {

		crsMatrix.di[grid.boundConds1[i].globNum1] = 1;
		crsMatrix.di[grid.boundConds1[i].globNum2] = 1;

		crsMatrix.F[grid.boundConds1[i].globNum1] = u1(grid.nodes[grid.boundConds1[i].globNum1]);
		crsMatrix.F[grid.boundConds1[i].globNum2] = u1(grid.nodes[grid.boundConds1[i].globNum2]);

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

void Output(CRSMatrix crsmatrix) {
	ofstream fOut("output.txt");

	for (int i = 0; i < crsmatrix.r.size(); i++) {
		fOut << crsmatrix.x[i] << '\t';
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

int main()
{
	Grid grid;
	DenseMatrix denseMatrix = DenseMatrix();
	CRSMatrix crsMatrix;

	Input(grid);

	Portrait(grid, crsMatrix);

	GlobalMatrix(grid, crsMatrix, denseMatrix);

	CalcboundCond2(grid, crsMatrix, denseMatrix);

	CalcboundCond3(grid, crsMatrix, denseMatrix);
	
	CalcboundCond1(grid, crsMatrix, denseMatrix);
	
	procLU(crsMatrix);

	LOS_LU(crsMatrix);

	Output(crsMatrix);
}

