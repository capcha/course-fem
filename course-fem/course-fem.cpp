#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

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

};

struct Grid {

	vector<Node> nodes;	
	vector<FinitElement> finitElements;
	vector<double> x, f;
	vector<BoundCond1> boundConds1;
	vector<BoundCond2> boundConds2;
	vector<BoundCond3> boundConds3;

};

struct Node {

	int x, y;

};

struct FinitElement {

	vector<Node> nodes;
	int areaNumber;

};

struct BoundCond1 {

	int formulaNumber;

};

struct BoundCond2 {

	int vertex1, vertex2, formulaNumber;

};

struct BoundCond3 {

	int vertex1, vertex2, formulaNumber;

};

// Ввод данных
void input(Grid& grid) {
	ifstream fInNodes("nodes.txt");
	ifstream fInFinitElems("finit_elems.txt");
	ifstream fInBoundConds1("boudary_conds1.txt");
	ifstream fInBoundConds2("boudary_conds2.txt");
	ifstream fInBoundConds3("boudary_conds3.txt");

	int count, temp;

	fInNodes >> count;
	
	grid.nodes.resize(count);

	for (int i = 0; i < count; i++) {
		fInNodes >> grid.nodes[i].x;
		fInNodes >> grid.nodes[i].y;
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

		grid.finitElements[i].areaNumber = i;
	}

	fInFinitElems.close();

	fInBoundConds1 >> count;
	
	grid.boundConds1.resize(count);

	for (int i = 0; i < count; i++) {
		fInBoundConds1 >> grid.boundConds1[i].formulaNumber;
	}

	fInBoundConds1.close();

	fInBoundConds2 >> count;

	grid.boundConds2.resize(count);

	for (int i = 0; i < count; i++) {
		fInBoundConds2 >> grid.boundConds2[i].vertex1;
		fInBoundConds2 >> grid.boundConds2[i].vertex2;
		fInBoundConds2 >> grid.boundConds2[i].formulaNumber;
	}

	fInBoundConds2.close();

	fInBoundConds2 >> count;

	grid.boundConds3.resize(count);

	for (int i = 0; i < count; i++) {
		fInBoundConds3 >> grid.boundConds3[i].vertex1;
		fInBoundConds3 >> grid.boundConds3[i].vertex2;
		fInBoundConds3 >> grid.boundConds3[i].formulaNumber;
	}

	fInBoundConds3.close();

}

// LU Факторизация
void calcLU(CRSMatrix& matrix) {
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
void calcDir(CRSMatrix& matrix, vector<double>& y, vector<double>& F) {
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
void calcRev(CRSMatrix& matrix, vector<double>& x, vector<double>& y) {
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
void multMV(CRSMatrix& matrix, vector<double>& x, vector<double>& res) {
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
double scalarProd(vector<double>& x, vector<double>& y) {
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

	calcLU(matrix);
	// A * x0
	multMV(matrix, matrix.x, matrix.temp);

	// f - A * x0
	for (int i = 0; i < n; i++) {
		matrix.temp[i] = matrix.F[i] - matrix.temp[i];
	}

	// L * r0 = f - A * x0
	calcDir(matrix, matrix.r, matrix.temp);

	// U * z0 = r0
	calcRev(matrix, matrix.z, matrix.r);

	// A * z0
	multMV(matrix, matrix.z, matrix.temp);

	// L * p0 = A * z0
	calcDir(matrix, matrix.p, matrix.temp);

	norm = scalarProd(matrix.r, matrix.r);

	int k;

	for (k = 0; k < maxiter && norm > epsilon && norm != temp_nev; k++) {
		cout << "iter:\t" << k << "\tnev:\t" << setprecision(4) << norm << endl;

		// если невязка не изменилась, то выходим из итерационного процесса
		temp_nev = norm;

		alpha = scalarProd(matrix.p, matrix.r) / scalarProd(matrix.p, matrix.p);

		for (int i = 0; i < n; i++) {

			matrix.x[i] = matrix.x[i] + alpha * matrix.z[i];
			matrix.r[i] = matrix.r[i] - alpha * matrix.p[i];

		}

		// U * temp = r
		calcRev(matrix, matrix.temp, matrix.r);

		// A * U-1 * r = temp0
		multMV(matrix, matrix.temp, matrix.temp0);

		// L * temp = A * U-1 * r 
		calcDir(matrix, matrix.temp, matrix.temp0);

		beta = -1 * scalarProd(matrix.p, matrix.temp) / scalarProd(matrix.p, matrix.p);

		// U * temp0 = r
		calcRev(matrix, matrix.temp0, matrix.r);

		norm = norm - alpha * alpha * scalarProd(matrix.p, matrix.p);

		for (int i = 0; i < n; i++) {

			matrix.z[i] = matrix.temp0[i] + beta * matrix.z[i];
			matrix.p[i] = matrix.temp[i] + beta * matrix.p[i];

		}

	}

	cout << setprecision(4) << norm;
	cout << "\t" << k << endl;

}


double lambda(int areaNumber) {

	if (areaNumber == 0) {
		return 10;
	}

	return 0;

}

double gamma(int areaNumber) {

	if (areaNumber == 0) {
		return 2;
	}

	return 0;

}

// Построение локальной матрицы жесткости 
void GMatrix(FinitElement& finitElement, vector<vector<double>>& G) {
	
	double detD = (finitElement.nodes[1].x - finitElement.nodes[0].x) * (finitElement.nodes[2].y - finitElement.nodes[0].y) -
						(finitElement.nodes[2].x - finitElement.nodes[0].x) * (finitElement.nodes[1].y - finitElement.nodes[0].y);

	double coef = lambda(finitElement.areaNumber) * abs(detD) / 2;


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

// Построение локальной матрицы массы 
void MMatrix(FinitElement& finitElement, vector<vector<double>>& M) {

	double detD = (finitElement.nodes[1].x - finitElement.nodes[0].x) * (finitElement.nodes[2].y - finitElement.nodes[0].y) -
		(finitElement.nodes[2].x - finitElement.nodes[0].x) * (finitElement.nodes[1].y - finitElement.nodes[0].y);

	double coef = gamma(finitElement.areaNumber) * abs(detD) / 24;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			M[i][j] = (i == j) ? 2 * coef : coef;
		}
	}

}

void LocalMatrix(FinitElement& finitElement, DenseMatrix& denseMatrix, vector<vector<double>>& A) {

	GMatrix(finitElement, denseMatrix.G);
	MMatrix(finitElement, denseMatrix.M);

	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			denseMatrix.A[j][k] = denseMatrix.G[j][k] + denseMatrix.M[j][k];
		}
	}
}

int main()
{
	Grid grid;
	DenseMatrix denseMatrix;
	CRSMatrix crsMatrix;

	input(grid);
	
}

