#include "utils.hpp"
#include <vector>

// FCH settings
int FCHlen;      // 链的长度
double interAng; // 内转角大小

// 输出文件存储地址
char dirname[40];

// 球面剖分
const int pSphere = 2 * pAng * (pAng - 1) + 2;

// ODF, ODF导数, 暂存ODF
// f, f', g
double f[pSphere];
double g[pSphere];
double fDer[pAng + 1];

// Irep表, Iatt表
double Irep[pSphere];
double Iatt[pSphere];

// A表, V表, M表
double A[pSphere][pSphere];
double V[pSphere][pSphere];
double M[2][pSphere][pSphere];

double eta;

// --------------------角度换算相关--------------------
Vector3d Sphere[pSphere]; // 球上的坐标
int ind[pAng + 1];        // 球坐标转一维的数组
double calPhi(int N) {
	if (N == 0)
		return 0;
	int i = (N - 1) / (2 * pAng) + 1;
	return i * dAng;
}
double calTheta(int N) {
	if (N == 0 | N == pSphere - 1)
		return 0;
	int j = (N - 1) % (2 * pAng);
	return j * dAng;
}
void initCoord() {
	for (int i = 0; i < pSphere; i++) {
		double phi = calPhi(i);
		double theta = calTheta(i);
		Sphere[i] =
		    Vector3d(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi));
	}

	// 0, {1+2pAng0}, {1+2pAng1}, ..., 1+2pAng(pAng-1) 是不同的phi角的整体编号
	ind[0] = 0;
	for (int i = 0; i < pAng; i++)
		ind[i + 1] = 1 + 2 * pAng * i;
}
inline double angle(int p, int q) {
	double phiP = calPhi(p);
	double phiQ = calPhi(q);

	double thetaP = calTheta(p);
	double thetaQ = calTheta(q);

	Vector3d A(sin(phiP) * cos(thetaP), sin(phiP) * sin(thetaP), cos(phiP));
	Vector3d B(sin(phiQ) * cos(thetaQ), sin(phiQ) * sin(thetaQ), cos(phiQ));

	double dotProd = A.dot(B);
	// 防止数值精度不够溢出
	if (dotProd > 1 || dotProd < -1)
		dotProd = dotProd > 0 ? 1.0 : -1.0;

	double res = acos(dotProd);
	return res;
}

// --------------------文件操作模块--------------------
inline bool read1d(char name[], double dat[pSphere]) {
	FILE *file = fopen(name, "r");
	if (file == NULL)
		return false;

	for (int p = 0; p < pSphere; p++)
		fscanf(file, "%lf", &dat[p]);
	fclose(file);

	return true;
}
inline bool read2d(char name[], double dat[pSphere][pSphere]) {
	FILE *file = fopen(name, "r");
	if (file == NULL)
		return false;

	for (int p = 0; p < pSphere; p++)
		for (int q = 0; q < pSphere; q++)
			fscanf(file, "%lf", &dat[p][q]);
	fclose(file);

	return true;
}
inline void dump1d(char name[], double dat[pSphere]) {
	FILE *file = fopen(name, "w");
	for (int p = 0; p < pSphere; p++)
		fprintf(file, "%.10lf\n", dat[p]);
	fclose(file);
}
inline void dump2d(char name[], double dat[pSphere][pSphere]) {
	FILE *file = fopen(name, "w");
	for (int p = 0; p < pSphere; p++) {
		for (int q = 0; q < pSphere; q++)
			fprintf(file, "%.10lf ", dat[p][q]);
		fprintf(file, "\n");
	}
	fclose(file);
}
inline void outODF(char name[]) {
	FILE *file = fopen(name, "w");
	for (int i = 0; i <= pAng; i++)
		fprintf(file, "%lf\n", f[ind[i]]);
	fclose(file);
}

// --------------------特殊函数--------------------
// Onsager trial function
inline double onsagerTrial(double a, double x) {
	return a * cosh(a * cos(x)) / (4 * M_PI * sinh(a));
}
// Parsons-Lee/Onsager's G(\eta)
#define OS 0
#define PL 1
inline double G(double eta, int type = PL) {
	// Onsager's G
	if (type == 0) {
		return 4 * eta;
	}
	// PL's G
	return eta * (4 - 3 * eta) / pow(1 - eta, 2);
}

// --------------------打表用的队列优化模块--------------------
typedef pair<double, double> pDouble;
// 二分插入
void insertx(vector<pDouble> &que, pDouble x) {
	int st = 0;
	int ed = que.size() - 1;
	int md;

	if (ed < 0) {
		que.insert(que.begin(), x);
		return;
	}

	while (st + 1 < ed) {
		md = (st + ed) / 2;
		if (x.first < que[md].first)
			ed = md;
		else
			st = md;
	}
	if (st == ed) // 必为st == 0
		que.insert(que.end(), x);
	else {
		// st + 1 == ed
		if (x.first < que[ed].first)
			que.insert(que.begin() + ed, x);
		else
			que.insert(que.begin() + ed + 1, x);
	}
}
// 二分寻找
pair<bool, double> findx(vector<pDouble> &que, pDouble x) {
	int st = 0;
	int ed = que.size() - 1;
	int md;

	if (ed < 0)
		return make_pair(false, 0.);

	while (st + 1 < ed) {
		md = (st + ed) / 2;
		if (x.first < que[md].first)
			ed = md;
		else
			st = md;
	}

	if (st == ed) {
		// 必为st == 0
		if (fabs(x.first - que[st].first) < EPS)
			return make_pair(true, que[st].second);
	} else {
		// st + 1 == ed
		if (fabs(x.first - que[st].first) < EPS)
			return make_pair(true, que[st].second);
		if (fabs(x.first - que[ed].first) < EPS)
			return make_pair(true, que[ed].second);
	}
	return make_pair(false, 0.);
}
// 测试队列优化
void testVEC() {
	char fname[120];
	{
		sprintf(fname, "%sData.Vex.txt", dirname);
		read2d(fname, V);

		sprintf(fname, "%sData.500K20Aatt.txt", dirname);
		read2d(fname, A);

		sprintf(fname, "%sData.1000K10M1.txt", dirname);
		read2d(fname, M[0]);

		sprintf(fname, "%sData.1000K10M2.txt", dirname);
		read2d(fname, M[1]);
	}

	vector<pDouble> que[3];
	for (int k = 0; k < 3; k++) {
		for (int p = 0; p < pSphere; p++) {
			for (int q = 0; q < pSphere; q++) {
				double res;
				double ang = angle(p, q);
				pair<bool, double> status = findx(que[k], make_pair(ang, 0.));
				if (status.first)
					res = status.second;
				else {
					if (k < 2)
						res = M[k][p][q];
					else if (k == 2)
						res = A[p][q];
					insertx(que[k], make_pair(ang, res));
				}
			}
		}
		char name[100];
		FILE *file;

		if (k < 2)
			sprintf(name, "%sM%d.txt", dirname, k + 1);
		else if (k == 2)
			sprintf(name, "%sA.txt", dirname);

		file = fopen(name, "w");
		for (int i = 0; i < que[k].size(); i++)
			fprintf(file, "%lf %lf\n", que[k][i].first, que[k][i].second);
		fclose(file);
	}

	/*
	insertx(que, make_pair(1.11, 1));
	insertx(que, make_pair(0.11, 2));
	insertx(que, make_pair(1.112, 3));
	insertx(que, make_pair(2.11, 4));
	insertx(que, make_pair(1.11002, 8));
	insertx(que, make_pair(1.11001, 5));

	for (int i = 0; i < que.size(); i++)
	    printf("%lf, %lf\n", que[i].first, que[i].second);

	pair<bool, double> tmp = findx(que, make_pair(1.110027, 0));
	if (tmp.first)
	    printf("%lf\n", tmp.second);
	else
	    printf("NO\n");
	*/
}

// --------------------表计算模块--------------------
// excluded volume using fitting data, (phi(*), theta(*)) denotes the direction
// type = 0,1,2,3
// V, HSC
inline double RXVofHSC(int p, int q) {
	double x = angle(p, q);

	// RXV for spherocylinder
	double k[4] = {3.0, 3.2, 4.0, 5.0};
	int type = 0;
	// According to JCP_1995_Jackson
	double exc =
	    4.0 / 3 * M_PI + 2 * M_PI * k[type] + 2 * k[type] * k[type] * sin(x);
	double hsc8 = 4.0 / 3 * M_PI + 2 * M_PI * k[type];
	return exc / hsc8;
}
// V, (from Joachim)
inline double RXVofLCH(int p, int q) {
	double x = angle(p, q);
	/*
	    // RXV for LCH chain.
	    double a = 1.052;
	    double b = 3.664;
	    double c = -1.168;
	    double d = 0.238;
	    double e = -0.077;

	    // According to Simulation results
	    double res = (a + b * x + c * x * x) / (1 + d * x + e * x * x);
	*/
	int m = FCHlen;
	double res = 11 - 3. / m + 3.5339 * (m + 1. / m - 2) * sin(x);

	return res / 8.;
}
// V表
inline void VexTablet() {
	for (int p = 0; p < pSphere; p++)
		for (int q = 0; q < pSphere; q++)
			V[p][q] = RXVofLCH(p, q);

	char fname[120];
	sprintf(fname, "%sData.Vex.txt", dirname);
	dump2d(fname, V);
}
// M表, t=1或2
inline void MTablet(int step1, int step2, int t) {
	// 队列优化
	// vector<pDouble> que;

	char filename[120];
	sprintf(filename, "%sData.%dK%dM%d.txt", dirname, step1 / 1000, step2, t);
	FILE *f = fopen(filename, "w");

	// printf("here\n");

	for (int pi = 0; pi < pAng; pi++) {
		int p0 = ind[pi];

		for (int qi = 0; qi < pAng; qi++) {
			for (int q = ind[qi]; q < ind[qi + 1]; q++) {
				double res = 0.;
				res = mpiAvgForMk(step1, step2, t, FCHlen, interAng, Sphere[p0],
				                  Sphere[q], fch, false);

				M[t - 1][p0][q] = res;
			}
		}

		for (int p = p0 + 1; p < ind[pi + 1]; p++) {
			for (int qi = 0; qi < pAng; qi++) {
				for (int q = ind[qi]; q < ind[qi + 1]; q++) {
					int qshift = q - (p - p0);

					if (qshift < ind[qi])
						qshift += ind[qi + 1] - ind[qi];

					M[t - 1][p][q] = M[t - 1][p0][qshift];
				}
			}
		}
	}

	if (rk == 0) {
		for (int p = 0; p < pSphere; p++) {
			for (int q = 0; q < pSphere; q++)
				fprintf(f, "%lf ", M[t - 1][p][q]);
			fprintf(f, "\n");
		}
	}
}
// A表
inline void AattTablet(int step1, int step2) {
	vector<pDouble> que;

	char filename[120];
	sprintf(filename, "%sData.%dK%dAatt.txt", dirname, step1 / 1000, step2);
	FILE *f = fopen(filename, "w");
	for (int p = 0; p < pSphere; p++) {
		for (int q = 0; q < pSphere; q++) {
			double res;
			double ang = angle(p, q);
			pair<bool, double> status = findx(que, make_pair(ang, 0.));
			if (status.first)
				res = status.second;
			else {
				// 这里应该使用mut作为相互作用能
				res = mpiAvgPropWrtPot(step1, step2, FCHlen, interAng, ang,
				                       Boltzmann, InvBoltz, att, false);
				insertx(que, make_pair(ang, res));
			}
			A[p][q] = res;
			if (rk == 0)
				fprintf(f, "%lf ", res);
		}
		if (rk == 0)
			fprintf(f, "\n");
	}
}

// --------------------具体积分--------------------
// 依赖: f, V表
// 用于: 迭代, 性质
inline void calcIrep() {
	for (int i = 0; i < pSphere; i++) {
		Irep[i] = 0.0;
		for (int j = 0; j < pSphere; j++)
			Irep[i] += f[j] * sin(calPhi(j)) * dAng * dAng * V[i][j];
	}
}
// 依赖: f, A表
// 用于: 迭代, 性质
inline void calcIatt() {
	for (int i = 0; i < pSphere; i++) {
		Iatt[i] = 0.0;
		for (int j = 0; j < pSphere; j++)
			Iatt[i] += f[j] * sin(calPhi(j)) * dAng * dAng * A[i][j];
	}
}
// 依赖: f
// 用于: 迭代, 性质
inline double orderParameter() {
	double res = 0.;
	for (int i = 0; i < pSphere; i++) {
		double phi = calPhi(i);
		res += .5 * (3. * cos(phi) * cos(phi) - 1.) * f[i] * sin(phi) * dAng *
		       dAng;
	}
	return res;
}
// 依赖: f
// 用于: 性质
inline void calcfDer() {
	fDer[0] = fDer[pAng] = 0.;
	for (int i = 1; i < pAng; i++) // 在相邻的phi角处对f数值求导
		fDer[i] = (f[ind[i + 1]] - f[ind[i - 1]]) / (2. * dAng);
}

void allProperties(char name[]) {
	calcIrep();
	calcIatt();
	calcfDer();
	double IdInt = 0;
	for (int i = 0; i < pSphere; i++) {
		if (fabs(f[i]) < EPS)
			continue;
		IdInt += f[i] * log(f[i]) * sin(calPhi(i)) * dAng * dAng;
	}

	double VInt = 0;
	for (int i = 0; i < pSphere; i++) {
		VInt += Irep[i] * f[i] * sin(calPhi(i)) * dAng * dAng;
	}

	double AInt = 0;
	for (int i = 0; i < pSphere; i++) {
		AInt += Iatt[i] * f[i] * sin(calPhi(i)) * dAng * dAng;
	}

	double KtInt = 0;
	for (int u1 = 0; u1 < pSphere; u1++) {
		double res = 0;
		for (int u2 = 0; u2 < pSphere; u2++) {
			double u2y = sin(calPhi(u2)) * sin(calTheta(u2));
			int u2phi = ((u2 == 0) ? 0 : (u2 - 1) / (2 * pAng) + 1);
			double tmp = f[u1] * fDer[u2phi] * M[0][u1][u2] * u2y *
			             sin(calPhi(u1)) * sin(calPhi(u2)) * pow(dAng, 4);
			KtInt += tmp;
			// if (fabs(tmp) > 1e3)
			printf("%6.6lf %d %6.6lf %6.6lf %6.6lf\n", f[u1], u2phi,
			       fDer[u2phi], M[0][u1][u2], tmp);
			// res += tmp;
		}
		// printf("%lf\n", res);
	}
	KtInt *= -.5 * eta;
	printf("KtInt %6.6lf\n", KtInt);

	double K2Int = 0;
	for (int u1 = 0; u1 < pSphere; u1++)
		for (int u2 = 0; u2 < pSphere; u2++) {
			double u1y = sin(calPhi(u1)) * sin(calTheta(u1));
			double u2y = sin(calPhi(u2)) * sin(calTheta(u2));
			int u1phi = ((u1 == 0) ? 0 : (u1 - 1) / (2 * pAng) + 1);
			int u2phi = ((u2 == 0) ? 0 : (u2 - 1) / (2 * pAng) + 1);

			K2Int += fDer[u1phi] * fDer[u2phi] * M[1][u1][u2] * u1y * u2y *
			         sin(calPhi(u1)) * sin(calPhi(u2)) * pow(dAng, 4);
		}
	K2Int *= -.5 * eta;

	// for (int i = 0; i < pAng; i++)
	// 	printf("%lf ", fDer[i]);
	// printf("\n");
	// printf("%lf %lf %lf %lf %lf\n", IdInt, VInt, AInt, KtInt, K2Int);

	double helixPitch = (fabs(KtInt) < 1e-8) ? INF : 2 * M_PI * K2Int / KtInt;
	double q;
	if (fabs(K2Int) < EPS && fabs(KtInt) < EPS) {
		q = 0.;
	} else if (fabs(K2Int) < EPS) {
		q = INF;
	} else {
		q = KtInt / K2Int;
	}

	double freeEneCho =
	    log(eta) - 1 + IdInt + G(eta) * VInt + .5 * eta * AInt - .5 * KtInt * q;
	double freeEneNem = log(eta) - 1 + IdInt + G(eta) * VInt + .5 * eta * AInt;
	double pressure = eta + VInt * 2 * eta * eta * (2 - eta) / pow(1 - eta, 3) +
	                  .5 * eta * eta * AInt - .5 * eta * KtInt * q;

	double chemPotential =
	    log(eta) + IdInt +
	    VInt * eta * (3 * eta * eta - 9 * eta + 8) / pow(1 - eta, 3) +
	    eta * AInt - KtInt * q;

	double orderParam = orderParameter();

	// OUTPUT!

	bool flag = true;
	FILE *file = fopen(name, "r");
	if (file == NULL)
		flag = false;

	file = fopen(name, "a+");
	if (!flag)
		fprintf(file, "%15s %15s %15s %15s %15s %15s %15s %15s %15s\n",
		        var2str(eta), var2str(freeEneCho), var2str(freeEneNem),
		        var2str(pressure), var2str(chemPotential), var2str(orderParam),
		        var2str(K2Int), var2str(KtInt), var2str(helixPitch));
	fprintf(file,
	        "%15.6lf %15.6lf %15.6lf %15.6lf %15.6lf %15.6lf %15.6lf %15.6lf "
	        "%15.6lf\n",
	        eta, freeEneCho, freeEneNem, pressure, chemPotential, orderParam,
	        K2Int, KtInt, helixPitch);
	fclose(file);
}

// --------------------实际操作调用--------------------
// DFT process, type = 0: isotropic, type = 1: nematic
int iteration(int type) {
	char logname[100];
	sprintf(logname, "%siter.log", dirname);

	FILE *file = fopen(logname, "a+");
	fprintf(file, "ITERATION START\n");

	double dev;

	// ODF 初始化
	for (int i = 0; i < pSphere; i++) {
		if (type == 0) {
			// average trial
			f[i] = 1.0 / (4 * M_PI);
		} else if (type == 1) {
			// Onsager trial
			f[i] = onsagerTrial(4.0, calPhi(i));
		}
	}

	calcIrep();
	calcIatt();

	int count = 0; // 迭代次数
	// printf("OP: %lf\pAng", orderParameter());
	while (true) {
		count++;
		fprintf(file, "ETA=%lf, NUM: %d, ", eta, count);

		// 对指数做位移, 防止直接零化, 方便归一化运算
		double mi = 1e6, ma = 0;
		for (int i = 0; i < pSphere; i++) {
			double tmp = 2 * G(eta) * Irep[i] + eta * Iatt[i];
			// double tmp = 2 * G(eta) * Irep[i];
			if (tmp > ma)
				ma = tmp;
			if (tmp < mi)
				mi = tmp;
		}
		double shift = (ma + mi) * .5f;
		// printf("%lf %lf %lf\n", ma, mi, shift);

		// break;

		for (int i = 0; i < pSphere; i++) {
			double tmp = -(2 * G(eta) * Irep[i] + eta * Iatt[i] - shift);
			// double tmp = -(2 * G(eta) * Irep[i] - shift);
			// printf("tmp %lf\n", tmp);
			g[i] = exp(tmp);
		}

		// 计算分母
		double denominator = 0.0;
		for (int i = 0; i < pSphere; i++)
			denominator += g[i] * sin(calPhi(i)) * dAng * dAng;
		// printf("denom: %10lf \pAng", denominator);

		// 归一化, 计算两次迭代误差
		double norm = 0.0;
		for (int i = 0; i < pSphere; i++) {
			g[i] /= denominator;
			norm += (f[i] - g[i]) * (f[i] - g[i]);
			f[i] = g[i];
		}

		calcIrep();
		calcIatt();

		dev = sqrt(norm);
		fprintf(file, "OP: %.5lf, ", orderParameter());
		fprintf(file, "DEV: %.5lf\n", dev);

		if (dev < 0.001)
			break;

		if (count > 250) {
			fprintf(file, "NOT Convergent\n");
			// 迭代不收敛
			break;
		}

		// char fname[200];
		// sprintf(fname, "%sODF/%.3lf_%d.txt", dirname, eta, count);
		// dump1d(fname, f);
	}

	fclose(file);
	return count;
}
// 依赖:
// 用于: M表, A表, V表
void preDataCal(int argc, char *argv[]) {
	// int step1 = char2num(argv[1]);
	// int step2 = char2num(argv[2]);
	// 经过测试之后，就用这个参数
	int step1 = 200;
	int step2 = 100;

	if (rk == 0) {
		t1 = MPI_Wtime();
	}

	// 计算M, A, V表
	// printf("here\n");
	MTablet(step1, step2, 1);
	MTablet(step1, step2, 2);
	AattTablet(step1, step2);
	VexTablet();

	// 在主节点输出数据
	if (rk == 0) {
		t2 = MPI_Wtime();
		printf("clock time = %f\n", t2 - t1);
	}
}
// 依赖: M表, A表, V表, iteration
// 用于: 迭代并计算性质
void mainProcess(int argc, char *argv[]) {
	// double minEta = char2num(argv[1]);
	// double maxEta = char2num(argv[2]);
	// double stepEta = char2num(argv[3]);
	// 测试之后, 用这个参数
	double minEta = 0.2;
	double maxEta = 0.2;
	double stepEta = 0.001;
	// printf("here%lf", minEta);
	// read A,V,M data
	char fname[120];
	{
		sprintf(fname, "%sData.Vex.txt", dirname);
		read2d(fname, V);

		// sprintf(fname, "%sData.500K20Aatt.txt", dirname);
		// sprintf(fname, "%sData.5K100Aatt.txt", dirname);
		sprintf(fname, "%sData.20K100Aatt.txt", dirname);
		read2d(fname, A);

		// sprintf(fname, "%sData.5K1M1.txt", dirname);
		// sprintf(fname, "%sData.5K100M1.txt", dirname);
		sprintf(fname, "%sData.20K100M1.txt", dirname);
		read2d(fname, M[0]);

		// sprintf(fname, "%sData.5K1M2.txt", dirname);
		// sprintf(fname, "%sData.5K100M2.txt", dirname);
		sprintf(fname, "%sData.20K100M2.txt", dirname);
		read2d(fname, M[1]);
	}

	double etalist[1000];
	int l = 0;
	for (double et = minEta; et <= maxEta; et += stepEta)
		etalist[l++] = et;

	char optname[120], isoname[120];
	sprintf(optname, "%sProp.optimal.txt", dirname);
	sprintf(isoname, "%sProp.isotropic.txt", dirname);

	for (int i = 0; i < l; i++) {
		eta = etalist[i];
		sprintf(fname, "%sODF/%.4lf.txt", dirname, eta);

		// Cholesteric相的数据
		bool flag = read1d(fname, f); // 直接读入ODF, 若不存在则重新计算
		// bool flag = false;
		if (flag)
			printf("eta:%.4lf ODF exists\n", eta);
		else {
			// 读入失败, 不存在文件, 则重新迭代计算
			int step = iteration(1); // 初始用Onsager trial, 记录收敛步数
			dump1d(fname, f);        // 暂存ODF到文件中
			printf("eta:%.4lf ODF non, recal step:%d\n", eta, step);
		}
		allProperties(optname);

		// Isotropic相的数据
		for (int j = 0; j < pSphere; j++)
			f[j] = 1.0 / (4 * M_PI);
		allProperties(isoname);
	}
}

// --------------------测试模块--------------------
// 测试f的积分归一性
void testIntegralValidacity() {
	double sum = 0.0;

	for (int i = 0; i < pSphere; i++)
		for (int j = 0; j < pSphere; j++)
			sum += pow(dAng, 4) * sin(calPhi(i)) * sin(calPhi(j)) * f[i] * f[j];

	printf("%.10lf\n", sum);
}
// 测试disSegSeg函数, 用csl.cpp对拍
void checkDisSegSeg() {
	int n = 10;
	double len;
	int k = 10;
	FILE *f = fopen("dis.txt", "w");
	fprintf(f, "%d\n", n);

	for (int i = 0; i < n; i++) {
		// int o1 = 1;
		// int o2 = 1;
		// int o3 = 0;

		// int s11 = 0;
		// int s12 = 0;
		// int s13 = 0;

		// int s21 = 1;
		// int s22 = 0;
		// int s23 = 0;

		int o1 = rand() % k + 1;
		int o2 = rand() % k + 1;
		int o3 = rand() % k + 1;

		int s11 = rand() % k;
		int s12 = rand() % k;
		int s13 = rand() % k;

		int s21 = rand() % k;
		int s22 = rand() % k;
		int s23 = rand() % k;

		int t11 = s11 + o1;
		int t12 = s12 + o2;
		int t13 = s13 + o3;

		int t21 = s21 + o2;
		int t22 = s22 + o3;
		int t23 = s23 + o1;

		Vector3d u1(o1, o2, o3);
		Vector3d u2(o2, o3, o1);
		len = u1.norm();
		u1 /= len;
		u2 /= len;

		fprintf(f, "%d %d %d %d %d %d\n", s11, s12, s13, t11, t12, t13);
		fprintf(f, "%d %d %d %d %d %d\n", s21, s22, s23, t21, t22, t23);

		Vector3d r1(s11, s12, s13);
		r1 += u1 / 2. * len;
		Vector3d r2(s21, s22, s23);
		r2 += u2 / 2. * len;

		printf("%lf\n", pow(disSegSeg(len, r1, u1, len, r2, u2), 2));
	}
}
// 用可视化测试M表, A表, V表
void visualTest(int argc, char *argv[]) {
	/**
	for (int i = 0; i < pSphere; i++) {
	    f[i] = onsagerTrial(4.0, calPhi(i));
	    // f[i] = 1 / 4.0 / M_PI;
	    // f[i] = sin(calPhi(i));
	    // f[i] = angle(0, i);
	}

	calcfDer();

	int ind[pAng + 1];
	ind[0] = 0;
	for (int i = 0; i < pAng; i++)
	    ind[i + 1] = 1 + 2 * pAng * i;

	f[0] = f[pSphere - 1] = fDer[0];

	for (int i = 1; i <= pAng - 1; i++) {
	    for (int j = ind[i]; j < ind[i + 1]; j++)
	        f[j] = fDer[i];
	}
	**/

	char fname[120];
	{
		sprintf(fname, "%sData.Vex.txt", dirname);
		read2d(fname, V);

		sprintf(fname, "%sData.500K20Aatt.txt", dirname);
		read2d(fname, A);

		sprintf(fname, "%sData.100000K1M1.txt", dirname);
		read2d(fname, M[0]);

		sprintf(fname, "%sData.500K20M2.txt", dirname);
		read2d(fname, M[1]);
	}

	double test[pSphere];

	for (int i = 0; i < pSphere; i++) {
		test[i] = 0.;
		for (int j = 0; j < pSphere; j++) {
			test[i] += M[0][i][j];
		}
	}

	sprintf(fname, "%sODF/test_1.txt", dirname);
	dump1d(fname, test);
}
// 测试Mk计算结果的可交换性
void AvgForMkTest(int argc, char *argv[]) {
	// for (int i = 0; i <= pAng; i++) {
	// 	double ang = i * dAng - M_PI / 2.;
	// 	Vector3d u1(0, 0, 0);
	// 	Vector3d u2(sin(ang), 0., cos(ang));
	// 	double res12 =
	// 	    mpiAvgForMk(125000, 1, 1, 10, M_PI / 6., u1, u2, fch, false);
	// 	double res21 =
	// 	    mpiAvgForMk(125000, 1, 1, 10, M_PI / 6., u2, u1, fch, false);
	// 	printf("%lf %lf %lf\n", ang, res12, res21);
	// }

	double phi = 0;
	for (int i = 0; i < 2 * pAng; i++) {
		double theta = i * dAng;
		Vector3d u(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi));
		double res = mpiAvgForMk(5000, 1, 1, 10, M_PI / 6., Vector3d(0, 0, 1.),
		                         u, fch, false);
		printf("%d %lf\n", i, res);
	}
}
// 测试fch势的对称性及和Goossen势的差异
void energyTest(int argc, char *argv[]) {
	FCH chain[2];

	chain[0].setInterAngle(M_PI / 6.);
	chain[1].setInterAngle(M_PI / 6.);

	double z = char2num(argv[1]);

	for (double r = 1.2; r < 1.6; r += 0.1) {
		chain[1].setCentroid(Vector3d(0, r, z));

		double ene[pAng + 1];

		double Goossen[pAng + 1];

		for (int i = 0; i <= pAng; i++) {
			ene[i] = 0.;
			double ang = i * dAng - M_PI / 2.;
			chain[1].setOrientation(Vector3d(sin(ang), 0., cos(ang)));
			ene[i] = nonBondPotential(chain[0], chain[1], fch);

			Goossen[i] = -cos(ang) * sin(ang) * r;
		}

		char name[100];
		sprintf(name, "%.3lf.txt", r);
		FILE *f = fopen(name, "w");
		for (int i = 0; i <= pAng; i++)
			fprintf(f, "%lf\n", ene[i]);
		fclose(f);
	}
}

void printODF(double ETA) {
	char fname[120], oname[120];
	sprintf(fname, "%sODF/%.4lf.txt", dirname, ETA);
	sprintf(oname, "%s/%.4lf.txt", dirname, ETA);
	read1d(fname, f);
	outODF(oname);
}

// 运行方法:
// 1.用preDataCal计算AVM的数据
// 2.读入预数据, 然后用mainProcess计算ODF
// 3.读入ODF, 然后计算各项性质
// **即mainProcess的数据都必须是外部读入的

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nm);
	MPI_Comm_rank(MPI_COMM_WORLD, &rk);

	// 初始化坐标
	initCoord();

	// 读入分子构型参数
	FCHlen = char2num(argv[1]);
	interAng = char2num(argv[2]) / 180. * M_PI;

	// 创建数据存储文件夹
	sprintf(dirname, "DFT/%d%2.1lf/", FCHlen, interAng * 180 / M_PI);
	char cmd[100];
	sprintf(cmd, "mkdir -p %s", dirname);
	system(cmd);
	sprintf(cmd, "mkdir -p %sODF", dirname);
	system(cmd);

	// preDataCal(argc, argv);
	mainProcess(argc, argv);

	MPI_Finalize();
	return 0;
}
