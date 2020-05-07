#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "Eigen/Eigen"
#include "mpi.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>

#define EPS 1e-5
#define sigma 1.0
#define INF 1e18
#define epsilon 1.0
#define var2str(x) (#x)

using namespace std;
using namespace Eigen;

// --------------------mpi params--------------------
int rk, nm;
double t1, t2;

// --------------------numerical setting--------------------
// 距离上1分为20份, 角度上Pi分为40份
const int pDis = 20;
const int pAng = 20;
const double dDis = 1. / pDis;
const double dAng = M_PI / pAng;

// --------------------random--------------------
random_device rd;
mt19937 generator(rd());
uniform_real_distribution<double> uniDis(0, 1);
normal_distribution<double> normDis(0, 1);
auto Uniform = bind(uniDis, generator);
auto Normal = bind(normDis, generator);

// 随机一个球面上的向量
Vector3d randSphereSurface() {
	Vector3d res;
	for (int i = 0; i < 3; i++)
		res(i) = Normal();
	res /= res.norm();
	return res;
}

// 随机一个盒子里的位置向量
Vector3d randBox(double length) {
	Vector3d res;
	for (int i = 0; i < 3; i++)
		res(i) = (Uniform() - .5) * length;
	return res;
}

// --------------------string--------------------
int sgn(double x) {
	if (fabs(x) < EPS)
		return 0;
	if (x > 0)
		return 1;
	return -1;
}

double char2num(char *s) {
	double res = 0;
	int l = strlen(s);
	int flag = 0;
	int i = 0;
	bool sign = false;

	if (s[0] == '-') {
		sign = true;
		i = 1;
	}

	for (; i < l; i++) {
		if (s[i] == '.') {
			flag = l - i - 1;
			continue;
		}
		res = res * 10 + s[i] - '0';
	}
	return res * pow(10, -flag) * (sign ? -1 : 1);
}

// --------------------Numerical Methods--------------------
// 点到线段最短距离
double disPointSeg(Vector3d p, double len, Vector3d r, Vector3d u) {
	double res = 0.;

	if (fabs(len) < EPS)
		res = (p - r).norm();
	Vector3d v1 = u;
	Vector3d v2 = p - (r - len / 2. * u);
	Vector3d v3 = p - (r + len / 2. * u);

	if (v1.dot(v2) < 0)
		res = v2.norm();
	else if (v1.dot(v3) > 0)
		res = v3.norm();
	else
		res = (v1.cross(v2)).norm();

	// printf("%lf\n", res);
	return res;
}

// 线段到线段最短距离
double disSegSeg(double len1, Vector3d r1, Vector3d u1, double len2,
                 Vector3d r2, Vector3d u2) {
	// length of u1 and u2 are already assumed to be unimodular
	double res = 0.;
	Vector3d s = r1 - r2;
	double d = u1.dot(u2);

	double s1 = s.dot(u1);
	double s2 = s.dot(u2);

	double denom = 1 - pow(d, 2);

	double lambda1 = 0.;
	double lambda2 = 0.;

	bool flag = true;

	if (fabs(denom) < EPS) {
		// 平行
		flag = false;
	} else {
		// 非平行
		lambda1 = -s1 + s2 * d;
		lambda2 = s2 - s1 * d;
		lambda1 /= denom;
		lambda2 /= denom;

		// 异面直线不相交
		if ((fabs(lambda1) > len1 / 2.) || (fabs(lambda2) > len2 / 2.))
			flag = false;
		else
			flag = true;
	}
	// printf("%lf %lf\n", lambda1, lambda2);
	// printf("%lf %lf %lf\n", p1[0], p1[1], p1[2]);
	// printf("%lf %lf %lf\n", p2[0], p2[1], p2[2]);

	if (flag) {
		// 异面相交, 直接计算点距离
		res = (s + lambda1 * u1 - lambda2 * u2).norm();
	} else {
		// 异面不相交, 或平行, 点到直线距离
		res = min(min(disPointSeg(r1 - len1 / 2. * u1, len2, r2, u2),
		              disPointSeg(r1 + len1 / 2. * u1, len2, r2, u2)),
		          min(disPointSeg(r2 - len2 / 2. * u2, len1, r1, u1),
		              disPointSeg(r2 + len2 / 2. * u2, len1, r1, u1)));
	}

	return res;
}

// --------------------Molecule and potentials--------------------
class FCH {
  private:
	int n;
	Vector3d orientation; // already normalized
	Vector3d centroid;
	double phi0;       //内旋角
	double interAngle; //质心内旋角

	double ibRadius; // radius of interaction bead
	double bbRadius; // radius of backbone bead
  public:
	FCH() {
		n = 10;
		orientation << 0, 0, 1.0;
		centroid << 0, 0, 0;

		phi0 = 0.0;
		interAngle = M_PI / 6;

		bbRadius = pow(2, 1 / 6.0) * sigma / 2;
		ibRadius = 0.1 * bbRadius;
	}
	void setRandom() {
		orientation = randSphereSurface();
		// phi0 = Uniform() * 2 * M_PI;
	}
	void translate(Vector3d r) { centroid += r; }
	void rotate(Vector3d axis, double angle) {
		axis /= axis.norm();
		AngleAxisd r(angle, axis);

		orientation = r.matrix() * orientation;
	}
	void print() {
		cout << "centroid:" << endl;
		cout << centroid.transpose() << endl << endl;

		cout << "orientation:" << endl;
		cout << orientation.transpose() << endl << endl;
	}

	void setN(int n_) { n = n_; }
	void setOrientation(Vector3d orientation_) { orientation = orientation_; }
	void setCentroid(Vector3d centroid_) { centroid = centroid_; }
	void setPhi0(double phi0_) { phi0 = phi0_; }
	void setInterAngle(double inter) { interAngle = inter; }

	// type = 0 for ib, type = 1 for bb
	Vector3d getCentroid(int type, int index) {
		// get the index bead from down to top
		// backbone bead
		double m = index - (n + 1.0) / 2;
		double rLength = m * bbRadius * 2;
		double phi = m * interAngle + phi0;

		if (type == 1)
			return rLength * orientation + centroid;

		// interaction bead
		Vector3d oriIb(cos(phi), sin(phi), 0);

		double x = orientation[0];
		double y = orientation[1];
		double l = sqrt(x * x + y * y);

		if (sgn(l) != 0) {
			double theta = acos(x / sqrt(x * x + y * y));
			if (y < 0)
				theta = -theta;
			double varphi = acos(orientation[2]);

			oriIb = AngleAxisd(theta, Vector3d(0, 0, 1.0)).matrix() *
			        AngleAxisd(varphi, Vector3d(0, 1.0, 0)).matrix() * oriIb;
		}
		return oriIb * bbRadius + rLength * orientation + centroid;
	}
	double getRadius(int type) {
		if (type)
			return bbRadius;
		return ibRadius;
	}
	double getLength() { return n * bbRadius * 2; }
	Vector3d getCen() { return centroid; }
	Vector3d getOri() { return orientation; }
	int getN() { return n; }
	double getVolume() {
		double h1 = ibRadius * ibRadius / (2 * bbRadius);
		double h2 = ibRadius - h1;
		double Vm = M_PI * (4.0 / 3 * (pow(bbRadius, 3) + pow(ibRadius, 3)) -
		                    h1 * h1 * (bbRadius - h1 / 3) -
		                    h2 * h2 * (ibRadius - h2 / 3));
		return Vm * n;
	}
};

double HSPotential(double sig, double r) {
	if (r >= sig)
		return 0;
	return INF;
}
double LJPotential(double sig, double r) {
	if (fabs(r) < EPS)
		return INF;
	return 4 * (pow(sig / r, 12) - pow(sig / r, 6)) * epsilon;
}
double sLJPotential(double sig, double r) {
	double co = 2.5 * sigma;
	if (r >= co)
		return 0;
	return LJPotential(sig, r) - LJPotential(sig, co);
}
double WCAPotential(double r) {
	if (r >= pow(2, 1 / 6.0) * sigma)
		return 0;
	return LJPotential(sigma, r) + epsilon;
}
#define HS 0
#define LJ 1
#define sLJ 2
#define WCA 3
double Potential(double sig, double r, int type) {
	if (type == HS)
		return HSPotential(sig, r);
	if (type == LJ)
		return LJPotential(sig, r);
	if (type == sLJ)
		return sLJPotential(sig, r);
	if (type == WCA)
		return WCAPotential(r);
	return 0;
}

#define fch 1
#define att 2
#define mut 3 // 大和小的相互作用

double nonBondPotential(FCH &A, FCH &B, int type) {
	int nA = A.getN();
	int nB = B.getN();
	double radA[2] = {A.getRadius(0), A.getRadius(1)};
	double radB[2] = {B.getRadius(0), B.getRadius(1)};

	double res = 0.0;
	double sig;
	double dis;
	// intermolecular
	// conclu: mut是主要的相互作用来源, 相比起来att和HS都是渺小的
	for (int i = 1; i <= nA; i++) {
		for (int j = 1; j <= nB; j++) {
			if (type == fch) {
				for (int s = 0; s < 2; s++)
					for (int t = 0; t < 2; t++) {
						sig = radA[s] + radB[t];
						dis =
						    (A.getCentroid(s, i) - B.getCentroid(t, j)).norm();
						if (s && t) // 都是大球用WCA势
							res += Potential(sig, dis, WCA);
						else // 否则一律用sLJ势
							res += Potential(sig, dis, sLJ);
					}
			} else if (type == HS) {
				// 大小都硬
				for (int s = 0; s < 2; s++)
					for (int t = 0; t < 2; t++) {
						sig = radA[s] + radB[t];
						dis =
						    (A.getCentroid(s, i) - B.getCentroid(t, j)).norm();
						res += Potential(sig, dis, HS);
					}
			} else if (type == att) {
				sig = radA[0] + radB[0];
				dis = (A.getCentroid(0, i) - B.getCentroid(0, j)).norm();
				res += Potential(sig, dis, sLJ);
			} else if (type == mut) {
				// 只算小大之间
				for (int s = 0; s < 2; s++) {
					int t = 1 - s;
					sig = radA[s] + radB[t];
					dis = (A.getCentroid(s, i) - B.getCentroid(t, j)).norm();
					res += Potential(sig, dis, sLJ);
				}
			}
		}
	}

	return res;
}

bool collide(FCH &A, FCH &B) {
	int nA = A.getN();
	int nB = B.getN();
	double bbA = A.getRadius(1);
	double bbB = B.getRadius(1);

	// intermolecular
	for (int i = 1; i <= nA; i++)
		for (int j = 1; j <= nB; j++) {
			double dis = (A.getCentroid(1, i) - B.getCentroid(1, j)).norm();
			if (dis < (bbA + bbB))
				return true;
		}
	return false;
}
double GoossenPotential(FCH &A, FCH &B) {
	Vector3d u1 = A.getOri();
	Vector3d u2 = B.getOri();

	Vector3d dis = B.getCen() - A.getCen();
	Vector3d crs = u1.cross(u2);

	double eps = 10.;
	double u = 1 / pow(dis.norm(), 6);
	double T = u1.dot(u2) * crs.dot(dis);

	return -eps * u * T;
}

function<double(double)> Mayer = [](double x) { return 1. - exp(-x); };
function<double(double)> Boltzmann = [](double x) { return exp(-x); };
function<double(double)> InvBoltz = [](double x) { return -log(x); };
function<double(double)> Identity = [](double x) { return x; };
function<double(double)> Nill = [](double x) { return 0.; };

// 基本废了, 把取向平均了就没用了
// --------------------PotRZFuncs--------------------
// 计算给定r和z情形下, 将取向平均掉的能量
#define BY_AVG 0
#define BY_MIN 1
#define BY_MC 2
pair<double, double> potRZbyAvg(double angle, double r, double z) {
	FCH ref, mov;
	ref.setInterAngle(angle);
	mov.setInterAngle(angle);
	mov.setCentroid(Vector3d(0., r, z));

	double potential[pAng / 2 * (pAng + 1)];
	double wPot[pAng / 2 * (pAng + 1)];
	double wCoef[pAng / 2 * (pAng + 1)];

	for (int k = 0; k < pAng / 2 * (pAng + 1); k++) {
		// parallel partition
		int i = k / (pAng / 2); // varphi = 0-pi
		int j = k % (pAng / 2); // theta = 0-pi/2

		double varphi = i * dAng;
		double theta = j * dAng;
		mov.setOrientation(Vector3d(sin(varphi) * cos(theta),
		                            sin(varphi) * sin(theta), cos(varphi)));
		//  change of coordinate angle
		// mov.setOrientation(Vector3d(cos(theta), sin(theta) * sin(varphi),
		// sin(theta) * cos(varphi)));
		double resPot = 0.;
		double resCoef = 0.;
		for (int r = 0; r < 4 * pAng * pAng; r++) {
			int k1 = r / (2 * pAng); // inter = 0-2pi
			int k2 = r % (2 * pAng); // inter = 0-2pi

			double inter1 = k1 * dAng;
			double inter2 = k2 * dAng;

			ref.setPhi0(inter1);
			mov.setPhi0(inter2);
			double tmp = nonBondPotential(ref, mov, fch);
			// printf("%lf\n", tmp);
			resPot += tmp * exp(-tmp);
			resCoef += exp(-tmp);
		}
		wPot[k] = resPot / (4 * pAng * pAng);
		wCoef[k] = resCoef / (4 * pAng * pAng);
	}

	// integral over c
	double resPot = 0.;
	double resCoef = 0.;
	for (int k = 0; k < pAng / 2 * (pAng + 1); k++) {
		int i = k / (pAng / 2); // varphi = 0-pi
		int j = k % (pAng / 2); // theta = 0-pi/2
		double varphi = i * dAng;
		double theta = j * dAng;

		resPot += wPot[k] * sin(varphi) / 2. / (pAng / 2);
		resCoef += wCoef[k] * sin(varphi) / 2. / (pAng / 2);
		//  change of coordinate angle
		// resPot += wPot[k] * sin(theta) / 2. / (pAng / 2);
		// resCoef += wCoef[k] * sin(theta) / 2. / (pAng / 2);
	}

	return make_pair(resPot, resCoef);
}
pair<double, double> potRZbyMin(double angle, double r, double z) {
	// double r = ir * dDis;
	// double z = iz * dDis;
	FCH ref, mov;
	ref.setInterAngle(angle);
	mov.setInterAngle(angle);
	mov.setCentroid(Vector3d(0., r, z));

	double potential[pAng / 2 * (pAng + 1)];
	double wPot[pAng / 2 * (pAng + 1)];
	double wCoef[pAng / 2 * (pAng + 1)];

	for (int k = 0; k < pAng / 2 * (pAng + 1); k++) {
		// parallel partition
		int i = k / (pAng / 2); // varphi = 0-pi
		int j = k % (pAng / 2); // theta = 0-pi/2

		double varphi = i * dAng;
		double theta = j * dAng;
		mov.setOrientation(Vector3d(sin(varphi) * cos(theta),
		                            sin(varphi) * sin(theta), cos(varphi)));

		double minPot = INF;
		double minCoef = 0.;
		for (int r = 0; r < 4 * pAng * pAng; r++) {

			int k1 = r / (2 * pAng); // inter = 0-2pi
			int k2 = r % (2 * pAng); // inter = 0-2pi

			double inter1 = k1 * dAng;
			double inter2 = k2 * dAng;

			ref.setPhi0(inter1);
			mov.setPhi0(inter2);
			double tmp = nonBondPotential(ref, mov, fch);

			if (tmp < minPot) {
				minPot = tmp;
			}
		}
		// printf("%d, %lf\n", k, minPot);

		wPot[k] = minPot * exp(-minPot);
		wCoef[k] = exp(-minPot);
	}

	// integral over c
	double resPot = 0.;
	double resCoef = 0.;
	for (int k = 0; k < pAng / 2 * (pAng + 1); k++) {
		int i = k / (pAng / 2); // varphi = 0-pi
		int j = k % (pAng / 2); // theta = 0-pi/2
		double varphi = i * dAng;
		double theta = j * dAng;

		resPot += wPot[k] * sin(varphi) / 2. / (pAng / 2);
		resCoef += wCoef[k] * sin(varphi) / 2. / (pAng / 2);
	}

	// if (resPot < 0)
	// printf("%lf %lf %lf\n", r, z, resPot);
	return make_pair(resPot, resCoef);
}
pair<double, double> potRZbyMc(double angle, double r, double z) {
	// double r = ir * dDis;
	// double z = iz * dDis;
	FCH chain[2];
	chain[0].setInterAngle(angle);
	chain[1].setInterAngle(angle);
	chain[1].setCentroid(Vector3d(0., r, z));

	double wPot[pAng / 2 * (pAng + 1)];
	double wCoef[pAng / 2 * (pAng + 1)];

	for (int k = 0; k < pAng / 2 * (pAng + 1); k++) {
		// parallel partition
		int i = k / (pAng / 2); // varphi = 0-pi
		int j = k % (pAng / 2); // theta = 0-pi/2

		double varphi = i * dAng;
		double theta = j * dAng;
		chain[1].setOrientation(Vector3d(
		    sin(varphi) * cos(theta), sin(varphi) * sin(theta), cos(varphi)));

		double minPot = INF;
		double phi[2];
		for (int r = 0; r < 4 * pAng * pAng; r++) {
			int k1 = r / (2 * pAng); // inter = 0-2pi
			int k2 = r % (2 * pAng); // inter = 0-2pi

			double inter1 = k1 * dAng;
			double inter2 = k2 * dAng;

			chain[0].setPhi0(inter1);
			chain[1].setPhi0(inter2);
			double tmp = nonBondPotential(chain[0], chain[1], fch);

			if (tmp < minPot) {
				minPot = tmp;
				phi[0] = inter1;
				phi[1] = inter2;
			}
		}

		int step = 100;
		double cur[2] = {phi[0], phi[1]};
		double res = 0.;
		int cutoff = step / 2;

		for (int t = 0; t < 2; t++)
			chain[t].setPhi0(cur[t]);
		double curE = nonBondPotential(chain[0], chain[1], fch);

		for (int k = 0; k < step; k++) {
			double rdAng[2];
			for (int t = 0; t < 2; t++) {
				rdAng[t] = Uniform() * dAng * 2;
				cur[t] += rdAng[t];
				if (cur[t] > 2 * M_PI)
					cur[t] -= 2 * M_PI;
				chain[t].setPhi0(cur[t]);
			}
			double rdE = nonBondPotential(chain[0], chain[1], fch);
			// printf("%lf\n", exp(curE - rdE));
			double trial = Uniform();
			if (trial <= exp(curE - rdE)) {
				curE = rdE;
			} else {
				for (int t = 0; t < 2; t++) {
					cur[t] -= rdAng[t];
					if (cur[t] < 0)
						cur[t] += 2 * M_PI;
					chain[t].setPhi0(cur[t]);
				}
			}
			if (k >= cutoff)
				res += curE;
		}

		// printf("%d, %lf\n", k, minPot);
		double ene = res / (step - cutoff + 1);
		// printf("%lf\n", ene);

		wPot[k] = ene * exp(-ene);
		wCoef[k] = exp(-ene);
	}

	// integral over c
	double resPot = 0.;
	double resCoef = 0.;
	for (int k = 0; k < pAng / 2 * (pAng + 1); k++) {
		int i = k / (pAng / 2); // varphi = 0-pi
		int j = k % (pAng / 2); // theta = 0-pi/2
		double varphi = i * dAng;
		double theta = j * dAng;

		resPot += wPot[k] * sin(varphi) / 2. / (pAng / 2);
		resCoef += wCoef[k] * sin(varphi) / 2. / (pAng / 2);
	}

	// if (resPot < 0)
	// printf("%lf %lf %lf\n", r, z, resPot);
	return make_pair(resPot, resCoef);
}

// --------------------PotRZFuncs--------------------
// 计算给定FCh螺旋角, 给定夹角, 然后平均均掉位置和内转角的性质
// 把内转角平均掉了, 所以只和夹角有关
// 性质由dealE函数给出, dealPhi函数再处理一次, 两者均是能量的函数
// 注意, 必须在主线程rk = 0上拿结果
double mpiAvgPropWrtPot(int step1, int step2, int FCHlen, double interAng,
                        double gamma, function<double(double)> dealE,
                        function<double(double)> dealPhi, int potType,
                        bool useMC) {
	// chain[0] 放在原点不动
	// chain[1] 进行随机采样
	FCH chain[2];
	chain[0].setInterAngle(interAng);
	chain[1].setInterAngle(interAng);
	chain[0].setN(FCHlen);
	chain[1].setN(FCHlen);

	Vector3d u(0, cos(gamma), sin(gamma));
	chain[1].setOrientation(u);

	// 进行MC的盒子宽度
	double len = chain[0].getLength();
	double sideLength = 2 * len;
	double Vm = chain[0].getVolume();

	double prop = 0.;
	double myProp = 0.;

	if (useMC) {
		// mpi并行循环
		for (int i = rk; i < step1; i += nm) {
			// 随机选取chain[1]的位置
			Vector3d cur;

			// 柱坐标随机
			cur(0) = 0;
			for (int j = 1; j < 3; j++)
				cur(j) = (Uniform() - .5) * sideLength;
			chain[1].setCentroid(cur);

			if (disSegSeg(len, Vector3d(0, 0, 0), Vector3d(0, 0, 1), len, cur,
			              u) > 3.)
				continue;

			// 距离优化

			double avgOverPhi = 0.;
			for (int j = 0; j < step2; j++) {
				// 随机选取两个内转角

				double phi[2];
				for (int k = 0; k < 2; k++) {
					phi[k] = Uniform() * 2 * M_PI;
					chain[k].setPhi0(phi[k]);
				}

				// 计入需要的性质函数
				double ene = nonBondPotential(chain[0], chain[1], potType);
				avgOverPhi += dealE(ene);
			}

			avgOverPhi /= step2;
			myProp += dealPhi(avgOverPhi);
		}
	} else {
		// 这里也可以只在半平面内采样?
		int par = sqrt(step1 / 2);
		double ds = sideLength / (2 * par);

		for (int ir = rk; ir < par; ir += nm) {
			Vector3d cur;
			double ene;
			for (int iz = -par; iz < par; iz++) {
				cur(0) = 0.;
				cur(1) = ir * ds;
				cur(2) = iz * ds;

				// 距离优化
				if (disSegSeg(len, Vector3d(0, 0, 0), Vector3d(0, 0, 1), len,
				              cur, u) > 3.)
					continue;

				chain[1].setCentroid(cur);

				double avgOverPhi = 0.;
				int par2 = sqrt(step2);
				double da = 2 * M_PI / par2;
				for (int i = 0; i < par2; i++)
					for (int j = 0; j < par2; j++) {
						double ang[2] = {i * da, j * da};
						for (int k = 0; k < 2; k++)
							chain[k].setPhi0(ang[k]);
						double ene =
						    nonBondPotential(chain[0], chain[1], potType);
						avgOverPhi += dealE(ene);
					}
				avgOverPhi /= step2;
				double tmp = dealPhi(avgOverPhi);
				myProp += tmp;
				// myProp /= 10.;
			}
		}
	}

	MPI_Reduce(&myProp, &prop, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	// printf("%lf %lf\n", prop,
	//        2 * M_PI * (prop * pow(sideLength, 3) / step1) / Vm);
	// printf("%lf\n", prop);

	// 对dr的积分归一化(积分外面只差一个eta)
	return 2 * M_PI * (prop * pow(sideLength, 3) / step1) / Vm;
}

// t = 1或2
double mpiAvgForMk(int step1, int step2, int t, int FCHlen, double interAng,
                   Vector3d u1, Vector3d u2, int potType, bool useMC) {
	FCH chain[2];
	chain[0].setInterAngle(interAng);
	chain[1].setInterAngle(interAng);
	chain[0].setN(FCHlen);
	chain[1].setN(FCHlen);

	chain[0].setOrientation(u1);
	chain[1].setOrientation(u2);

	// chain[1].rotate(Vector3d(0, 1.0, 0), gamma);

	// 进行MC的盒子宽度
	double len = chain[0].getLength();
	double sideLength = 2 * len;
	double Vm = chain[0].getVolume();

	double prop = 0.;
	double myProp = 0.;

	int cnt = 0;
	int myCnt = 0;

	if (useMC) {
		// mpi并行循环
		for (int i = rk; i < step1; i += nm) {
			// 随机选取chain[1]的位置
			Vector3d cur;

			// 柱坐标随机, cur(1) = r, cur(2) = z
			cur(0) = 0.;
			for (int j = 1; j < 3; j++)
				cur(j) = (Uniform() - .5) * sideLength;
			chain[1].setCentroid(cur);

			// 以step2的数量对内转角取平均
			/*
			double avgOverPhi = 0.;
			for (int j = 0; j < step2; j++) {
			    // 随机选取两个内转角
			    double phi[2];
			    for (int k = 0; k < 2; k++) {
			        phi[k] = Uniform() * 2 * M_PI;
			        chain[k].setPhi0(phi[k]);
			    }

			    // 计入需要的性质函数
			    double ene = nonBondPotential(chain[0], chain[1], potType);
			    avgOverPhi += (exp(-ene) - 1) * pow(cur(2), t);
			}
			avgOverPhi /= step2;
			myProp += cur(1) * avgOverPhi;
			*/

			// 考虑直接不去平均内转角, 只算一个构型
			double ene = nonBondPotential(chain[0], chain[1], potType);
			myProp += cur(1) * (exp(-ene) - 1) * pow(cur(2), t);
		}
	} else {
		// 在半平面内平均采样
		int par = sqrt(step1 / 2);
		double ds = sideLength / (2 * par);
		// printf("%d\n", par);
		// int cnt = 0;
		for (int ir = rk; ir < par; ir += nm) {
			Vector3d cur;
			double ene;
			for (int iz = -par; iz < par; iz++) {
				cur(0) = 0.;
				cur(1) = ir * ds;
				cur(2) = iz * ds;

				// 距离优化
				if (disSegSeg(len, Vector3d(0, 0, 0), u1, len, cur, u2) > 3.)
					continue;

				chain[1].setCentroid(cur);

				double ene = nonBondPotential(chain[0], chain[1], potType);
				double tmp = cur(1) * (exp(-ene / 10) - 1) * pow(cur(2), t);
				myProp += tmp;
				myCnt++;

				// if (fabs(ene) > EPS)
				// 	cnt++;
				// printf("%lf %lf\n", ene, tmp);
			}
		}
		// printf("%lf\n", double(cnt / (par + 1.) / (2 * par + 1.)));
	}
	// printf("%lf\n", prop);

	MPI_Reduce(&myProp, &prop, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&myCnt, &cnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	// 实际计入的构型数量
	// printf("%d ", cnt);
	// 对dr的积分归一化(积分外面只差一个eta)
	return 2 * M_PI * (prop * pow(sideLength, 2) / step1) / Vm;
}

#endif