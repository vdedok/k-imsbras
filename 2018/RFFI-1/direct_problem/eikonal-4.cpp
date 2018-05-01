#include <iostream>
#include <stdio.h>
#include <math.h>

#define SIZE 1.0
const double ZERO = 0.0000000001;
const double INFTY = 1000000000000.0;

struct TraceData {
	double x0, y0, z0;
	double vx0, vy0, vz0;
	double x1, y1, z1;
	double vx1, vy1, vz1;
	double tau;
	double length;
};

class DirectEikonal4 {
public:
	DirectEikonal4();
	virtual ~DirectEikonal4();

	TraceData tracePath(double x0, double y0, double z0, double x1, double y1, double z1);

	TraceData tracePathDebug(double x0, double y0, double z0, double vx, double vy, double vz,
				double x1, double y1, double z1, int debugMode = 0);
protected:
	virtual double getN(double x, double y, double z) = 0;
	virtual double getLogNx(double x, double y, double z) = 0;
	virtual double getLogNy(double x, double y, double z) = 0;
	virtual double getLogNz(double x, double y, double z) = 0;

private:
	double f1(double t, double x, double y, double z, double px, double py, double pz, double n2);
	double f2(double t, double x, double y, double z, double px, double py, double pz, double n2);
	double f3(double t, double x, double y, double z, double px, double py, double pz, double n2);
	double f4(double t, double x, double y, double z, double px, double py, double pz, double n);
	double f5(double t, double x, double y, double z, double px, double py, double pz, double n);
	double f6(double t, double x, double y, double z, double px, double py, double pz, double n);

	bool isInCube(double x, double y, double z);

	TraceData traceTraectory(double x, double y, double z, double vx, double vy, double vz);
	TraceData traceTraectoryTest(double x, double y, double z, double vx, double vy, double vz);

	double sX, sY, sZ;
};

DirectEikonal4::DirectEikonal4() {
	fprintf(stderr, "DirectEikonal4::DirectEikonal4()\n");

	sX = SIZE;
	sY = SIZE;
	sZ = SIZE;
}

DirectEikonal4::~DirectEikonal4() {
	fprintf(stderr, "DirectEikonal4::~DirectEikonal4()\n");
}

double DirectEikonal4::f1(double t, double x, double y, double z, double px, double py, double pz, double n2) {
	return px / n2;
}

double DirectEikonal4::f2(double t, double x, double y, double z, double px, double py, double pz, double n2) {
	return py / n2;
}

double DirectEikonal4::f3(double t, double x, double y, double z, double px, double py, double pz, double n2) {
	return pz / n2;
}

double DirectEikonal4::f4(double t, double x, double y, double z, double px, double py, double pz, double n) {
	return getLogNx(x, y, z);
}

double DirectEikonal4::f5(double t, double x, double y, double z, double px, double py, double pz, double n) {
	return getLogNy(x, y, z);
}

double DirectEikonal4::f6(double t, double x, double y, double z, double px, double py, double pz, double n) {
	return getLogNz(x, y, z);
}

bool DirectEikonal4::isInCube(double x, double y, double z) {
	if (x > sX)
		return false;
	if (x < -sX)
		return false;
	if (y > sY)
		return false;
	if (y < -sY)
		return false;
	if (z > sZ)
		return false;
	if (z < -sZ)
		return false;

	return true;
}

TraceData DirectEikonal4::tracePath(double x0, double y0, double z0,
			double x1, double y1, double z1) {
	TraceData result;
	double distance = INFTY;
	double vx0, vy0, vz0;

	double length = 0;

	vx0 = x1 - x0;
	vy0 = y1 - y0;
	vz0 = z1 - z0;
	length = sqrt(vx0 * vx0 + vy0 * vy0 + vz0 * vz0);
	vx0 = vx0 / length;
	vy0 = vy0 / length;
	vz0 = vz0 / length;

	return traceTraectory(x0, y0, z0, vx0, vy0, vz0);
}

TraceData DirectEikonal4::tracePathDebug(double x0, double y0, double z0, double vx, double vy, double vz,
				double x1, double y1, double z1, int debugMode) {
	switch (debugMode) {
	case 1: return traceTraectoryTest(x0, y0, z0, vx, vy, vz);
	}

	return traceTraectory(x0, y0, z0, vx, vy, vz);
}

TraceData DirectEikonal4::traceTraectory(double x, double y, double z, double vx, double vy, double vz) {
	TraceData result;

	double K1x, K1y, K1z, K1px, K1py, K1pz;
	double K2x, K2y, K2z, K2px, K2py, K2pz;
	double K3x, K3y, K3z, K3px, K3py, K3pz;
	double K4x, K4y, K4z, K4px, K4py, K4pz;

	double h = 0.001, h2 = h / 2.0;
	double t = 0;
	double x1 = x, y1 = y, z1 = z;
	double n = getN(x, y, z), n2;
	double px = n * vx, py = n * vy, pz = n * vz;

	fprintf(stderr, "S: (%f, %f, %f), P: (%f, %f, %f)\n", x, y, z, px, py, pz);
	int cnt = 0;

	result.x0 = x;
	result.y0 = y;
	result.z0 = z;
	result.vx0 = vx;
	result.vy0 = vy;
	result.vz0 = vz;

	do {
		cnt++;

		x = x1;
		y = y1;
		z = z1;

		n = getN(x, y, z);
		n2 = n * n;

		K1x = h * f1(t, x, y, z, px, py, pz, n2);
		K1y = h * f2(t, x, y, z, px, py, pz, n2);
		K1z = h * f3(t, x, y, z, px, py, pz, n2);
		K1px = h * f4(t, x, y, z, px, py, pz, n2);
		K1py = h * f5(t, x, y, z, px, py, pz, n2);
		K1pz = h * f6(t, x, y, z, px, py, pz, n2);

		K2x = h * f1(t + h2, x + 0.5 * K1x, y + 0.5 * K1y, z + 0.5 * K1z,
				px  + 0.5 * K1px, py  + 0.5 * K1py, pz +  + 0.5 * K1pz, n2);
		K2y = h * f2(t + h2, x + 0.5 * K1x, y + 0.5 * K1y, z + 0.5 * K1z,
				px  + 0.5 * K1px, py  + 0.5 * K1py, pz +  + 0.5 * K1pz, n2);
		K2z = h * f3(t + h2, x + 0.5 * K1x, y + 0.5 * K1y, z + 0.5 * K1z,
				px  + 0.5 * K1px, py  + 0.5 * K1py, pz +  + 0.5 * K1pz, n2);
		K2px = h * f4(t + h2, x + 0.5 * K1x, y + 0.5 * K1y, z + 0.5 * K1z,
				px  + 0.5 * K1px, py  + 0.5 * K1py, pz +  + 0.5 * K1pz, n2);
		K2py = h * f5(t + h2, x + 0.5 * K1x, y + 0.5 * K1y, z + 0.5 * K1z,
				px  + 0.5 * K1px, py  + 0.5 * K1py, pz +  + 0.5 * K1pz, n2);
		K2pz = h * f6(t + h2, x + 0.5 * K1x, y + 0.5 * K1y, z + 0.5 * K1z,
				px  + 0.5 * K1px, py  + 0.5 * K1py, pz +  + 0.5 * K1pz, n2);

		K3x = h * f1(t + h2, x + 0.5 * K2x, y + 0.5 * K2y, z + 0.5 * K2z,
				px  + 0.5 * K2px, py  + 0.5 * K2py, pz +  + 0.5 * K2pz, n2);
		K3y = h * f2(t + h2, x + 0.5 * K2x, y + 0.5 * K2y, z + 0.5 * K2z,
				px  + 0.5 * K2px, py  + 0.5 * K2py, pz +  + 0.5 * K2pz, n2);
		K3z = h * f3(t + h2, x + 0.5 * K2x, y + 0.5 * K2y, z + 0.5 * K2z,
				px  + 0.5 * K2px, py  + 0.5 * K2py, pz +  + 0.5 * K2pz, n2);
		K3px = h * f4(t + h2, x + 0.5 * K2x, y + 0.5 * K2y, z + 0.5 * K2z,
				px  + 0.5 * K2px, py  + 0.5 * K2py, pz +  + 0.5 * K2pz, n2);
		K3py = h * f5(t + h2, x + 0.5 * K2x, y + 0.5 * K2y, z + 0.5 * K2z,
				px  + 0.5 * K2px, py  + 0.5 * K2py, pz +  + 0.5 * K2pz, n2);
		K3pz = h * f6(t + h2, x + 0.5 * K2x, y + 0.5 * K2y, z + 0.5 * K2z,
				px  + 0.5 * K2px, py  + 0.5 * K2py, pz +  + 0.5 * K2pz, n2);

		K4x = h * f1(t + h, x + K3x, y + K3y, z + K3z, px + K3px, py + K3py, pz + K3pz, n2);
		K4y = h * f2(t + h, x + K3x, y + K3y, z + K3z, px + K3px, py + K3py, pz + K3pz, n2);
		K4z = h * f3(t + h, x + K3x, y + K3y, z + K3z, px + K3px, py + K3py, pz + K3pz, n2);
		K4px = h * f4(t + h, x + K3x, y + K3y, z + K3z, px + K3px, py + K3py, pz + K3pz, n2);
		K4py = h * f5(t + h, x + K3x, y + K3y, z + K3z, px + K3px, py + K3py, pz + K3pz, n2);
		K4pz = h * f6(t + h, x + K3x, y + K3y, z + K3z, px + K3px, py + K3py, pz + K3pz, n2);

		t = t + h;
		x1 = x + 1 / 6.0 * (K1x + 2.0 * (K2x + K3x) + K4x);
		y1 = y + 1 / 6.0 * (K1y + 2.0 * (K2y + K3y) + K4y);
		z1 = z + 1 / 6.0 * (K1z + 2.0 * (K2z + K3z) + K4z);
		px = px + 1 / 6.0 * (K1px + 2.0 * (K2px + K3px) + K4px);
		py = py + 1 / 6.0 * (K1py + 2.0 * (K2py + K3py) + K4py);
		pz = pz + 1 / 6.0 * (K1pz + 2.0 * (K2pz + K3pz) + K4pz);

		//fprintf(stderr, "(%f, %f, %f)\n", x, y, z);
	} while (isInCube(x1, y1, z1));
	//} while (false);

	fprintf(stderr, "D: (%f, %f, %f), T: %f, V: (%f, %f, %f)\n", x1, y1, z1, t, px, py, pz);
	fprintf(stderr, "Cnt: %d. Time: %f\n", cnt, t);

	result.x1 = x1;
	result.y1 = y1;
	result.z1 = z1;
	result.vx1 = px;
	result.vy1 = py;
	result.vz1 = pz;
	result.tau = 0;

	return result;
}

TraceData DirectEikonal4::traceTraectoryTest(double x, double y, double z, double vx, double vy, double vz) {
	TraceData result;

	double h = 0.000001;
	double t = 0;
	double x1 = x, y1 = y, z1 = z;
	double n = getN(x, y, z), n2;
	double px = n * vx, py = n * vy, pz = n * vz;

	fprintf(stderr, "S: (%f, %f, %f), P: (%f, %f, %f)\n", x, y, z, px, py, pz);

	int cnt = 0;

	result.x0 = x;
	result.y0 = y;
	result.z0 = z;
	result.vx0 = vx;
	result.vy0 = vy;
	result.vz0 = vz;

	do {
		cnt++;

		x = x1;
		y = y1;
		z = z1;

		n = getN(x, y, z);
		n2 = n * n;

		t = t + h;
		x1 = x + h * f1(t, x, y, z, px, py, pz, n2);
		y1 = y + h * f2(t, x, y, z, px, py, pz, n2);
		z1 = z + h * f3(t, x, y, z, px, py, pz, n2);
		px = px + h * f4(t, x, y, z, px, py, pz, n2);
		py = py + h * f5(t, x, y, z, px, py, pz, n2);
		pz = pz + h * f6(t, x, y, z, px, py, pz, n2);

		//fprintf(stderr, "(%f, %f, %f)\n", x, y, z);
	} while (isInCube(x1, y1, z1));

	fprintf(stderr, "D: (%f, %f, %f), T: %f, V: (%f, %f, %f)\n", x1, y1, z1, t, px, py, pz);
	fprintf(stderr, "Cnt: %d. Time: %f\n", cnt, t);

	result.x1 = x1;
	result.y1 = y1;
	result.z1 = z1;
	result.vx1 = px;
	result.vy1 = py;
	result.vz1 = pz;

	return result;
}

class DirectEikonal4Exp : public DirectEikonal4 {
public:
	DirectEikonal4Exp(double _a, double _b, double _c, double _d,
			double _x0, double _y0, double _z0);
protected:
	virtual double getN(double x, double y, double z);
	virtual double getLogNx(double x, double y, double z);
	virtual double getLogNy(double x, double y, double z);
	virtual double getLogNz(double x, double y, double z);
private:
	DirectEikonal4Exp();
	double a;
	double b;
	double c;
	double d;
	double x0;
	double y0;
	double z0;
};

DirectEikonal4Exp::DirectEikonal4Exp() {
	a = 0;
	b = 0;
	c = 0;
	d = 0;
	x0 = 0;
	y0 = 0;
	z0 = 0;
}

DirectEikonal4Exp::DirectEikonal4Exp(double _a, double _b, double _c, double _d,
		double _x0, double _y0, double _z0) {
	a = _a;
	b = _b;
	c = _c;
	d = _d;
	x0 = _x0;
	y0 = _y0;
	z0 = _z0;
}

double DirectEikonal4Exp::getN(double x, double y, double z) {
	return 1.0 + a * exp(-b * (x - x0) * (x - x0) - c * (y - y0) * (y - y0) - d * (z - z0) * (z - z0));
}

double DirectEikonal4Exp::getLogNx(double x, double y, double z) {
	double nom, denom;
	double e = a * exp(-b * (x - x0) * (x - x0) - c * (y - y0) * (y - y0) - d * (z - z0) * (z - z0) );
	nom = -2.0 * b * (x - x0) * e;
	denom = 1.0 + e;
	return nom / denom;
}

double DirectEikonal4Exp::getLogNy(double x, double y, double z) {
	double nom, denom;
	double e = a * exp(-b * (x - x0) * (x - x0) - c * (y - y0) * (y - y0) - d * (z - z0) * (z - z0) );
	nom = -2.0 * c * (y - y0) * e;
	denom = 1.0 + e;
	return nom / denom;
}

double DirectEikonal4Exp::getLogNz(double x, double y, double z) {
	double nom, denom;
	double e = a * exp(-b * (x - x0) * (x - x0) - c * (y - y0) * (y - y0) - d * (z - z0) * (z - z0) );
	nom = -2.0 * d * (z - z0) * e;
	denom = 1.0 + e;
	return nom / denom;
}

int main() {
	fprintf(stderr, "Running eikonal-4 solutioner...\n");
	double a = 0.3, b = 5.0, c = 5.0, d = 5.0;
	double x0 = 0, y0 = 0, z0 = 0;
	DirectEikonal4Exp directProblem(a, b, c, d, x0, y0, z0);
/*
	Running eikonal-4 solutioner...
	DirectEikonal4::DirectEikonal4()
	S: (-1.000000, 0.200000, 0.000000), P: (1.001655, 0.000000, 0.000000)
	D: (1.000001, -0.108818, 0.000000), T: 2.253099, V: (0.954605, -0.304209, 0.000000)
	Cnt: 2253099. Time: 2.253099
	S: (-1.000000, 0.200000, 0.000000), P: (1.001655, 0.000000, 0.000000)
	D: (1.000397, -0.109010, 0.000000), T: 2.254000, V: (0.954067, -0.304213, 0.000000)
	Cnt: 2254. Time: 2.254000
	Finish...
	DirectEikonal4::~DirectEikonal4()
*/
	directProblem.tracePathDebug(1.00000, -0.108818, 0.000000, -0.954605, 0.304209, 0.000000, 1.0, 0.1, 0.0, 1);
	directProblem.tracePathDebug(1.00000, -0.108818, 0.000000, -0.954605, 0.304209, 0.000000, 1.0, 0.1, 0.0, 0);

	fprintf(stderr, "Finish...\n");
	return 0;
}
