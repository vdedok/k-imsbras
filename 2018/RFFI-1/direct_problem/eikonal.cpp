#include <iostream>
#include <stdio.h>
#include <math.h>

//#define DEBUG
//#define RUN_TESTS

const int SIZE = 100;
const double DELTA = 1.0;
const double ZERO = 0.0000000001;
const double INFTY = 1000000000000;
const double TEST_ZERO = 0.00001;
const double REGION_SIZE = SIZE * DELTA - ZERO;
const int CHECK_STEP = 10;
const double SHOT_EPS = 0.01;
const double SPEED = 5.0;

struct TraceData {
	double x0, y0, z0;
	double x1, y1, z1;
	double nx, ny, nz;
	double length;
	double time;
};

struct DirectProblemData {
	double x0, y0, z0;
	double x1, y1, z1;
	double length;
	double time;
};

struct TestData {
	double x0, y0, z0;
	double vx, vy, vz;
	double x1, y1, z1;
	double nx, ny, nz;
	double length;
};


class DirectEikonal {
public:
	DirectEikonal();
	DirectProblemData tracePath(double x0, double y0, double z0, double x1, double y1, double z1);
	DirectProblemData tracePath2(double x0, double y0, double z0, double x1, double y1, double z1);
	DirectProblemData traceTrajectory(double x0, double y0, double z0, double vx, double vy, double vz);
private:
	double NData[SIZE][SIZE][SIZE];

	void init();
	void setupNData();
	bool runTest();

	TraceData trace(double x0, double y0, double z0, double vx, double vy, double vz);
};

DirectEikonal::DirectEikonal() {
	init();
#ifdef RUN_TESTS
	if (!runTest())
		fprintf(stderr, "Tests failed. Exit...\n");
	else
		fprintf(stderr, "Tests passed\n");
#endif
	setupNData();
}

void DirectEikonal::init() {
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			for (int k = 0; k < SIZE; k++)
				NData[i][j][k] = 1.0;
		}
	}
}

void DirectEikonal::setupNData() {
#ifdef DEBUG
	for (int i = 0; i < SIZE; i++) {
		NData[5][i][0] = 2.0;
		NData[6][i][0] = 2.0;
		NData[7][i][0] = 2.0;
		NData[8][i][0] = 2.0;
		NData[9][i][0] = 2.0;
	}
#endif
	int x0, y0, z0;
	int x1, y1, z1;
	int r0, r1;

	x0 = SIZE / 3;
	y0 = SIZE / 3;
	z0 = SIZE / 2;
	r0 = SIZE / 6;
	r0 = r0 * r0;

	x1 = 2 * SIZE / 3;
	y1 = 2 * SIZE / 3;
	z1 = 2 * SIZE / 3;
	r1 = SIZE / 6;
	r1 = r1 * r1;

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			for (int k = 0; k < SIZE; k++) {
				int r;
				r = (i - x0) * (i - x0) + (j - y0) * (j - y0) + (k - z0) * (k - z0);
				if (r < r0) {
					NData[i][j][k] = 1.2;
					continue;
					//fprintf(stderr, "1.2\n");
				}
				r = (i - x1) * (i - x1) + (j - y1) * (j - y1) + (k - z1) * (k - z1);
				if (r < r1) {
					NData[i][j][k] = 1.3;
					//fprintf(stderr, "1.2: %f, %f, %f\n", (double)i / SIZE, (double)j / SIZE, (double)k / SIZE);
				}
 			}
		}
	}

}

TraceData DirectEikonal::trace(double x0, double y0, double z0, double vx, double vy, double vz) {
	TraceData traceData;

	int pX, pY, pZ;
	pX = x0 / DELTA;
	pY = y0 / DELTA;
	pZ = z0 / DELTA;

#ifdef DEBUG
	fprintf(stderr, "x0: %f, y0: %f, z0: %f\n", x0, y0, z0);
	fprintf(stderr, "INT: x: %d, y: %d, z: %d\n", pX, pY, pZ);
#endif

	//x = x0 + t * vx;
	int dx, dy, dz;
	if (vx > 0.0)
		dx = 1;
	else
		dx = -1;

	if (vy > 0.0)
		dy = 1;
	else
		dy = -1;

	if (vz > 0.0)
		dz = 1;
	else
		dz = -1;

	double planeX0 = DELTA * pX;
	double planeX1 = DELTA * (pX + dx);

	double planeY0 = DELTA * pY;
	double planeY1 = DELTA * (pY + dy);

	double planeZ0 = DELTA * pZ;
	double planeZ1 = DELTA * (pZ + dz);

#ifdef DEBUG
	fprintf(stderr, "planeX0: %f, planeX1: %f\n", planeX0, planeX1);
	fprintf(stderr, "planeY0: %f, planeY1: %f\n", planeY0, planeY1);
	fprintf(stderr, "planeZ0: %f, planeZ1: %f\n", planeZ0, planeZ1);
#endif

	double v = sqrt(vx * vx + vy * vy + vz * vz);
	vx = vx / v;
	vy = vy / v;
	vz = vz / v;

#ifdef DEBUG
	fprintf(stderr, "vx: %f, vy: %f, vz: %f\n", vx, vy, vz);
#endif

	double t = INFTY;
	double nx = 0, ny = 0, nz = 0;
	double px = 0, py = 0, pz = 0;

	if (fabs(vx) > ZERO) {
		double t0 = (planeX0 - x0) / vx;
		double t1 = (planeX1 - x0) / vx;

#ifdef DEBUG
		fprintf(stderr, "X: t0: %f, t1: %f\n", t0, t1);
#endif

		if ((t0 > 0.0) && (t0 < t)) {
			t = t0;
			px = x0 + t * vx;
			py = y0 + t * vy;
			pz = z0 + t * vz;

			if (vx > 0.0)
				nx = -1.0;
			else
				nx = 1.0;
			ny = 0.0;
			nz = 0.0;
		}

		if ((t1 > 0.0) && (t1 < t)) {
			t = t1;
			px = x0 + t * vx;
			py = y0 + t * vy;
			pz = z0 + t * vz;

			if (vx > 0.0)
				nx = -1.0;
			else
				nx = 1.0;
			ny = 0.0;
			nz = 0.0;
		}
	}

	if (fabs(vy) > ZERO) {
		double t0 = (planeY0 - y0) / vy;
		double t1 = (planeY1 - y0) / vy;

#ifdef DEBUG
		fprintf(stderr, "Y: t0: %f, t1: %f\n", t0, t1);
#endif

		if ((t0 > 0.0) && (t0 < t)) {
			t = t0;
			px = x0 + t * vx;
			py = y0 + t * vy;
			pz = z0 + t * vz;

			nx = 0.0;
			if (vy > 0.0)
				ny = -1.0;
			else
				ny = 1.0;
			nz = 0.0;
		}

		if ((t1 > 0.0) && (t1 < t)) {
			t = t1;
			px = x0 + t * vx;
			py = y0 + t * vy;
			pz = z0 + t * vz;

			nx = 0.0;
			if (vy > 0.0)
				ny = -1.0;
			else
				ny = 1.0;
			nz = 0.0;
		}
	}

	if (fabs(vz) > ZERO) {
		double t0 = (planeZ0 - z0) / vz;
		double t1 = (planeZ1 - z0) / vz;

#ifdef DEBUG
		fprintf(stderr, "Z: t0: %f, t1: %f\n", t0, t1);
#endif

		if ((t0 > 0.0) && (t0 < t)) {
			t = t0;
			px = x0 + t * vx;
			py = y0 + t * vy;
			pz = z0 + t * vz;

			nx = 0.0;
			ny = 0.0;
			if (vz > 0.0)
				nz = -1.0;
			else
				nz = 1.0;
		}

		if ((t1 > 0.0) && (t1 < t)) {
			t = t1;
			px = x0 + t * vx;
			py = y0 + t * vy;
			pz = z0 + t * vz;

			nx = 0.0;
			ny = 0.0;
			if (vz > 0.0)
				nz = -1.0;
			else
				nz = 1.0;
		}
	}

#ifdef DEBUG
	fprintf(stderr, "x0: %f, y0: %f, z0: %f\n", x0, y0, z0);
	fprintf(stderr, "px: %f, py: %f, pz: %f\n", px, py, pz);
	fprintf(stderr, "nx: %f, ny: %f, nz: %f\n", nx, ny, nz);
#endif

	traceData.x0 = x0;
	traceData.y0 = y0;
	traceData.z0 = z0;
	traceData.x1 = px;
	traceData.y1 = py;
	traceData.z1 = pz;
	traceData.nx = nx;
	traceData.ny = ny;
	traceData.nz = nz;
	traceData.length = sqrt(
			(px - x0) * (px - x0) + (py - y0) * (py - y0)
					+ (pz - z0) * (pz - z0));

	int pX0, pY0, pZ0;
	pX0 = (x0 + px) / 2.0 / DELTA;
	pY0 = (y0 + py) / 2.0 / DELTA;
	pZ0 = (z0 + pz) / 2.0 / DELTA;

	// v = c / n = 1 / n
	// t = s / v = s * n
	traceData.time = traceData.length * NData[pX0][pY0][pZ0];

#ifdef DEBUG
	fprintf(stderr, "len: %f\n", traceData.length);
#endif

	return traceData;
}

bool DirectEikonal::runTest() {
	TestData tests[]  = {
			{0.0, 0.2, 0.0, 1.0, 0.2, 0.0,    1.0, 0.4, 0.0, -1.0, 0.0, 0.0, 1.019804},
			{1.0, 0.4, 0.0, -1.0, -0.2, 0.0,  0.0, 0.2, 0.0, 1.0, 0.0, 0.0, 1.019804},
			{0.0, 0.2, 0.0, 0.2, 0.8, 0.0,    0.2, 1.0, 0.0, 0.0, -1.0, 0.0, 0.824621},
			{0.2, 1.0, 0.0, -0.2, -0.8, 0.0,  0.0, 0.2, 0.0, 1.0, 0.0, 0.0,  0.824621}
	};

	TraceData result;
	int testSize = sizeof(tests) / sizeof(TestData);

	for (int i = 0; i < testSize; i++) {
		result = trace(tests[i].x0, tests[i].y0, tests[i].z0, tests[i].vx, tests[i].vy, tests[i].vz);
		if ((fabs(result.x1 - tests[i].x1) > TEST_ZERO) ||
			(fabs(result.y1 - tests[i].y1) > TEST_ZERO) ||
			(fabs(result.z1 - tests[i].z1) > TEST_ZERO) ||
			(fabs(result.nx - tests[i].nx) > TEST_ZERO) ||
			(fabs(result.ny - tests[i].ny) > TEST_ZERO) ||
			(fabs(result.nz - tests[i].nz) > TEST_ZERO) ||
			(fabs(result.length - tests[i].length) > TEST_ZERO)) {
			fprintf(stderr, "Test %d failed\n", i);
			return false;
		}
	}

	return true;
}

DirectProblemData DirectEikonal::traceTrajectory(double x0, double y0, double z0, double vx, double vy, double vz) {
	DirectProblemData result;
	double length = 0;
	double time = 0;
	TraceData traceRay;

	result.x0 = x0;
	result.y0 = y0;
	result.z0 = z0;

	while ((x0 < REGION_SIZE) && (y0 < REGION_SIZE) && (z0 < REGION_SIZE)) {
		traceRay = trace(x0, y0, z0, vx, vy, vz);

		int pX0, pY0, pZ0;
		pX0 = (x0 + traceRay.x1) / 2.0 / DELTA;
		pY0 = (y0 + traceRay.y1) / 2.0 / DELTA;
		pZ0 = (z0 + traceRay.z1) / 2.0 / DELTA;

		int pX1, pY1, pZ1;
		pX1 = (traceRay.x1 - traceRay.nx / 2.0) / DELTA;
		pY1 = (traceRay.y1 - traceRay.ny / 2.0) / DELTA;
		pZ1 = (traceRay.z1 - traceRay.nz / 2.0) / DELTA;

#ifdef DEBUG
		fprintf(stderr, "Sector0: %d, %d, %d\n", pX0, pY0, pZ0);
		fprintf(stderr, "Sector1: %d, %d, %d\n", pX1, pY1, pZ1);
#endif

		x0 = traceRay.x1;
		y0 = traceRay.y1;
		z0 = traceRay.z1;
		length += traceRay.length;
		time += traceRay.time;

#ifdef DEBUG
		fprintf(stderr, "Next point: %f, %f, %f\n", x0, y0, z0);
#endif

		double n1, n2;
		n1 = NData[pX0][pY0][pZ0];
/*
		fprintf(stderr, "Sector0: %d, %d, %d: %f\n", pX0, pY0, pZ0, NData[pX0][pY0][pZ0]);
		fprintf(stderr, "Sector1: %d, %d, %d: %f\n", pX1, pY1, pZ1, NData[pX1][pY1][pZ1]);
*/

		if ((pX1 < 0) || (pY1 < 0) || (pZ1 < 0) ||
				(pX1 >= SIZE) || (pY1 >= SIZE) || (pZ1 >= SIZE)) {
#ifdef DEBUG
			fprintf(stderr, "Trace finished\n");
#endif
			break;
		}
		double vLength = sqrt(vx * vx + vy * vy + vz * vz);
		vx = vx / vLength * n1;
		vy = vy / vLength * n1;
		vz = vz / vLength * n1;

		n2 = NData[pX1][pY1][pZ1];

		double scalarProd = vx * traceRay.nx + vy * traceRay.ny + vz * traceRay.nz;
		double coeff = sqrt((n2 * n2 - n1 * n1) / (scalarProd * scalarProd) + 1.0) - 1.0;
		coeff = coeff * scalarProd;

		vx = vx + coeff * traceRay.nx;
		vy = vy + coeff * traceRay.ny;
		vz = vz + coeff * traceRay.nz;
#ifdef DEBUG
		fprintf(stderr, "n1: %f, n2: %f\n", n1, n2);
#endif
	}
#ifdef DEBUG
	fprintf(stderr, "Total length: %f, time: %f\n", length, time);
#endif

	result.x1 = x0;
	result.y1 = y0;
	result.z1 = z0;
	result.length = length;
	result.time = time;

	return result;
}

DirectProblemData DirectEikonal::tracePath2(double x0, double y0, double z0,
		double x1, double y1, double z1) {
	DirectProblemData result;
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

	//
	double vx = vx0;
	double vy = vy0;
	double vz = vz0;

	int steps = 0;

	while ((distance > SHOT_EPS) && (steps++ < CHECK_STEP)) {
		DirectProblemData local = traceTrajectory(x0, y0, z0, vx, vy, vz);
		double localDistance = sqrt(
				  (x1 - local.x1) * (x1 - local.x1)
				+ (y1 - local.y1) * (y1 - local.y1)
				+ (z1 - local.z1) * (z1 - local.z1));
		fprintf(stderr, "local: length: %f, time: %f\n", local.length, local.time);
		if (localDistance < distance) {
			distance = localDistance;
			result = local;
		}

		double diff_x = x1 - local.x1;
		double diff_y = y1 - local.y1;
		double diff_z = z1 - local.z1;

		fprintf(stderr, "dest_x: %f, real_x: %f, diff_x: %f\n", x1, local.x1, diff_x);
		fprintf(stderr, "dest_y: %f, real_y: %f, diff_y: %f\n", y1, local.y1, diff_y);
		fprintf(stderr, "dest_z: %f, real_z: %f, diff_x: %f\n", z1, local.z1, diff_z);

//
		fprintf(stderr, "v_old x: %f, y: %f, z: %f\n", vx, vy, vz);
		vx = (x1 - x0) + SPEED * localDistance * diff_x;
		vy = (y1 - y0) + SPEED * localDistance * diff_y;
		vz = (z1 - z0) + SPEED * localDistance * diff_z;

		length = sqrt(vx * vx + vy * vy + vz * vz);
		vx = vx / length;
		vy = vy / length;
		vz = vz / length;
		fprintf(stderr, "v_new x: %f, y: %f, z: %f\n", vx, vy, vz);

		fprintf(stderr, "steps: %d, localDistance: %f\n", steps, localDistance);
	}

	//

	return result;
}

DirectProblemData DirectEikonal::tracePath(double x0, double y0, double z0,
		double x1, double y1, double z1) {
	DirectProblemData result;
	double vx0, vy0, vz0;
	double vx1, vy1, vz1;
	double vx2, vy2, vz2;
	double distance = INFTY;

	double length = 0;

	vx0 = x1 - x0;
	vy0 = y1 - y0;
	vz0 = z1 - z0;
	length = sqrt(vx0 * vx0 + vy0 * vy0 + vz0 * vz0);
	vx0 = vx0 / length;
	vy0 = vy0 / length;
	vz0 = vz0 / length;

	// TODO: FIX iff vx0 and vy0 == 0 both
	vx1 = vy0;
	vy1 = - vx0;
	vz1 = 0;
	length = sqrt(vx1 * vx1 + vy1 * vy1 + vz1 * vz1);
	vx1 = vx1 / length / 100.0;
	vy1 = vy1 / length / 100.0;
	vz1 = vz1 / length / 100.0;

	vx2 = vy0 * vz1 - vz0 * vy1;
	vy2 = vz0 * vx1 - vx0 * vz1;
	vz2 = vx0 * vy1 - vy0 * vx1;
	length = sqrt(vx2 * vx2 + vy2 * vy2 + vz2 * vz2);
	vx2 = vx2 / length / 100.0;
	vy2 = vy2 / length / 100.0;
	vz2 = vz2 / length / 100.0;

	fprintf(stderr, "x0: %f, y0: %f, z0: %f\n", x0, y0, z0);
	fprintf(stderr, "x1: %f, y1: %f, z1: %f\n", x1, y1, z1);
	fprintf(stderr, "vx0: %f, vy0: %f, vz0: %f\n", vx0, vy0, vz0);
	fprintf(stderr, "vx1: %f, vy1: %f, vz1: %f\n", vx1, vy1, vz1);
	fprintf(stderr, "vx2: %f, vy2: %f, vz2: %f\n", vx2, vy2, vz2);

	for (int i = -CHECK_STEP; i <= CHECK_STEP; i++) {
		for (int j = -CHECK_STEP; j <= CHECK_STEP; j++) {
			double vx = vx0 + i * vx1 + j * vx2;
			double vy = vy0 + i * vy1 + j * vy2;
			double vz = vz0 + i * vz1 + j * vz2;

			DirectProblemData local = traceTrajectory(x0, y0, z0, vx, vy, vz);
			double localDistance = sqrt((x1 - local.x1) * (x1 - local.x1) +
										(y1 - local.y1) * (y1 - local.y1) +
										(z1 - local.z1) * (z1 - local.z1));
			if (localDistance < distance) {
				distance = localDistance;
				result = local;
			}
		}
	}

	return result;
}

int main(int argc, char** argv) {
	fprintf(stderr, "Running eikonal solutioner...\n");
	DirectEikonal directProblem;
	fprintf(stderr, "Running direct problem\n");
//	0.480000, 0.350000
	DirectProblemData d = directProblem.traceTrajectory(0.450000 * SIZE, 0.350000 * SIZE, 0.0, 0.0, 0.0, 1.0);
	fprintf(stderr, "x0: %f, y0: %f, z0: %f\n", d.x0, d.y0, d.z0);
	fprintf(stderr, "x1: %f, y1: %f, z1: %f\n", d.x1, d.y1, d.z1);
	fprintf(stderr, "length: %f, time: %f\n", d.length, d.time);
/*
	DirectProblemData d2 = directProblem.tracePath(0, 0.2, 0, 10, 7, 0);
	fprintf(stderr, "D2:\n");
	fprintf(stderr, "x0: %f, y0: %f, z0: %f\n", d2.x0, d2.y0, d2.z0);
	fprintf(stderr, "x1: %f, y1: %f, z1: %f\n", d2.x1, d2.y1, d2.z1);
	fprintf(stderr, "length: %f, time: %f\n", d2.length, d2.time);
*/
	DirectProblemData d2 = directProblem.tracePath2(0.450000 * SIZE, 0.350000 * SIZE, 0.0, 0.400000 * SIZE, 0.320000 * SIZE, 100.0);
		fprintf(stderr, "D2:\n");
		fprintf(stderr, "x0: %f, y0: %f, z0: %f\n", d2.x0, d2.y0, d2.z0);
		fprintf(stderr, "x1: %f, y1: %f, z1: %f\n", d2.x1, d2.y1, d2.z1);
		fprintf(stderr, "length: %f, time: %f\n", d2.length, d2.time);
//	return 0;
	FILE *f = fopen("data.txt", "w");
	if (f == NULL) {
		fprintf(stderr, "Error creating file\n");
		return 0;
	}
	char delim = ' ';
	for (int i = 0; i < SIZE; i++) {
		delim = ' ';
		for (int j = 0; j < SIZE; j++) {
			double px, py;
			px = (double) i * DELTA;
			py = (double) j * DELTA;

			DirectProblemData d = directProblem.traceTrajectory(px, py, 0.0, 0.0, 0.0, 1.0);
			fprintf(f, "%c %f", delim, (double)d.time - 100.0);
			delim = ',';
		}
		fprintf(f, "\n");
	}
	fclose(f);

	fprintf(stderr, "Finish!\n");
	return 0;
}
