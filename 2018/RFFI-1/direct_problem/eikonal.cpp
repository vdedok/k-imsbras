#include <iostream>
#include <stdio.h>
#include <math.h>

//#define DEBUG
//#define RUN_TESTS

const int SIZE = 2;
const double DELTA = 1.0;
const double ZERO = 0.0000000001;
const double INFTY = 1000000000000.0;
const double TEST_ZERO = 0.00001;
const double REGION_SIZE = SIZE * DELTA + ZERO;
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

struct Point {
       double x, y, z;
       int p;
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
#ifndef DEBUG
				NData[i][j][k] = 1.0;
#endif

#ifdef DEBUG
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
#endif
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
	fprintf(stderr, "DirectEikonal::trace: x0: %f, y0: %f, z0: %f\n", x0, y0, z0);
	fprintf(stderr, "DirectEikonal::trace: vx: %f, vy: %f, vz: %f\n", vx, vy, vz);
	fprintf(stderr, "DirectEikonal::trace: INT: x: %d, y: %d, z: %d\n", pX, pY, pZ);
#endif

	//x = x0 + t * vx;
	int dx = 0, dy = 0, dz = 0;
	if (vx > ZERO)
		dx = 1;
	else if (vx < -ZERO)
		dx = -1;

	if (vy > ZERO)
		dy = 1;
	else if (vy < -ZERO)
		dy = -1;

	if (vz > ZERO)
		dz = 1;
	else if (vz < -ZERO)
		dz = -1;

#ifdef DEBUG
	fprintf(stderr, "DirectEikonal::trace: dx: %i, dy: %i, dz: %i\n", dx, dy, dz);
#endif

	/*
	double planeX0 = DELTA * pX;
	double planeX1 = DELTA * (pX + dx);

	double planeY0 = DELTA * pY;
	double planeY1 = DELTA * (pY + dy);

	double planeZ0 = DELTA * pZ;
	double planeZ1 = DELTA * (pZ + dz);
	*/

	double planeX0 = floor(x0);
	double planeX1 = ceil(x0);
	if (fabs(round(x0) - x0) < ZERO) {
		planeX0 = round(x0);
		planeX1 = planeX0 + DELTA * dx;
	}
/*
	if (fabs(planeX0 - planeX1) < ZERO)
		planeX1 = DELTA * (pX + dx);
		*/

	double planeY0 = floor(y0);
	double planeY1 = ceil(y0);
	if (fabs(round(y0) - y0) < ZERO) {
		planeY0 = round(y0);
		planeY1 = planeY0 + DELTA * dy;
	}
	/*
	if (fabs(planeY0 - planeY1) < ZERO)
		planeY1 = DELTA * (pY + dy);
*/
	double planeZ0 = floor(z0);
	double planeZ1 = ceil(z0);
	if (fabs(round(z0) - z0) < ZERO) {
		planeZ0 = round(z0);
		planeZ1 = planeZ0 + DELTA * dz;
	}
	/*
	if (fabs(planeZ0 - planeZ1) < ZERO)
		planeZ1 = DELTA * (pZ + dz);
		*/

/*
	if (planeX1 > SIZE)
		planeX1 = SIZE - 1;
	if (planeY1 > SIZE)
		planeY1 = SIZE - 1;
	if (planeZ1 > SIZE)
		planeZ1 = SIZE - 1;
*/

#ifdef DEBUG
	fprintf(stderr, "DirectEikonal::trace: planeX0: %f, planeX1: %f\n", planeX0, planeX1);
	fprintf(stderr, "DirectEikonal::trace: planeY0: %f, planeY1: %f\n", planeY0, planeY1);
	fprintf(stderr, "DirectEikonal::trace: planeZ0: %f, planeZ1: %f\n", planeZ0, planeZ1);
#endif

	double v = sqrt(vx * vx + vy * vy + vz * vz);
	vx = vx / v;
	vy = vy / v;
	vz = vz / v;

#ifdef DEBUG
	fprintf(stderr, "DirectEikonal::trace: vx: %f, vy: %f, vz: %f\n", vx, vy, vz);
#endif

	double t = INFTY;
	double nx = 0, ny = 0, nz = 0;
	double px = 0, py = 0, pz = 0;

	if (fabs(vx) > ZERO) {
		double t0 = (planeX0 - x0) / vx;
		double t1 = (planeX1 - x0) / vx;

#ifdef DEBUG
		fprintf(stderr, "DirectEikonal::trace: X: t0: %f, t1: %f, t: %f\n", t0, t1, t);
#endif

		if ((t0 > ZERO) && (t0 < (t - ZERO))) {
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

		if ((t1 > ZERO) && (t1 < (t - ZERO))) {
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
		fprintf(stderr, "DirectEikonal::trace: Y: t0: %f, t1: %f, t: %f\n", t0, t1, t);
#endif

		if ((t0 > ZERO) && (t0 < (t - ZERO))) {
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

		if ((t1 > ZERO) && (t1 < (t - ZERO))) {
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
		fprintf(stderr, "DirectEikonal::trace: Z: t0: %f, t1: %f, t: %f\n", t0, t1, t);
#endif

		if ((t0 > ZERO) && (t0 < (t - ZERO))) {
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

		if ((t1 > ZERO) && (t1 < (t - ZERO))) {
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
	fprintf(stderr, "DirectEikonal::trace: x0: %f, y0: %f, z0: %f\n", x0, y0, z0);
	fprintf(stderr, "DirectEikonal::trace: px: %f, py: %f, pz: %f\n", px, py, pz);
	fprintf(stderr, "DirectEikonal::trace: nx: %f, ny: %f, nz: %f\n", nx, ny, nz);
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
	fprintf(stderr, "DirectEikonal::trace: ***** len: %f\n", traceData.length);
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

//	fprintf(stderr, "(%f, %f, %f) : %f\n", x0, y0, z0, REGION_SIZE);

	while ((x0 < REGION_SIZE) && (y0 < REGION_SIZE) && (z0 < REGION_SIZE)) {
//		fprintf(stderr, "STEP: (%f, %f, %f) : %f\n", x0, y0, z0, REGION_SIZE);
		traceRay = trace(x0, y0, z0, vx, vy, vz);

		int pX0, pY0, pZ0;
		pX0 = (x0 + traceRay.x1) / 2.0 / DELTA;
		pY0 = (y0 + traceRay.y1) / 2.0 / DELTA;
		pZ0 = (z0 + traceRay.z1) / 2.0 / DELTA;

		int pX1, pY1, pZ1;
		pX1 = (traceRay.x1 - traceRay.nx / 2.0) / DELTA;
		pY1 = (traceRay.y1 - traceRay.ny / 2.0) / DELTA;
		pZ1 = (traceRay.z1 - traceRay.nz / 2.0) / DELTA;
		if ((traceRay.x1 - traceRay.nx / 2.0) < 0)
			pX1 = -1;
		if ((traceRay.y1 - traceRay.ny / 2.0) < 0)
			pY1 = -1;
		if ((traceRay.z1 - traceRay.nz / 2.0) < 0)
			pZ1 = -1;

#ifdef DEBUG
		fprintf(stderr, "DirectEikonal::traceTrajectory: traceRay.x1 - traceRay.nx / 2.0: %f. pX1: %d\n", traceRay.x1 - traceRay.nx / 2.0, pX1);
		fprintf(stderr, "DirectEikonal::traceTrajectory: traceRay.y1 - traceRay.ny / 2.0: %f. pY1: %d\n", traceRay.y1 - traceRay.ny / 2.0, pY1);
		fprintf(stderr, "DirectEikonal::traceTrajectory: traceRay.z1 - traceRay.nz / 2.0: %f. pZ1: %d\n", traceRay.z1 - traceRay.nz / 2.0, pZ1);
#endif

#ifdef DEBUG
		fprintf(stderr, "DirectEikonal::traceTrajectory: Sector0: %d, %d, %d\n", pX0, pY0, pZ0);
		fprintf(stderr, "DirectEikonal::traceTrajectory: Sector1: %d, %d, %d\n", pX1, pY1, pZ1);
#endif

		x0 = traceRay.x1;
		y0 = traceRay.y1;
		z0 = traceRay.z1;
		length += traceRay.length;
		time += traceRay.time;

#ifdef DEBUG
		fprintf(stderr, "DirectEikonal::traceTrajectory: Next point: %f, %f, %f\n", x0, y0, z0);
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
			fprintf(stderr, "DirectEikonal::traceTrajectory: Trace finished\n");
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
		fprintf(stderr, "DirectEikonal::traceTrajectory: n1: %f, n2: %f\n", n1, n2);
#endif
	}
#ifdef DEBUG
	fprintf(stderr, "DirectEikonal::traceTrajectory: Total length: %f, time: %f\n", length, time);
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

#ifdef DEBUG
	fprintf(stderr, "DirectEikonal::tracePath2: (%f, %f, %f) -> (%f, %f, %f). Dir: (%f, %f, %f)\n", x0, y0, z0, x1, y1, z1, vx0, vy0, vz0);
#endif
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
#ifdef DEBUG
		fprintf(stderr, "DirectEikonal::tracePath2: local: length: %f, time: %f\n", local.length, local.time);
#endif
		if (localDistance < distance) {
			distance = localDistance;
			result = local;
		}

		double diff_x = x1 - local.x1;
		double diff_y = y1 - local.y1;
		double diff_z = z1 - local.z1;

#ifdef DEBUG
		fprintf(stderr, "DirectEikonal::tracePath2: dest_x: %f, real_x: %f, diff_x: %f\n", x1, local.x1, diff_x);
		fprintf(stderr, "DirectEikonal::tracePath2: dest_y: %f, real_y: %f, diff_y: %f\n", y1, local.y1, diff_y);
		fprintf(stderr, "DirectEikonal::tracePath2: dest_z: %f, real_z: %f, diff_x: %f\n", z1, local.z1, diff_z);
//
		fprintf(stderr, "DirectEikonal::tracePath2: v_old x: %f, y: %f, z: %f\n", vx, vy, vz);
#endif
		vx = (x1 - x0) + SPEED * localDistance * diff_x;
		vy = (y1 - y0) + SPEED * localDistance * diff_y;
		vz = (z1 - z0) + SPEED * localDistance * diff_z;

		length = sqrt(vx * vx + vy * vy + vz * vz);
		vx = vx / length;
		vy = vy / length;
		vz = vz / length;
#ifdef DEBUG
		fprintf(stderr, "DirectEikonal::tracePath2: v_new x: %f, y: %f, z: %f\n", vx, vy, vz);

		fprintf(stderr, "DirectEikonal::tracePath2: steps: %d, localDistance: %f\n", steps, localDistance);
#endif
	}

#ifdef DEBUG
	fprintf(stderr, "DirectEikonal::tracePath2: steps: %d\n", steps);
#endif
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
	if ((fabs(vx0) < ZERO) && (fabs(vy0) < ZERO))
		fprintf(stderr, "error\n");
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

#ifdef DEBUG
	fprintf(stderr, "x0: %f, y0: %f, z0: %f\n", x0, y0, z0);
	fprintf(stderr, "x1: %f, y1: %f, z1: %f\n", x1, y1, z1);
	fprintf(stderr, "vx0: %f, vy0: %f, vz0: %f\n", vx0, vy0, vz0);
	fprintf(stderr, "vx1: %f, vy1: %f, vz1: %f\n", vx1, vy1, vz1);
	fprintf(stderr, "vx2: %f, vy2: %f, vz2: %f\n", vx2, vy2, vz2);
#endif

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

void fillPointsTest1(Point *points) {
	int p = 0;

	points[p].x = 0.5;
	points[p].y = 0;
	points[p].z = 0;
	points[p].p = 0;
	p++;

	points[p].x = 1.5;
	points[p].y = 2.0;
	points[p].z = 0;
	points[p].p = 1;
	p++;
}

void fillPoints(Point *points) {
	int p = 0;

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			double x, y, z;

			x = DELTA * i + DELTA / 2;
			y = DELTA * j + DELTA / 2;
			z = 0;
			points[p].x = x;
			points[p].y = y;
			points[p].z = z;
			points[p].p = 0;
			p++;

			x = DELTA * i + DELTA / 2;
			y = DELTA * j + DELTA / 2;
			z = SIZE * DELTA;
			points[p].x = x;
			points[p].y = y;
			points[p].z = z;
			points[p].p = 1;
			p++;

			x = DELTA * i + DELTA / 2;
			y = 0;
			z = DELTA * j + DELTA / 2;
			points[p].x = x;
			points[p].y = y;
			points[p].z = z;
			points[p].p = 2;
			p++;

			x = DELTA * i + DELTA / 2;
			y = SIZE * DELTA;
			z = DELTA * j + DELTA / 2;
			points[p].x = x;
			points[p].y = y;
			points[p].z = z;
			points[p].p = 3;
			p++;

			x = 0;
			y = DELTA * i + DELTA / 2;
			z = DELTA * j + DELTA / 2;
			points[p].x = x;
			points[p].y = y;
			points[p].z = z;
			points[p].p = 4;
			p++;

			x = SIZE * DELTA;
			y = DELTA * i + DELTA / 2;
			z = DELTA * j + DELTA / 2;
			points[p].x = x;
			points[p].y = y;
			points[p].z = z;
			points[p].p = 5;
			p++;
		}
	}
}

int main(int argc, char** argv) {
	fprintf(stderr, "Running eikonal solutioner...\n");

	int pointsCount = 6 * SIZE * SIZE;
	//int pointsCount = 2;
	Point *points = new Point[pointsCount];
	fillPoints(points);
	//fillPointsTest1(points);

#ifdef DEBUG2
	for (int i = 0; i < 6 * SIZE * SIZE; i++) {
		fprintf(stderr, "(%f, %f, %f)\n", points[i].x, points[i].y, points[i].z);
	}
#endif

	DirectEikonal directProblem;
	double length = 0, time = 0;
	int cnt = 0;

#ifdef DEBUG
	double tl = 0;
#endif
	for (int i = 0; i < pointsCount; i++) {
		for (int j = i + 1; j < pointsCount; j++) {
			Point p1, p2;
			p1 = points[i];
			p2 = points[j];

			if (p1.p == p2.p)
				continue;

#ifdef DEBUG
			double local_len = sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
			tl += local_len;
#endif
			DirectProblemData d = directProblem.tracePath2(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
			length += d.length;
			time += d.time;

#ifdef DEBUG
			if (fabs(local_len - d.length) > ZERO) {
				fprintf(stderr, "ERROR: (%d, %d) : (%f, %f, %f)\n", i, j, d.length, d.time, local_len);
				exit(0);
			}
#endif
			cnt++;
		}
	}

	fprintf(stderr, "length: %f, time: %f, cnt: %d\n", length, time, cnt);

#ifdef DEBUG
	fprintf(stderr, "test len: %f\n", tl);
#endif

/*

//	DirectEikonal directProblemSrc, directProblemTest;
	fprintf(stderr, "Running direct problem\n");
	DirectEikonal *directProblem2 = new DirectEikonal();
	fprintf(stderr, "SizeOf: %lu\n", sizeof(directProblem));
//	0.480000, 0.350000
	delete directProblem2;
	DirectProblemData d = directProblem.traceTrajectory(0.450000 * SIZE, 0.350000 * SIZE, 0.0, 0.0, 0.0, 1.0);
	fprintf(stderr, "x0: %f, y0: %f, z0: %f\n", d.x0, d.y0, d.z0);
	fprintf(stderr, "x1: %f, y1: %f, z1: %f\n", d.x1, d.y1, d.z1);
	fprintf(stderr, "length: %f, time: %f\n", d.length, d.time);

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
*/
	delete [] points;
	fprintf(stderr, "Finish!\n");
	return 0;
}
