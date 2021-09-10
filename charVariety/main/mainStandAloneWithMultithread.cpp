// X[x, y, z] = {x, z, xz - y}; Y[x, y, z] = {z, y, yz - x}
#include <cstdio>
#include <cmath>
#include <time.h>
#include <omp.h>

#define N 10
#define Pi 3.14159
#define sqrttwo 0.70710678118
//#define iterate 1
#define eps 0.000001
#define coMax 30
#define r 0.0001
#define rr 0.01

int cMax = 0;
int bad_coMax = 0;
double ab(double x){
	if (x > 0) return x;
	else return (-x);
}
double ext(double x){
	return (sqrt(x/2 + 1) + sqrt(x/2 - 1)) * sqrttwo;
}
bool compmax(double a, double b, double c){
	double aa = a, bb = b, cc = c;
	if (a < 0) aa = -a;
	if (b < 0) bb = -b;
	if (c < 0) cc = -c;
	if (aa >= bb && aa >= cc) return true;
	else return false;
}

struct pair{
	double x; double y; double z; double v1; double v2; double v3;
};

pair DTX(pair pr){
	pair prpr;
	prpr.x = pr.x; 			prpr.y = pr.z;	 		prpr.z = (pr.x)*(pr.z) - (pr.y);
	prpr.v1 = pr.v1; 		prpr.v2 = pr.v3; 		prpr.v3 = (pr.z)*(pr.v1) - pr.v2 + (pr.x)*(pr.v3);
	return prpr;
}

pair DTY(pair pr){
	pair prpr;
	prpr.x = pr.z; 			prpr.y = pr.y; 			prpr.z = (pr.y)*(pr.z) - (pr.x);
	prpr.v1 = pr.v3; 		prpr.v2 = pr.v2; 		prpr.v3 = -pr.v1 + (pr.z)*(pr.v2) + (pr.y)*(pr.v3);
	return prpr;
}

pair DTX1(pair pr){
	pair prpr;
	prpr.x = pr.x; 			prpr.y = (pr.x)*(pr.y) - (pr.z);	 					prpr.z = pr.y;
	prpr.v1 = pr.v1; 		prpr.v2 = (pr.y)*(pr.v1) + (pr.x)*(pr.v2) - pr.v3; 		prpr.v3 = pr.v2;
	return prpr;
}

pair DTY1(pair pr){
	pair prpr;
	prpr.x = (pr.x)*(pr.y) - (pr.z); 						prpr.y = pr.y; 			prpr.z = pr.x;
	prpr.v1 = (pr.y)*(pr.v1) + (pr.x)*(pr.v2) - pr.v3; 		prpr.v2 = pr.v2; 		prpr.v3 = pr.v1;
	return prpr;
}

double sqr(double a){
	return a*a;
}

double minn(double a, double b){
	return (a>b?b:a);
}

// Basic 2x2 matrix operations
struct matrix{
	double p11, p12, p21, p22;
};
matrix rotation(double t){
	matrix ans = {cos(t), -sin(t), sin(t), cos(t)};
	return ans;
}
matrix diagonal(double lambda){
	matrix ans = {lambda, 0, 0, 1/lambda};
	return ans;
}
matrix operator*(matrix const &A, matrix const &B){
	matrix C;
	C.p11 = A.p11 * B.p11 + A.p12 * B.p21;
	C.p12 = A.p11 * B.p12 + A.p12 * B.p22;
	C.p21 = A.p21 * B.p11 + A.p22 * B.p21;
	C.p22 = A.p21 * B.p12 + A.p22 * B.p22;
	return C;
}
double min_angle(matrix A){
	// return the contracting direction (angle) of A in the range [0, pi]
	double m = 0.5 * (sqr(A.p11) + sqr(A.p21) - sqr(A.p12) - sqr(A.p22));
	double n = A.p11 * A.p12 + A.p21 * A.p22;
	
	double angle;
	if (ab(m) < eps) 
		angle = Pi/2;
	else 
		angle = atan(n/m);
	
	if (m < 0) 
		angle += Pi;
	
	return fmod(0.5 * (angle + Pi), Pi);
}
double coefficient1(matrix A){
	return 0.5 * (sqr(A.p11) + sqr(A.p21) - sqr(A.p12) - sqr(A.p22));
}
double coefficient2(matrix A){
	return (A.p11 * A.p12 + A.p21 * A.p22);
}

double expansion_sq(matrix A, double angle){
	return sqr(A.p11 * cos(angle) + A.p12 * sin(angle)) + sqr(A.p21 * cos(angle) + A.p22 * sin(angle));
}
double avg_expansion(matrix* A, int matNum, double angle){
	double ans = 0.0;
	for (int i = 0; i < matNum; i++){
		ans += 0.5 * log(expansion_sq(A[i], angle));
	}
	return ans / double(matNum);
}

// Methods to find minimum of average expansion
double min_avg_expansion_newton(matrix* A, int matNum){
	
	double miin = 100000.0, ans;
	double angle_new, angle_curr;
	int co;
	
	double der_total = 0.0, dder_total = 0.0;
	double exp_sq, der_exp_sq, dder_exp_sq;
	
	double coeff0[100], coeff1[100], coeff2[100];
	
	for (int i = 0; i < matNum; i++){
		coeff0[i] = 0.5 * (sqr(A[i].p11) + sqr(A[i].p12) + sqr(A[i].p21) + sqr(A[i].p22));
		coeff1[i] = coefficient1(A[i]);
		coeff2[i] = coefficient2(A[i]);
	}
	
	for (int j = 0; j < matNum; j++){
		angle_new = min_angle(A[j]);
		angle_curr = angle_new - 1;
		co = 0;
	
		while (ab(angle_new - angle_curr) > eps && co < coMax){
			angle_curr = angle_new;
			der_total = 0.0; 
			dder_total = 0.0;
						
			for (int i = 0; i < matNum; i++){
				exp_sq = coeff0[i] + coeff1[i] * cos(2 * angle_curr) + coeff2[i] * sin(2 * angle_curr);
				der_exp_sq = 2 * (-coeff1[i] * sin(2 * angle_curr) + coeff2[i] * cos(2 * angle_curr));
				dder_exp_sq = 4 * (-coeff1[i] * cos(2 * angle_curr) - coeff2[i] * sin(2 * angle_curr));
				
				der_total += 0.5 * der_exp_sq / exp_sq;
				dder_total += 0.5 * (dder_exp_sq / exp_sq - sqr(der_exp_sq / exp_sq));
			}
		
			angle_new = angle_curr - der_total / dder_total;
			co++;
		}
		cMax = (cMax < co ? co:cMax);
		if (co == coMax) bad_coMax++;
		miin = minn(miin, avg_expansion(A, matNum, angle_new));
	}
	return miin;
}

int dum;
int n = 16, l[] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5}, count, co = 0, badco = 0, tot = 0;
int code[][5] = {
{0, 0, 0, 0, 1}, 
{0, 0, 0, 1, 1}, 
{0, 0, 1, 1, 1}, 
{0, 1, 1, 1, 1}, 
{1, 0, 0, 0, 0}, 
{1, 1, 0, 0, 0}, 
{1, 1, 1, 0, 0}, 
{1, 1, 1, 1, 0}, 
{2, 3, 3, 3, 3}, 
{2, 2, 3, 3, 3}, 
{2, 2, 2, 3, 3},
{2, 2, 2, 2, 3}, 
{3, 2, 2, 2, 2}, 
{3, 3, 2, 2, 2}, 
{3, 3, 3, 2, 2}, 
{3, 3, 3, 3, 2}};

FILE *ofp = fopen("character_variety_test.txt", "a");

//Apply the map (id) on aa = (x, y, z, vv1, vv2, vv3)
pair action(pair aa, int id){
	pair a = aa;
	for (int j = l[id] - 1; j >= 0; j--){
		if (code[id][j] == 0) a = DTX(a); 
		else if (code[id][j] == 1) a = DTY(a); 
		else if (code[id][j] == 2) a = DTX1(a);
		else if (code[id][j] == 3) a = DTY1(a);
	}
	return a;
}


int main(){
	double k = 3.99, c = 10000000000000;
	double ans_x, ans_y, ans_z, rt[10], cand, tem;
	double t1, t2, j;
	int ii, i_bound = 4.0/r, thread_id;
	matrix act[50];

	
	time_t tttt;
	time(&tttt);
	printf("Program started on: %s", ctime(&tttt));
	
//	clock_t time_taken;
	
//	time_taken = clock();

	#pragma omp parallel private(thread_id, t1, t2) num_threads(2)
	{
		t1 = omp_get_wtime();
		#pragma omp for private(j, rt, cand, tem, act) reduction(min: c) reduction(+: bad_coMax) reduction(+: tot) nowait
		for (ii = 0; ii < i_bound; ii++){
			double i = double(ii)*r - 2;
			for (j = -2; j < 2; j += r){
				double st = i*i*j*j - 4*(i*i + j*j - k);
				if (st > 0.000001){
					st = sqrt(st);
					rt[0] = (i*j + st) * 0.5;
					rt[1] = (i*j - st) * 0.5;
					for (int s = 0; s < 2; s++){
						double A = 2*i - j*rt[s], B = 2*j - i*rt[s], C = 2*rt[s] - i*j;
						if (compmax(C, A, B)){
							double x = i, 		y = j, 		z = rt[s];
							double vv1 = 0, 	vv2 = C, 	vv3 = -B;
							double ww1 = -C, 	ww2 = 0, 	ww3 = A;
							for (int ss = 0; ss < 3; ss++){	
								bool bad_count = false;		
								
								for (int map_id = 0; map_id < n; map_id++){
									
									pair pr, ps;
									pr.x = x; pr.y = y; pr.z = z; pr.v1 = vv1; pr.v2 = vv2; pr.v3 = vv3; 
									ps.x = x; ps.y = y; ps.z = z; ps.v1 = ww1; ps.v2 = ww2; ps.v3 = ww3;
									
									pr = action(pr, map_id);
									ps = action(ps, map_id);
																	
									double AA = 2*(pr.x) - (pr.y)*(pr.z), BB = 2*(pr.y) - (pr.x)*(pr.z), CC = 2*(pr.z) - (pr.x)*(pr.y);
									
									if (compmax(CC, AA, BB)){
										double sst = 1 / sqrt(ab(C * CC));
										act[map_id].p11 = pr.v2 * sst;
										act[map_id].p21 = - pr.v1 * sst;
										act[map_id].p12 = ps.v2 * sst;
										act[map_id].p22 = - ps.v1 * sst;
									} 
									else if (compmax(BB, AA, CC)){
										double sst = 1 / sqrt(ab(B * BB));
										act[map_id].p11 = pr.v1 * sst;
										act[map_id].p21 = - pr.v3 * sst;
										act[map_id].p12 = ps.v1 * sst;
										act[map_id].p22 = - ps.v3 * sst;
									} 
									else {
										double sst = 1 / sqrt(ab(A * AA));
										act[map_id].p11 = pr.v3 * sst;
										act[map_id].p21 = - pr.v2 * sst;
										act[map_id].p12 = ps.v3 * sst;
										act[map_id].p22 = - ps.v2 * sst;
									}
								}
								
								cand = min_avg_expansion_newton(act, n);
								if (c > cand){
									c = cand;
									ans_x = x;
									ans_y = y;
									ans_z = z;
								}
								tot++;
	
								tem = z; z = y; y = x; x = tem ;
								tem = vv3; vv3 = vv2; vv2 = vv1; vv1 = tem;
								tem = ww3; ww3 = ww2; ww2 = ww1; ww1 = tem;	
							}
						}
					}
				}
			}
		}
		t2 = omp_get_wtime();
		thread_id = omp_get_thread_num();
		printf("Thread: %d, Time taken: %f\n", thread_id, t2 - t1);
		
	}
	
//	time_taken = clock() - time_taken;
//	double time_taken_sec = (double)time_taken / CLOCKS_PER_SEC;
	
	time_t ttt;
	time(&ttt);
	
	fprintf(ofp, "Newton\n");
	fprintf(ofp, "Program run on: %s", ctime(&ttt));
	fprintf(ofp, "k-2: %f c: %f r: %f rr: %f tot: %d\n", k-2, c, r, rr, tot);
	fprintf(ofp, "P = (%f, %f, %f)\n", ans_x, ans_y, ans_z);
//	fprintf(ofp, "Time taken: %f seconds\n", time_taken_sec);
	fprintf(ofp, "cMax = %d, bad_coMax = %d\n", cMax, bad_coMax);
	fprintf(ofp, "\n\n\n");
	return 0;
}
