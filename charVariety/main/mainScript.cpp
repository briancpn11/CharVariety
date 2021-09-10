// charVariety/main/mainScript.cpp

// X[x, y, z] = {x, z, xz - y}; Y[x, y, z] = {z, y, yz - x}

#include "../utils/helpers.h" 
#include <cmath>
#include <time.h>
#include <omp.h>

#define r 0.01
#define rr 0.01

int n = 16, nIter[] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5}, tot = 0;
int maps[][5] ={{0, 0, 0, 0, 1}, 
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

// int n = 2, nIter[] = {2, 2}, tot = 0;
// int maps[][5] ={{0, 1}, {1, 0}};

FILE *ofp = fopen("result/minAvgExpansionGrid.txt", "a");

int main(){
	double k = 3.99, c = 10000000000000;
	double ans_x, ans_y, ans_z, rt[10], cand, tem;
	int i_bound = 4.0/r, thread_id;
	Matrix act[50];

	
	time_t tttt;
	time(&tttt);
	printf("Program started on: %s", ctime(&tttt));
	
	clock_t time_taken;
	
	time_taken = clock();

	for (int ii = 0; ii < i_bound; ii++){
		double i = double(ii)*r - 2;
		for (double j = -2; j < 2; j += r){
			double st = i*i*j*j - 4*(i*i + j*j - k);
			if (st > 0.000001){
				st = sqrt(st);
				rt[0] = (i*j + st) * 0.5;
				rt[1] = (i*j - st) * 0.5;
				for (int s = 0; s < 2; s++){
					double A = 2*i - j*rt[s], B = 2*j - i*rt[s], C = 2*rt[s] - i*j;
					if (compMax(C, A, B)){
						double x = i, 		y = j, 		z = rt[s];
						double vv1 = 0, 	vv2 = C, 	vv3 = -B;
						double ww1 = -C, 	ww2 = 0, 	ww3 = A;
						for (int ss = 0; ss < 3; ss++){	
							bool bad_count = false;		
							
							for (int map_id = 0; map_id < n; map_id++){
								
								Pair pr, ps;
								pr.x = x; pr.y = y; pr.z = z; pr.v1 = vv1; pr.v2 = vv2; pr.v3 = vv3; 
								ps.x = x; ps.y = y; ps.z = z; ps.v1 = ww1; ps.v2 = ww2; ps.v3 = ww3;
								
								pr = action(pr, maps[map_id], nIter[map_id]);
								ps = action(ps, maps[map_id], nIter[map_id]);
																
								double AA = 2*(pr.x) - (pr.y)*(pr.z), BB = 2*(pr.y) - (pr.x)*(pr.z), CC = 2*(pr.z) - (pr.x)*(pr.y);
								
								if (compMax(CC, AA, BB)){
									double sst = 1 / sqrt(abs(C * CC));
									act[map_id].p11 = pr.v2 * sst;
									act[map_id].p21 = - pr.v1 * sst;
									act[map_id].p12 = ps.v2 * sst;
									act[map_id].p22 = - ps.v1 * sst;
								} 
								else if (compMax(BB, AA, CC)){
									double sst = 1 / sqrt(abs(B * BB));
									act[map_id].p11 = pr.v1 * sst;
									act[map_id].p21 = - pr.v3 * sst;
									act[map_id].p12 = ps.v1 * sst;
									act[map_id].p22 = - ps.v3 * sst;
								} 
								else {
									double sst = 1 / sqrt(abs(A * AA));
									act[map_id].p11 = pr.v3 * sst;
									act[map_id].p21 = - pr.v2 * sst;
									act[map_id].p12 = ps.v3 * sst;
									act[map_id].p22 = - ps.v2 * sst;
								}
							}
							
							cand = minAvgExpansionNewton(act, n);
							if (c > cand){
								c = cand;
								ans_x = x;
								ans_y = y;
								ans_z = z;
							}
							tot++;
	
							tem = z; z = y; y = x; x = tem;
							tem = vv3; vv3 = vv2; vv2 = vv1; vv1 = tem;
							tem = ww3; ww3 = ww2; ww2 = ww1; ww1 = tem;	
						}
					}
				}
			}
		}
	}
	
	time_taken = clock() - time_taken;
	double time_taken_sec = (double)time_taken / CLOCKS_PER_SEC;
	
	time_t ttt;
	time(&ttt);
	
	fprintf(ofp, "Program run on: %s", ctime(&ttt));
	fprintf(ofp, "Minimum average expansion: %f\n", c);
	fprintf(ofp, "Other parameters: \n");
	fprintf(ofp, " k-2: %f\n", k-2);
	fprintf(ofp, " Spacing in the surface directions: r: %f\n", r);
	fprintf(ofp, " Spacing in the tangent direction: rr: %f\n", rr);
	fprintf(ofp, " Total number of grid points on the surface: %d\n", tot);
	fprintf(ofp, "The point where the min average expansion occur: P = (%f, %f, %f)\n", ans_x, ans_y, ans_z);
	fprintf(ofp, "Time taken: %f seconds\n", time_taken_sec);
	// fprintf(ofp, "cMax = %d, bad_coMax = %d\n", cMax, bad_coMax);
	fprintf(ofp, "\n\n\n");
	return 0;
}
