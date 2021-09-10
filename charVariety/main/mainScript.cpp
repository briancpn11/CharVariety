// charVariety/main/mainScript.cpp

// X[x, y, z] = {x, z, xz - y}; Y[x, y, z] = {z, y, yz - x}

#include "../utils/helpers.h" 
#include <cmath>
#include <time.h>
#include <omp.h>

#define r 0.01
#define rr 0.01

int nMap = 16, nIter[] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
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
    int tot = 0;
    double k = 3.99, c = 10000000000000;
    double ans_x, ans_y, ans_z, m_z[10], cand, tem;
    Matrix act[50];

    time_t timeNow;
    time(&timeNow);
    printf("Program started on: %s", ctime(&timeNow));

    clock_t startTime;
    startTime = clock();

    int thread_id, ii, i_bound=4.0/r; double t1, t2;
    #pragma omp parallel private(thread_id, t1, t2) num_threads(1)
    {
        t1 = omp_get_wtime();
        #pragma omp for private(m_z, cand, tem, act) reduction(min: c) reduction(+: tot) nowait
        for (ii = 0; ii < i_bound; ii++){
            double m_x = double(ii)*r - 2;
        // for (double m_x = -2; m_x < 2; m_x += r){            // swapped from this to facilitate the omp syntax
            for (double m_y = -2; m_y < 2; m_y += r){
                double disc = sqr(m_x*m_y) - 4*(sqr(m_x) + sqr(m_y) - k);
                if (disc > EPS){
                    disc = sqrt(disc);
                    m_z[0] = (m_x*m_y + disc) * 0.5;
                    m_z[1] = (m_x*m_y - disc) * 0.5;
                    for (int s = 0; s < 2; s++){
                        double A = 2*m_x - m_y*m_z[s], B = 2*m_y - m_x*m_z[s], C = 2*m_z[s] - m_x*m_y;
                        if (compMax(C, A, B)){
                            double x = m_x, 	y = m_y, 	z = m_z[s];
                            double vv1 = 0, 	vv2 = C, 	vv3 = -B;
                            double ww1 = -C, 	ww2 = 0, 	ww3 = A;
                            for (int ss = 0; ss < 3; ss++){						
                                for (int map_id = 0; map_id < nMap; map_id++){
                                    Pair pr(x, y, z, vv1, vv2, vv3), ps(x, y, z, ww1, ww2, ww3);
                                    pr = action(pr, maps[map_id], nIter[map_id]);
                                    ps = action(ps, maps[map_id], nIter[map_id]);
                                                                    
                                    double AA = 2*(pr.x) - (pr.y)*(pr.z), BB = 2*(pr.y) - (pr.x)*(pr.z), CC = 2*(pr.z) - (pr.x)*(pr.y);
                                    if (compMax(CC, AA, BB)){
                                        double sst = 1 / sqrt(absVal(C * CC));
                                        act[map_id] = {pr.v2 * sst, ps.v2 * sst, -pr.v1 * sst, -ps.v1 * sst};
                                    } 
                                    else if (compMax(BB, AA, CC)){
                                        double sst = 1 / sqrt(absVal(B * BB));
                                        act[map_id] = {pr.v1 * sst, ps.v1 * sst, -pr.v3 * sst, -ps.v3 * sst};
                                    } 
                                    else {
                                        double sst = 1 / sqrt(absVal(A * AA));
                                        act[map_id] = {pr.v3 * sst, ps.v3 * sst, -pr.v2 * sst, -ps.v2 * sst};
                                    }
                                }
                                cand = minAvgExpansionNewton(act, nMap);
                                if (c > cand){
                                    c = cand; ans_x = x; ans_y = y; ans_z = z;
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
        t2 = omp_get_wtime();
        thread_id = omp_get_thread_num();
        printf("Thread: %d, Time taken: %.3f sec\n", thread_id, t2 - t1);
    }
	
    double timeTaken = clock() - startTime;
    double time_taken_sec = (double)timeTaken / CLOCKS_PER_SEC;

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
