#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>

const int maxq = 1<<16;

double failure_rate(int q, int p, int n, int m, int g, int reclen, double noise[], int noise_width) {
    int d = 0;

    // find the largest integer d such that 2md < p(1 - 1/g)
    while (m*(d+1)*2. < p*(1. - 1./g)) d++;
    printf("d = %d\n", d);

    // Let F[i][a][c] = Pr[{A X_1}_p^T X_2 - epsilon^T X_2 = c, (A X_1)^T X_2 mod (q/p) = a], where A is a i x i matrix
    // Then F[i+1][*][*] can be derived from F[i][*][*].
    // Hence, we only need to store F[i][*][*] in memory when compute F[i+1][*][*].
    // In specific, we store F[i][*][*] in F[i mod 2][*][*] to save memory cost during the computation.
    // The result is stored in F[n mod 2][*][*]
    static double F[2][1<<5][maxq];
    memset(F, 0, sizeof(F));
    F[0][0][0] = 1.;

    for (int i = 0;i < n;i++) {
        memset(F[(i+1) & 1], 0, sizeof(F[(i+1) & 1]));
        for (int a = 0;a < q/p;a++) {
            for (int c = 0;c < q;c++) {
                // let e1 be the n-th element of the vector {A X_1}_p, 
                // x2 be the n-th element of the vector X_2, 
                // and eps be the n-th element of the vector epsilon.
                // Then f[i+1][a + e1*x2][c + (e1 - eps)*x2] += f[i][a][c] * Pr[e1, x2, eps], 
                // where Pr[e1, x2, eps] is the joint distribution of (e1, x2, eps)
                if (F[i & 1][a][c] == 0.) continue ;
                for (int e1 = -q/(2*p);e1 < q/(2*p);e1++) {
                    for (int x2 = -noise_width;x2 <= noise_width;x2++) {
                        for (int eps = -q/(2*p);eps < q/(2*p);eps++) {
                            F[(i+1) & 1][(a + e1*x2 + q) % (q/p)][(c + (e1 - eps)*x2 + q) % q] += 
                                F[i & 1][a][c] * noise[x2 + noise_width] / (q/p * q/p);
                        }
                    }
                }
            }
        }
    }

    static double G[2][maxq];
    memset(G, 0, sizeof(G));
    G[0][0] = 1.;

    // Let G[i][c] = Pr[X_1^T {A^T X_2}_p = c], where A is a i x i matrix.
    // Then G[i+1][*] can be derived from G[i][*].
    // Hence, we only need to store G[i][*][*] in memory when compute G[i+1][*][*].
    // In specific, we store G[i][*] in G[i mod 2][*] to save memory cost during the computation.
    // The result is stored in G[n mod 2][*]
    for (int i = 0;i < n;i++) {
        memset(G[(i+1) & 1], 0, sizeof(G[(i + 1) & 1]));
        for (int c = 0;c < q;c++) {
            // Let x1 be the n-th element of the vector X_1, 
            // and e2 be the n-th element of the vector {A^T X_2}_p
            // Then, G[i+1][c + x1*e2] += G[i][c] * Pr[x1, e2]
            for (int e2 = -q/(2*p);e2 < q/(2*p);e2++) {
                for (int x1 = -noise_width;x1 <= noise_width;x1++) {
                    G[(i+1) & 1][(c + x1*e2 + q) % q] += G[i & 1][c] * noise[x1 + noise_width] / (q/p);
                }
            }
        }
    }

    // See the formulas in P14 of the paper for the detail
    static double distr[maxq];
    memset(distr, 0, sizeof(distr));
    for (int a = 0;a < q/p;a++) {
        double pr_a = 1./(q/p);
        for (int c1 = 0;c1 < q;c1++) if (c1 % (q/p) == a) {
            for (int c2 = 0;c2 < q;c2++) {
                distr[(c1 - c2 + q) % q] += G[n & 1][c1] * F[n & 1][a][c2] / pr_a;
            }
        }
    }

    // derive the distribution of \Sigma_2 - \Sigma_1 in page 12
    static double distr_p[maxq];
    memset(distr_p, 0, sizeof(distr_p));
    for (int c = 0;c < q;c++) {
        distr_p[ (int)(1.*p/q*c + .5) % p ] += distr[c];
    }

    double failure_pr = 0.;
    for (int x = 0;x < p;x++) {
        if (std::min(x, p - x) > d) {
            failure_pr += distr_p[x];
        }
    }

    // apply the union bound
    failure_pr *= reclen;

    return failure_pr;
}

double failure_rate_simple(int q, int p, int n, int m, int g, int reclen, double noise[], int noise_width) {
    int d = 0;

    // find the largest integer d such that 2md < p(1 - 1/g)
    while (m+2*(d + 1) < g) d++;
    printf("d = %d\n", d);

    // Let F[i][a][c] = Pr[{A X_1}_p^T X_2 - epsilon^T X_2 = c, (A X_1)^T X_2 mod (q/p) = a], where A is a i x i matrix
    // Then F[i+1][*][*] can be derived from F[i][*][*].
    // Hence, we only need to store F[i][*][*] in memory when compute F[i+1][*][*].
    // In specific, we store F[i][*][*] in F[i mod 2][*][*] to save memory cost during the computation.
    // The result is stored in F[n mod 2][*][*]
    static double F[2][1<<5][maxq];
    memset(F, 0, sizeof(F));
    F[0][0][0] = 1.;

    for (int i = 0;i < n;i++) {
        memset(F[(i+1) & 1], 0, sizeof(F[(i+1) & 1]));
        for (int a = 0;a < q/p;a++) {
            for (int c = 0;c < q;c++) {
                // let e1 be the n-th element of the vector {A X_1}_p, 
                // x2 be the n-th element of the vector X_2, 
                // and eps be the n-th element of the vector epsilon.
                // Then f[i+1][a + e1*x2][c + (e1 - eps)*x2] += f[i][a][c] * Pr[e1, x2, eps], 
                // where Pr[e1, x2, eps] is the joint distribution of (e1, x2, eps)
                if (F[i & 1][a][c] == 0.) continue ;
                for (int e1 = -q/(2*p);e1 < q/(2*p);e1++) {
                    for (int x2 = -noise_width;x2 <= noise_width;x2++) {
                        for (int eps = -q/(2*p);eps < q/(2*p);eps++) {
                            F[(i+1) & 1][(a + e1*x2 + q) % (q/p)][(c + (e1 - eps)*x2 + q) % q] += 
                                F[i & 1][a][c] * noise[x2 + noise_width] / (q/p * q/p);
                        }
                    }
                }
            }
        }
    }

    static double G[2][maxq];
    memset(G, 0, sizeof(G));
    G[0][0] = 1.;

    // Let G[i][c] = Pr[X_1^T {A^T X_2}_p = c], where A is a i x i matrix.
    // Then G[i+1][*] can be derived from G[i][*].
    // Hence, we only need to store G[i][*][*] in memory when compute G[i+1][*][*].
    // In specific, we store G[i][*] in G[i mod 2][*] to save memory cost during the computation.
    // The result is stored in G[n mod 2][*]
    for (int i = 0;i < n;i++) {
        memset(G[(i+1) & 1], 0, sizeof(G[(i + 1) & 1]));
        for (int c = 0;c < q;c++) {
            // Let x1 be the n-th element of the vector X_1, 
            // and e2 be the n-th element of the vector {A^T X_2}_p
            // Then, G[i+1][c + x1*e2] += G[i][c] * Pr[x1, e2]
            for (int e2 = -q/(2*p);e2 < q/(2*p);e2++) {
                for (int x1 = -noise_width;x1 <= noise_width;x1++) {
                    G[(i+1) & 1][(c + x1*e2 + q) % q] += G[i & 1][c] * noise[x1 + noise_width] / (q/p);
                }
            }
        }
    }

    // See the formulas in P14 of the paper for the detail
    static double distr[maxq];
    memset(distr, 0, sizeof(distr));
    for (int a = 0;a < q/p;a++) {
        double pr_a = 1./(q/p);
        for (int c1 = 0;c1 < q;c1++) if (c1 % (q/p) == a) {
            for (int c2 = 0;c2 < q;c2++) {
                distr[(c1 - c2 + q) % q] += G[n & 1][c1] * F[n & 1][a][c2] / pr_a;
            }
        }
    }

    // derive the distribution of \Sigma_2 - \Sigma_1 in page 12
    static double distr_p[maxq];
    memset(distr_p, 0, sizeof(distr_p));
    for (int c = 0;c < q;c++) {
        distr_p[ (int)(1.*p/q*c + .5) % p ] += distr[c];
    }

    double failure_pr = 0.;
    for (int x = 0;x < p;x++) {
        if (std::min(x, p - x) > d) {
            failure_pr += distr_p[x];
        }
    }

    // apply the union bound
    failure_pr *= reclen;

    return failure_pr;
}


double bandwidth(int n, int q, int p, int g, int l) {
    return (2*(int)(ceil(log2(p)))*n*l + l*l*(int)(ceil(log2(g))))*1./8000.;
}

void LWR_recommended()
{
    int n = 680, q = 1<<15, p = 1<<12, g = 1<<8, l = 8, m = 1<<4, reclen = l*l;
    double noise[] = {3./65536, 44./65536, 389./65536, 2090./65536, 
                     6938./65536, 14249./65536, 18110./65536, 14249./65536, 
                     6938./65536, 2090./65536, 389./65536, 44./65536, 3./65536}; 
//    double noise[] = { 
//         17./65536, 220./65536, 1570./65536, 6383./65536, 14792./65536,
//         19572./65536, 14792./65536, 6383./65536, 1570./65536, 220./65536, 17./65536
//    };

    int noise_width = 6;
    double err = failure_rate(q, p, n, m, g, reclen, noise, noise_width);
    printf("error rate = %le log2 = %lf, bandwidth = %lf\n", err, log2(err), bandwidth(n, q, p, g, l));
}

void LWR_paranoid()
{
    int n = 832, q = 1<<15, p = 1<<12, g = 1<<8, l = 8, m = 1<<4, reclen = l*l;
    double noise[] = {4./65536, 97./65536, 1033./65536, 5580./65536, 
                      15326./65536, 21456./65536, 15326./65536, 5580./65536, 1033./65536, 97./65536, 4./65536};
    int noise_width = 5;
    double err = failure_rate(q, p, n, m, g, reclen, noise, noise_width);
    printf("error rate = %le log2 = %lf, bandwidth = %lf\n", err, log2(err), bandwidth(n, q, p, g, l));
}

void LWR_AKCN_recommended_simple()
{
    //int n = 680, q = 1<<15, p = 1<<12, g = 1<<8, l = 8, m = 1<<4, reclen = l*l;
    int n = 672, q = 1<<15, p = 1<<12, g = 1<<8, l = 8, m = 1<<4, reclen = l*l;
    double noise[] = {3./65536, 44./65536, 389./65536, 2090./65536, 
                     6938./65536, 14249./65536, 18110./65536, 14249./65536, 
                     6938./65536, 2090./65536, 389./65536, 44./65536, 3./65536}; 
//    double noise[] = { 
//         17./65536, 220./65536, 1570./65536, 6383./65536, 14792./65536,
//         19572./65536, 14792./65536, 6383./65536, 1570./65536, 220./65536, 17./65536
//    };

    int noise_width = 6;
    double err = failure_rate_simple(q, p, n, m, g, reclen, noise, noise_width);
    printf("error rate = %le log2 = %lf, bandwidth = %lf\n", err, log2(err), bandwidth(n, q, p, g, l));
}

void LWR_AKCN_paranoid_simple()
{
    int n = 832, q = 1<<15, p = 1<<12, g = 1<<8, l = 8, m = 1<<4, reclen = l*l;
    double noise[] = {4./65536, 97./65536, 1033./65536, 5580./65536, 
                      15326./65536, 21456./65536, 15326./65536, 5580./65536, 1033./65536, 97./65536, 4./65536};
    int noise_width = 5;
    double err = failure_rate_simple(q, p, n, m, g, reclen, noise, noise_width);
    printf("error rate = %le log2 = %lf, bandwidth = %lf\n", err, log2(err), bandwidth(n, q, p, g, l));
}

int main() {
    //LWR_recommended();
    //LWR_paranoid();

    LWR_AKCN_recommended_simple();
    LWR_AKCN_paranoid_simple();
    return 0;
}
