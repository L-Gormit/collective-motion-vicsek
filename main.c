#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

/*
 * Vicsek model for collective motion
 * ---------------------------------
 * Periodic domain, metric interactions, angular noise
 *
 * Outputs:
 *  - Traj.txt  : particle trajectories
 *  - phase.txt : global order parameter Φ(t)
 */

/* ---------------- Data structures ---------------- */

typedef struct {
    double pos[2];   // particle position (x,y)
    double angle;    // orientation angle θ
    int index;       // particle index
} Particle;

typedef struct {
    Particle *p;     // particle array
    double L;        // domain size
    int N;           // number of particles
} Domain;

/* ---------------- Utility functions ---------------- */

/* Periodic wrapping of coordinates */
void wrap(double *x, double L){
    if (*x >= L) *x -= L;
    if (*x < 0)  *x += L;
}

/* Periodic distance in one dimension */
double pbc_distance(double dx, double L){
    if (dx >  0.5 * L) dx -= L;
    if (dx < -0.5 * L) dx += L;
    return dx;
}

/* ---------------- Initialization ---------------- */

/* Random uniform initialization of particles */
void initialize_particle(Domain *dom){
    for (int i = 0; i < dom->N; i++){
        dom->p[i].index = i;
        dom->p[i].pos[0] = dom->L * rand() / RAND_MAX;
        dom->p[i].pos[1] = dom->L * rand() / RAND_MAX;
        dom->p[i].angle  = 2.0 * M_PI * rand() / RAND_MAX;
    }
}

/* ---------------- Measurements ---------------- */

/* Global order parameter Φ */
double order_parameter(Domain *dom)
{
    double vx = 0.0, vy = 0.0;

    for (int i = 0; i < dom->N; i++){
        vx += cos(dom->p[i].angle);
        vy += sin(dom->p[i].angle);
    }

    return sqrt(vx*vx + vy*vy) / dom->N;
}

/* ---------------- Vicsek alignment ---------------- */

/*
 * Compute mean direction of neighbors within radius R0
 * using complex representation: exp(iθ)
 */
double mean_neighbor_angle(Domain *dom, int i, double R0)
{
    double complex sum = 0.0;
    double xi = dom->p[i].pos[0];
    double yi = dom->p[i].pos[1];

    for (int j = 0; j < dom->N; j++){
        double dx = pbc_distance(dom->p[j].pos[0] - xi, dom->L);
        double dy = pbc_distance(dom->p[j].pos[1] - yi, dom->L);

        if (dx*dx + dy*dy < R0*R0){
            sum += cexp(I * dom->p[j].angle);
        }
    }

    /* If no neighbors, keep current direction */
    if (cabs(sum) < 1e-3)
        return dom->p[i].angle;

    double angle = carg(sum);
    if (angle < 0) angle += 2.0 * M_PI;

    return angle;
}

/* Update particle orientations with angular noise */
void update_angles(Domain *dom, double R0, double eta)
{
    double *new_angle = malloc(dom->N * sizeof(double));

    for (int i = 0; i < dom->N; i++){
        double theta = mean_neighbor_angle(dom, i, R0);
        double noise = eta * (rand() / (double)RAND_MAX - 0.5);

        new_angle[i] = fmod(theta + noise, 2.0 * M_PI);
        if (new_angle[i] < 0) new_angle[i] += 2.0 * M_PI;
    }

    for (int i = 0; i < dom->N; i++)
        dom->p[i].angle = new_angle[i];

    free(new_angle);
}

/* ---------------- Time integration ---------------- */

/*
 * Time loop:
 *  - move particles
 *  - apply periodic BC
 *  - compute order parameter
 *  - update orientations
 */
void motion(Domain *dom, double dt, double T,
            double v0, double R0, double eta)
{
    FILE *fp  = fopen("Traj.txt",  "w");
    FILE *fp2 = fopen("phase.txt", "w");

    for (double t = 0; t < T; t += dt){
        for (int i = 0; i < dom->N; i++){
            dom->p[i].pos[0] += dt * v0 * cos(dom->p[i].angle);
            dom->p[i].pos[1] += dt * v0 * sin(dom->p[i].angle);

            wrap(&dom->p[i].pos[0], dom->L);
            wrap(&dom->p[i].pos[1], dom->L);

            fprintf(fp, "%g %d %g %g %g\n",
                    t, i,
                    dom->p[i].pos[0],
                    dom->p[i].pos[1],
                    dom->p[i].angle);
        }

        fprintf(fp2, "%g %g\n", t, order_parameter(dom));
        update_angles(dom, R0, eta);
    }

    fclose(fp);
    fclose(fp2);
}

/* ---------------- Main ---------------- */

int main(void)
{
    srand(time(NULL));

    double L, R0, eta, dt, T, v0;
    int N;

    FILE *fptr = fopen("param.txt", "r");
    if (!fptr){
        perror("param.txt");
        return 1;
    }

    char label[64];
    fscanf(fptr, "%s %lf", label, &L);
    fscanf(fptr, "%s %lf", label, &R0);
    fscanf(fptr, "%s %lf", label, &eta);
    fscanf(fptr, "%s %d",  label, &N);
    fscanf(fptr, "%s %lf", label, &dt);
    fscanf(fptr, "%s %lf", label, &T);
    fscanf(fptr, "%s %lf", label, &v0);
    fclose(fptr);

    printf("L=%.2f R0=%.2f eta=%.2f N=%d dt=%.2f T=%.2f v0=%.2f\n",
           L, R0, eta, N, dt, T, v0);

    Domain dom;
    dom.L = L;
    dom.N = N;
    dom.p = malloc(N * sizeof(Particle));

    initialize_particle(&dom);
    motion(&dom, dt, T, v0, R0, eta);

    free(dom.p);
    return 0;
}

