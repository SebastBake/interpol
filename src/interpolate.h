/***************************************************************************
 * Written by Sebastian Baker for Numerical Algorithms at Unimelb
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "thomas_alg.h"

#ifndef INTERPOLATE_H

#define INTERP_INIT_ARRLEN 20
#define QUADRATIC_NUM_PTS 3
#define CUB_SPLINE_COMPUTE_SUCCESS 1
#define CUB_SPLINE_COMPUTE_FAIL -1
#define TINY(x) fabs(x) < 1e-11
#define CUB_SPLINE_B(h_1, a_1, a_2, c_1, c_2) (a_2-a_1)/h_1 - h_1*(2.0*c_1+c_2)/3.0
#define CUB_SPLINE_C_RHS(h_0, h_1, a_0, a_1, a_2) 3.0*(a_2-a_1)/h_1 + 3.0*(a_0-a_1)/h_0
#define CUB_SPLINE_C_a(h_1, h_0) 2.0*(h_1+h_0)
#define CUB_SPLINE_D(h_1, c_1, c_2) (c_2-c_1)/(3.0*h_1)
#define EVAL_CUB_SPLINE(a,b,c,d,x_i,x) a + b*(x-x_i) + c*pow(x-x_i,2) + d*pow(x-x_i,3)

// simple struct to hold an x,fx tuple
typedef struct interp_pt {

	double x;
	double fx;

} interp_pt_t;

// dynamic array to hold a set of points to interpolate
typedef struct interp_set {

	interp_pt_t** pts;
	int N;
	int arrLen;

} interp_set_t;

// represents a lagrange term
typedef struct lagrange_term {

	int index;
	int order;
	double* root;
	double fx_i;
	double denominator;

} lagrange_term_t;

// represents a lagrange equation
typedef struct lagrange_eqn {

	lagrange_term_t** terms;
	int num_terms;

} lagrange_eqn_t;

// represents a cubic spline segment
typedef struct cub_spline_segment {

	int index;
	double x_lo;
	double a;
	double b;
	double c;
	double d;

} cub_spline_seg_t;

// represents a cubic spline
typedef struct cub_spline {

	cub_spline_seg_t** segs;
	int num_segs; // not including the unused/incomplete end segment eg: 4 points means 3 segments

} cub_spline_t;

// manage set dynamic array
interp_set_t* newInterpSet();
interp_pt_t* newInterpPt(double x, double fx);
void appendPtToSet(interp_set_t* set, interp_pt_t* pt);
void freeInterpSet(interp_set_t* set);
void freeInterpPt(interp_pt_t* pt);

// create, evaluate lagrange equations
lagrange_eqn_t* newLagrangeEqn(interp_set_t* set);
void freeLagrangeEqn(lagrange_eqn_t* eqn);
interp_pt_t* evaluateLagrangeEqn(lagrange_eqn_t* eqn, double x);

// create, evaluate cubic splines
cub_spline_t* newCubSpline(interp_set_t* set);
void freeCubSpline(cub_spline_t* spline);
interp_pt_t* evaluateCubSpline(cub_spline_t* spline, double x);

#endif