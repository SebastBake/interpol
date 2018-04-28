/***************************************************************************
 * Written by Sebastian Baker for Numerical Algorithms at Unimelb
 ***************************************************************************/

#include "thomas_alg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int computeTridiag_a_s(tridiag_t *m);
int computeTridiag_Q_s(tridiag_t *m);
int computeTrigiag_x(tridiag_t *m);

// create a new tridiagonal matrix linear system structure
tridiag_t *newTridiag() {
	
	tridiag_t *m = (tridiag_t*)malloc(sizeof(tridiag_t));
	assert(m != NULL);

	m->arrLen = THOMAS_INIT_ARR_LEN;
	m->N = 0;

	m->rows = (tridiag_row_t*)calloc( m->arrLen, sizeof(tridiag_row_t));
	assert(m->rows != NULL);

	return m;
}

// free a tridiagonal matrix structure
void freeTridiag(tridiag_t *m) {
	assert(m!=NULL);
	free(m->rows);
	free(m);
}

// add a row onto the bottom of a tridiagonal matrix
void appendTridiagRow(tridiag_t *m, double a, double b, double c, double Q) {

	assert(m!=NULL);	

	m->N++;
	if(m->N == m->arrLen) {
		m->arrLen += THOMAS_INIT_ARR_LEN;
		m->rows = (tridiag_row_t*)realloc( m->rows, m->arrLen*sizeof(tridiag_row_t));
		assert(m->rows!=NULL);
	}

	tridiag_row_t *r = getTridiagRow(m, m->N);
	r->a = a;
	r->b = b;
	r->c = c;
	r->Q = Q;
}

// get a row of a tridiagonal matrix
tridiag_row_t *getTridiagRow(tridiag_t *m, int i) {
	assert(m!=NULL);
	assert(i<=m->N);
	return &(m->rows[i-1]);
}

// Use thomas algorithm to solve a tridiagonal linear system
int solveTridiag(tridiag_t *m) {

	assert(m!=NULL);
	
	if(
		computeTridiag_a_s(m) == SOLVER_SUCCESS &&
		computeTridiag_Q_s(m) == SOLVER_SUCCESS &&
		computeTrigiag_x(m) == SOLVER_SUCCESS
	) {
		return SOLVER_SUCCESS;
	}

	return SOLVER_FAIL;
}

// compute the a_s values for the tridiagonal linear system
int computeTridiag_a_s(tridiag_t *m) {

	assert(m!=NULL);

	// case for row 1 : a_s = a
	int i=COMPUTE_A_S_START;
	tridiag_row_t *r = getTridiagRow(m, i);	
	tridiag_row_t *rPrev = r;
	r->a_s = r->a;

	// case for row 2,3,4...N : a_s = a − c * bPrev / a_sPrev
	while(++i <= m->N) {
		rPrev = r;
		r = getTridiagRow(m, i);
		if(TINY(rPrev->a_s)) { return SOLVER_FAIL; }
		r->a_s = r->a - r->c * rPrev->b / rPrev->a_s;
	}

	return SOLVER_SUCCESS;
}

// compute the Q_s values for the tridiagonal linear system
int computeTridiag_Q_s(tridiag_t *m) {
	
	assert(m!=NULL);

	// case for row 1 : Q_s = Q
	int i=COMPUTE_Q_S_START;
	tridiag_row_t *r = getTridiagRow(m, i);	
	tridiag_row_t *rPrev = r;
	r->Q_s = r->Q;

	// case for row 2,3,4...N : Q − c * Q_sPrev /a_sPrev
	while(++i <= m->N) {
		rPrev = r;
		r = getTridiagRow(m, i);
		if(TINY(rPrev->a_s)) { return SOLVER_FAIL; }
		r->Q_s = r->Q - r->c * rPrev->Q_s / rPrev->a_s;
	}

	return SOLVER_SUCCESS;
}

// compute the x values for the tridiagonal linear system
int computeTrigiag_x(tridiag_t *m) {
	
	assert(m!=NULL);
	
	// case for row N : x = Q_s/a_s
	int i=m->N;
	tridiag_row_t *r = getTridiagRow(m, i);	
	tridiag_row_t *rPrev = r;
	r->x = r->Q_s/r->a_s;

	// case for row N-1, N-2, N-3...1 : x = ( Q_s − b * xPrev ) / a_s
	while(--i >= COMPUTE_x_END) {
		rPrev = r;
		r = getTridiagRow(m, i);
		if(TINY(rPrev->a_s)) { return SOLVER_FAIL; }
		r->x = (r->Q_s - r->b * rPrev->x) / r->a_s;
	}

	return SOLVER_SUCCESS;
}