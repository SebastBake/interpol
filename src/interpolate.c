/***************************************************************************
 * Written by Sebastian Baker for Numerical Algorithms at Unimelb
 ***************************************************************************/

#include "interpolate.h"

double evaluateLagrangeTerm(lagrange_term_t* term, double x);
lagrange_term_t* newLagrangeTerm(interp_set_t* set, int index);
void freeLagrangeTerm(lagrange_term_t* term);

double splineH(int index, cub_spline_t* spline);
int computeCubSplineConstants(interp_set_t* set, cub_spline_t* spline);
int computeCubSplineAs(interp_set_t* set, cub_spline_t* spline);
int computeCubSplineBs(interp_set_t* set, cub_spline_t* spline);
int computeCubSplineCs(interp_set_t* set, cub_spline_t* spline);
int computeCubSplineDs(interp_set_t* set, cub_spline_t* spline);
interp_pt_t* evaluateCubSplineSegment(cub_spline_seg_t* seg, double x);
cub_spline_seg_t* newEmptyCubSplineSegment(int index, double x_lo);
void freeCubSplineSegment(cub_spline_seg_t* seg);

// Create a new data set
interp_set_t* newInterpSet() {

	interp_set_t* newSet = (interp_set_t*)malloc(sizeof(interp_set_t));
	assert( newSet!= NULL );

	newSet->pts = (interp_pt_t**)malloc(INTERP_INIT_ARRLEN*sizeof(interp_pt_t*));
	assert( newSet->pts != NULL);

	newSet->arrLen = INTERP_INIT_ARRLEN;
	newSet->N = 0;

	return newSet;
}

// Create a new point
interp_pt_t* newInterpPt(double x, double fx) {

	interp_pt_t* newPt = (interp_pt_t*)malloc(sizeof(interp_pt_t));
	assert(newPt != NULL);

	newPt->x = x;
	newPt->fx = fx;

	return newPt;
}

// Frees a set (including points)
void freeInterpSet(interp_set_t* set) {
	assert(set!=NULL);
	
	int i=0;
	for(i=0; i < set->N ; i++) {
		freeInterpPt(set->pts[i]);
	}
	free(set->pts);
	free(set);
}

// frees an interpolation point
void freeInterpPt(interp_pt_t* pt) {
	assert(pt!=NULL);
	free(pt);
}

// appends a new point to the set
void appendPtToSet(interp_set_t* set, interp_pt_t* pt) {

	assert(set!=NULL);
	assert(pt!=NULL);

	// extend array if necessary
	if( set->N >= set->arrLen) {
		set->arrLen += INTERP_INIT_ARRLEN;
		int newSize = set->arrLen * sizeof(interp_pt_t*);
		set->pts = (interp_pt_t**)realloc(set->pts, newSize );
	}

	// insert new pt
	set->pts[set->N] = pt;
	set->N++;
}

// generates a lagrange term based on a point in the set
lagrange_term_t* newLagrangeTerm(interp_set_t* set, int index) {

	assert(set!=NULL);

	// initialise new lgrange_term_t
	lagrange_term_t* newTerm = (lagrange_term_t*)malloc(sizeof(lagrange_term_t));
	assert(newTerm != NULL);

	newTerm->order = set->N-1;
	newTerm->index = index;
	newTerm->fx_i = (set->pts[index])->fx;

	newTerm->root = (double*)calloc(newTerm->order+1, sizeof(double));
	assert(newTerm->root!=NULL);

	// calculate denominator
	newTerm->denominator = 1;
	interp_pt_t** pts = set->pts;
	int i=0;
	for(i=0; i <= newTerm->order ; i++) {
		if (i == index) { continue; }
		newTerm->denominator *= ((pts[index])->x - (pts[i])->x);
	}

	// calculate roots
	for(i=0; i <= newTerm->order ; i++) {
		if (i == index) { continue; }
		newTerm->root[i] = (pts[i])->x;
	}

	return newTerm;
}

// generates a lagrange eqn based on a set
lagrange_eqn_t* newLagrangeEqn(interp_set_t* set) {

	assert(set!=NULL);

	// initialise new lagrange_eqn_t
	lagrange_eqn_t* newEqn = (lagrange_eqn_t*)malloc(sizeof(lagrange_eqn_t));
	assert(newEqn != NULL);

	newEqn->num_terms = set->N;

	int arrSize = newEqn->num_terms * sizeof(lagrange_term_t*);
	newEqn->terms = (lagrange_term_t**)malloc(arrSize);
	assert(newEqn->terms != NULL);

	// generate terms
	int i=0;
	for (i=0; i < newEqn->num_terms ; i++) {
		newEqn->terms[i] = newLagrangeTerm(set, i);
	}

	return newEqn;
}

// frees lagrange eqn structure
void freeLagrangeEqn(lagrange_eqn_t* eqn) {
	
	assert(eqn!=NULL);

	// generate terms
	int i=0;
	for (i=0; i < eqn->num_terms ; i++) {
		freeLagrangeTerm(eqn->terms[i]);
	}

	free(eqn->terms);
	free(eqn);
}

// frees lagrange term struct
void freeLagrangeTerm(lagrange_term_t* term) {
	
	assert(term!=NULL);
	free(term->root);
	free(term);
}

// evaluates a lagrange equation at x
interp_pt_t* evaluateLagrangeEqn(lagrange_eqn_t* eqn, double x) {

	assert(eqn!=NULL);

	double fx = 0;

	// sum of evaluated terms
	int i=0;
	for(i=0; i < eqn->num_terms ; i++) {
		fx += evaluateLagrangeTerm(eqn->terms[i],x);
	}

	return newInterpPt(x, fx);
}

// evaluate lagrange term at x
double evaluateLagrangeTerm(lagrange_term_t* term, double x) {

	assert(term!=NULL);
	
	double numerator = term->fx_i;

	int i=0;
	for(i=0; i <= term->order; i++) {
		if(i == term->index) { continue; }
		numerator *= (x - term->root[i]);
	}

	return numerator/term->denominator;
}

// generates a set of cubic splines from a set
cub_spline_t* newCubSpline(interp_set_t* set) {

	assert(set!=NULL);

	cub_spline_t* spline = (cub_spline_t*)malloc(sizeof(cub_spline_t));
	assert(spline!=NULL);

	int segsSize = (set->N)*sizeof(cub_spline_seg_t*);
	spline->segs = (cub_spline_seg_t**)malloc(segsSize);
	assert(spline->segs != NULL);

	spline->num_segs = set->N-1;

	int i=0;
	for(i=0; i < set->N; i++) {
		spline->segs[i] = newEmptyCubSplineSegment(i, (set->pts[i])->x);
	}

	// compute the segment constants
	assert(CUB_SPLINE_COMPUTE_SUCCESS == computeCubSplineConstants(set, spline));

	return spline;
}

// generate a new empty spline segment
cub_spline_seg_t* newEmptyCubSplineSegment(int index, double x_lo) {

	cub_spline_seg_t* seg = (cub_spline_seg_t*)malloc(sizeof(cub_spline_seg_t));
	assert(seg!=NULL);

	seg->index = index;
	seg->x_lo = x_lo;
	seg->a = 0;
	seg->b = 0;
	seg->c = 0;
	seg->d = 0;

	return seg;
}

// free a cubic spline structure including segments
void freeCubSpline(cub_spline_t* spline) {

	assert(spline!=NULL);

	int i=0;
	for (i=0; i<= spline->num_segs; i++) {
		freeCubSplineSegment(spline->segs[i]);
	}

	free(spline->segs);
	free(spline);
}

// free a cubic spline segment structure
void freeCubSplineSegment(cub_spline_seg_t* seg) {
	assert(seg!=NULL);
	free(seg);
}

// get the spline h_i value
double splineH(int index, cub_spline_t* spline) {

	assert(spline!=NULL);
	assert(index < spline->num_segs);

	double x_lo = (spline->segs[index])->x_lo;
	double x_hi = (spline->segs[index+1])->x_lo;

	return x_hi - x_lo;
}

// compute the spline segment constants
int computeCubSplineConstants(interp_set_t* set, cub_spline_t* spline) {

	assert(set!=NULL);
	assert(spline!=NULL);

	int successFlag = CUB_SPLINE_COMPUTE_SUCCESS;

	if( computeCubSplineAs(set, spline) == successFlag && 
		computeCubSplineCs(set, spline) == successFlag && 
		computeCubSplineBs(set, spline) == successFlag && 
		computeCubSplineDs(set, spline) == successFlag
	) {
		return successFlag;
	}
	
	return CUB_SPLINE_COMPUTE_FAIL;
}

// compute spline a's
int computeCubSplineAs(interp_set_t* set, cub_spline_t* spline) {

	assert(set!=NULL);
	assert(spline!=NULL);

	int i=0;
	for(i=0; i <=spline->num_segs; i++) {
		(spline->segs[i])->a = (set->pts[i])->fx;
	}

	return CUB_SPLINE_COMPUTE_SUCCESS;
}

// compute spline b's
int computeCubSplineBs(interp_set_t* set, cub_spline_t* spline) {

	assert(set!=NULL);
	assert(spline!=NULL);

	int i=0;
	double h_i, a_i, a_ip, c_i, c_ip;
	for(i=0; i <spline->num_segs; i++) {

		h_i = splineH(i,spline);
		a_i = (spline->segs[i])->a;
		a_ip = (spline->segs[i+1])->a;
		c_i = (spline->segs[i])->c;
		c_ip = (spline->segs[i+1])->c;

		if(TINY(h_i)) { return CUB_SPLINE_COMPUTE_FAIL; }

		(spline->segs[i])->b = CUB_SPLINE_B(h_i, a_i, a_ip, c_i, c_ip);
	}
	
	return CUB_SPLINE_COMPUTE_SUCCESS;
}

// compute spline c's by solving a tridiagonal matrix
int computeCubSplineCs(interp_set_t* set, cub_spline_t* spline) {

	assert(set!=NULL);
	assert(spline!=NULL);

	// generate tridiag matrix
	tridiag_t *m = newTridiag();
	appendTridiagRow(m, 1.0, 0, 0, 0);
	
	int i=0;
	double h_im, h_i,  a_im, a_i, a_ip;
	double a, b, c, Q;
	for(i=1; i<spline->num_segs; i++) {
		h_im = splineH(i-1,spline);
		h_i = splineH(i,spline);
		a_im = (spline->segs[i-1])->a;
		a_i = (spline->segs[i])->a;
		a_ip = (spline->segs[i+1])->a;

		a = CUB_SPLINE_C_a(h_i, h_im);
		b = h_i;
		c = h_im;
		Q = CUB_SPLINE_C_RHS(h_im, h_i, a_im, a_i, a_ip);

		appendTridiagRow(m, a, b, c, Q);
	}
	appendTridiagRow(m, 1.0, 0, 0, 0);

	// solve tridiag
	if (solveTridiag(m) != SOLVER_SUCCESS) {
		freeTridiag(m);
		return CUB_SPLINE_COMPUTE_FAIL;
	}

	// load c's into spline
	for(i=0; i<=spline->num_segs; i++) {
		(spline->segs[i])->c = getTridiagRow(m, i+1)->x;
	}

	freeTridiag(m);

	return CUB_SPLINE_COMPUTE_SUCCESS;
}

// compute spline d's
int computeCubSplineDs(interp_set_t* set, cub_spline_t* spline) {

	assert(set!=NULL);
	assert(spline!=NULL);

	int i=0;
	double h_i, c_i, c_ip;
	for(i=0; i<spline->num_segs; i++) {

		h_i = splineH(i,spline);
		c_i = (spline->segs[i])->c;
		c_ip = (spline->segs[i+1])->c;

		if(TINY(h_i)) { return CUB_SPLINE_COMPUTE_FAIL; }

		(spline->segs[i])->d = CUB_SPLINE_D(h_i, c_i, c_ip);
	}

	return CUB_SPLINE_COMPUTE_SUCCESS;
}

// evaluate a spline at x
interp_pt_t* evaluateCubSpline(cub_spline_t* spline, double x) {

	assert(spline!=NULL);

	double x_lo = (spline->segs[0])->x_lo;
	double x_hi = (spline->segs[spline->num_segs])->x_lo;
	assert(x>=x_lo && x<=x_hi);

	// find relevant segment
	int i=0;
	while(x>=x_lo) {
		i++;
		x_lo = (spline->segs[i])->x_lo;
	}

	// calculate value
	return evaluateCubSplineSegment(spline->segs[i-1], x);
}

// evaluate a cubic spline segment at x
interp_pt_t* evaluateCubSplineSegment(cub_spline_seg_t* seg, double x) {
	
	assert(seg!=NULL);
	double fx = EVAL_CUB_SPLINE(seg->a, seg->b, seg->c, seg->d, seg->x_lo, x);
	return newInterpPt(x, fx);
}

