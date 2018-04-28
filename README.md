# interpol
A C library for interpolating data.

Using the lagrange method will produce a polynomial of order n-1 through each of the n points provided.

Using the cubic spline method will produce a set of cubic polynomials where the gradient is continuous at each point.

## Usage

### Create a set of 2d points to interpolate

```C
interp_set_t* set = newInterpSet();

//Note: Ordering of points added to the set is not important
interp_pt_t* p1 = newInterpPt(x1, y1);
appendPtToSet(set, p1);

interp_pt_t* p2 = newInterpPt(x2, y2);
appendPtToSet(set, p2);

//...

interp_pt_t* pn = newInterpPt(xn, yn);
appendPtToSet(set, pn);
```

### Interpolate a set at ```a```

#### Lagrange

```C
// Generate a lagrange polynomial from an interp_set_t*
lagrange_eqn_t* lag_poly = newLagrangeEqn(set);

//  Evaluate at ```a```
interp_pt_t* a_y = evaluateLagrangeEqn(lag_poly, a_x);
```

#### Spline

```C
// Generate a cubic spline from an interp_set_t*
cub_spline_t* spline = newCubSpline(set);

//  Evaluate at ```a```
interp_pt_t* a_y = evaluateCubSpline(spline, a_x);
````




