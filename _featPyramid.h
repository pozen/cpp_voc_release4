#ifndef _FeatPyramid_H_
#define _FeatPyramid_H_

#include <iostream>
#include "Matrix.h"
using namespace std;

#define eps 0.0001
#define round(x) (((x)-floor(x))>=(0.5)?(long)(ceil(x)):(long)(floor(x)))

struct _alphainfo 
{
	int si;
	int di;
	double alpha;
};

class _FeatPyramid
{
public:

	_FeatPyramid ();

	_FeatPyramid(Matrix<double> *im, int *maxSize, int padX, int padY, double sbin, double interval);

	~_FeatPyramid ();

	void alphacopy(double *src, double *dst, struct _alphainfo *ofs, int n);

	void resize1dtran(double *src, int sheight, double *dst, int dheight, 
		int width, int chan);

	Matrix<double>* resize(const Matrix<double> *mxsrc, const double mxscale);

	int getPaddingX ();

	int getPaddingY ();

	Matrix<double>* process(const Matrix<double> *mximage, const double mxsbin);

	Matrix<double>* padArray (Matrix<double> *mat, int dimPad[3], double val);

	void _featPyramid (const Matrix<double> *im, int padX, int padY, double sbin, double interval);

	void destroyFeatPyramid();

	Matrix<double>** getFeat() const	{return _feat;}

	void setFeat (Matrix<double> **feat)
	{
		if (feat != NULL)
			_feat = feat;
	}

	void setFeat (Matrix<double>* feat, int i)
	{
		_feat[i] = feat;
	}

	Matrix<double>* getFeat(int i)const
	{
		return _feat[i];
	}

	double* getScales() const	{return _scales;}

	void setScales (double *scales)
	{
		if (scales != NULL)
			_scales = scales;
	}

	void setScales (double scales, int i)
	{
		_scales[i] = scales;
	}

	int getDim() const	{return _dim;}

	void setDim (int dim)
	{
		_dim = dim;
	}

	int* getImSize() const	{ return _imSize; }
	void setImSize (int *imSize);
	int getPadX() const	{ return _padX; }
	void setPadX (int padX)	{ _padX = padX; }
	int getPadY() const	{ return _padY; }
	void setPadY (int padY)	{ _padY = padY; }

private:
	Matrix<double> **_feat;
	double *_scales;
	int _dim;
	int *_imSize;
	int *_maxSize;
	int _padX;
	int _padY;
};

#endif