#include "_FeatPyramid.h"
#include "math.h"

#define bzero(buf, bytes) ((void) memset (buf, 0, bytes))

double uu[9] = { 1.0000, 0.9397, 0.7660, 0.500, 0.1736, -0.1736, -0.5000, -0.7660, -0.9397 };
double vv[9] = { 0.0000, 0.3420, 0.6428, 0.8660, 0.9848, 0.9848, 0.8660, 0.6428, 0.3420 };

_FeatPyramid::_FeatPyramid ()
{
	_feat = NULL;
	_scales = NULL;
	_dim = -1;
	_imSize = NULL;
}

_FeatPyramid::_FeatPyramid (Matrix<double> *im, int *maxSize, int padX, int padY, double sbin, double interval)
{
	_maxSize = maxSize;
	_featPyramid (im, padX, padY, sbin, interval);
}

_FeatPyramid::~_FeatPyramid()
{
	destroyFeatPyramid();
}

void _FeatPyramid::alphacopy (double *src, double *dst, struct _alphainfo *ofs, int n)
{
	struct _alphainfo *end = ofs + n;
	while (ofs != end)
	{
		dst[ofs->di] += ofs->alpha * src[ofs->si];
		ofs++;
	}
}

void _FeatPyramid::resize1dtran (double *src, int sheight, double *dst, int dheight, int width, int chan)
{
	double scale = (double)dheight/(double)sheight;
	double invscale = (double)sheight/(double)dheight;

	int len = (int)ceil(dheight*invscale) + 2*dheight;

	_alphainfo *ofs = new _alphainfo[len];//!
	int k = 0;

	double fsy1;
	double fsy2;
	int sy1;
	int sy2;

	for (int dy = 0; dy < dheight; dy++)
	{
		fsy1 = dy * invscale;
		fsy2 = fsy1 + invscale;
		sy1 = (int)ceil(fsy1);
		sy2 = (int)floor(fsy2);

		if (sy1 - fsy1 > 1e-3)
		{
			ofs[k].di = dy*width;
			ofs[k].si = sy1-1;
			ofs[k++].alpha = (sy1 - fsy1) * scale;
		}

		for (int sy = sy1; sy < sy2; sy++)
		{
			ofs[k].di = dy*width;
			ofs[k].si = sy;
			ofs[k++].alpha = scale;
		}

		if (fsy2 - sy2 > 1e-3)
		{
			ofs[k].di = dy*width;
			ofs[k].si = sy2;
			ofs[k++].alpha = (fsy2 - sy2) * scale;
		}
	}

	// Resize each column of each color channel
	bzero(dst, chan*width*dheight*sizeof(double));

	double *s;
	double *d;

	for (int c = 0; c < chan; c++)
	{
		for (int x = 0; x < width; x++)
		{
			s = src + c*width*sheight + x*sheight;
			d = dst + c*width*dheight + x;
			alphacopy(s, d, ofs, k);
		}
	}
	delete []ofs;//!
}

Matrix<double>* _FeatPyramid::resize (const Matrix<double> *mxsrc, const double mxscale)
{
	//double *src = mxsrc->data;
	double *src = mxsrc->GetData();
	//int *sdims = (int*)mxsrc->dims;
	int *sdims = mxsrc->GetDims();
	double scale = mxscale;

	int ddims[3];
	ddims[0] = (int)round(sdims[0]*scale);
	ddims[1] = (int)round(sdims[1]*scale);
	ddims[2] = sdims[2];

	Matrix<double>* mxdst = new Matrix<double>(3, ddims);
	double *dst = mxdst->GetData();
	for (int l = 0; l < ddims[0] * ddims[1] * ddims[2]; l++)
		dst[l] = 0;


	double *tmp = new double [ddims[0] * sdims[1] * sdims[2]];
	for (int l = 0; l < ddims[0] * sdims[1] * sdims[2]; l++)
		tmp[l] = 0;

	resize1dtran(src, sdims[0], tmp, ddims[0], sdims[1], sdims[2]);
	resize1dtran(tmp, sdims[1], dst, ddims[1], ddims[0], sdims[2]);

	delete[] tmp;
	return mxdst;
}

int _FeatPyramid::getPaddingX ()
{
	return _maxSize[1];//!
}
int _FeatPyramid::getPaddingY ()
{
	return _maxSize[0];//!
}

Matrix<double>* _FeatPyramid::process (const Matrix<double> *mximage, const double mxsbin)
{
	double *im = mximage->GetData();
	int *dims = mximage->GetDims();
	int sbin = (int) mxsbin;

	// memory for caching orientation histograms & their norms
	int blocks[2];
	blocks[0] = (int)round((double)dims[0]/(double)sbin);
	blocks[1] = (int)round((double)dims[1]/(double)sbin);

	double *hist = new double [blocks[0]*blocks[1]*18];

	for (int i = 0; i < blocks[0]*blocks[1]*18; i++)
		hist[i] = 0;

	double *norm = new double [blocks[0]*blocks[1]];

	for (int i = 0; i < blocks[0]*blocks[1]; i++)
		norm[i] = 0;
	int out[3];

	out[0] = max(blocks[0]-2, 0);
	out[1] = max(blocks[1]-2, 0);
	out[2] = 27+4+1;

	Matrix<double> *mxfeat;
	{
		mxfeat = new Matrix<double>(3, out);
		double *feat = mxfeat->GetData();

		int visible[2];
		visible[0] = blocks[0]*sbin;
		visible[1] = blocks[1]*sbin;

		double *s;
		double dy, dx;
		double v;

		double dy2, dx2;
		double v2;

		double dy3, dx3;
		double v3;

		double best_dot = 0;
		int best_o = 0;

		double dot;

		double xp, yp;

		int ixp, iyp;
		double vx0, vy0, vx1, vy1;

		double *src1, *src2;
		double *dst, *end;

		double *dst2;
		double *src, *p, n1, n2, n3, n4;

		double t1 = 0;
		double t2 = 0;
		double t3 = 0;
		double t4 = 0;

		double h1, h2, h3, h4;
		double sum = 0;
		for (int x = 1; x < visible[1]-1; x++)
		{
			for (int y = 1; y < visible[0]-1; y++)
			{
				// first color channel
				s = im + min(x, dims[1]-2)*dims[0] + min(y, dims[0]-2);
				dy = *(s+1) - *(s-1);
				dx = *(s+dims[0]) - *(s-dims[0]);
				v = dx*dx + dy*dy;

				// second color channel
				if(dims[2] > 1)
				{
					s += dims[0]*dims[1];
					dy2 = *(s+1) - *(s-1);
					dx2 = *(s+dims[0]) - *(s-dims[0]);
					v2 = dx2*dx2 + dy2*dy2;
					if (v2 > v)
					{
						v = v2;
						dx = dx2;
						dy = dy2;
					}
				}

				// third color channel
				if(dims[2] > 2)
				{
					s += dims[0]*dims[1];
					dy3 = *(s+1) - *(s-1);
					dx3 = *(s+dims[0]) - *(s-dims[0]);
					v3 = dx3*dx3 + dy3*dy3;
					if (v3 > v)
					{
						v = v3;
						dx = dx3;
						dy = dy3;
					}
				}

				best_dot = 0;
				best_o = 0;

				for (int o = 0; o < 9; o++)
				{
					dot = uu[o]*dx + vv[o]*dy;
					if (dot > best_dot)
					{
						best_dot = dot;
						best_o = o;
					}

					else if (-dot > best_dot)
					{
						best_dot = -dot;
						best_o = o+9;
					}
				}

				xp = ((double)x+0.5)/(double)sbin - 0.5;
				yp = ((double)y+0.5)/(double)sbin - 0.5;

				ixp = (int)floor(xp);
				iyp = (int)floor(yp);
				vx0 = xp-ixp;
				vy0 = yp-iyp;
				vx1 = 1.0-vx0;
				vy1 = 1.0-vy0;
				v = sqrt(v);
				if (ixp >= 0 && iyp >= 0)
				{
					*(hist + ixp*blocks[0] + iyp + best_o*blocks[0]*blocks[1]) += vx1*vy1*v;
				}

				if (ixp+1 < blocks[1] && iyp >= 0)
				{
					*(hist + (ixp+1)*blocks[0] + iyp + best_o*blocks[0]*blocks[1]) += vx0*vy1*v;
				}

				if (ixp >= 0 && iyp+1 < blocks[0])
				{
					*(hist + ixp*blocks[0] + (iyp+1) + best_o*blocks[0]*blocks[1]) += vx1*vy0*v;
				}

				if (ixp+1 < blocks[1] && iyp+1 < blocks[0])
				{
					*(hist + (ixp+1)*blocks[0] + (iyp+1) + best_o*blocks[0]*blocks[1]) += vx0*vy0*v;
				}
			}
		}
		for (int o = 0; o < 9; o++)
		{
			src1 = hist + o*blocks[0]*blocks[1];
			src2 = hist + (o+9)*blocks[0]*blocks[1];
			dst = norm;
			end = norm + blocks[1]*blocks[0];

			while (dst < end)
			{
				*(dst++) += (*src1 + *src2) * (*src1 + *src2);
				src1++;
				src2++;
			}
		}
		for (int x = 0; x < out[1]; x++)
		{
			for (int y = 0; y < out[0]; y++)
			{
				dst2 = feat + x*out[0] + y;

				p = norm + (x+1)*blocks[0] + y+1;
				n1 = 1.0 / sqrt(*p + *(p+1) + *(p+blocks[0]) + *(p+blocks[0]+1) + eps);
				p = norm + (x+1)*blocks[0] + y;
				n2 = 1.0 / sqrt(*p + *(p+1) + *(p+blocks[0]) + *(p+blocks[0]+1) + eps);
				p = norm + x*blocks[0] + y+1;
				n3 = 1.0 / sqrt(*p + *(p+1) + *(p+blocks[0]) + *(p+blocks[0]+1) + eps);
				p = norm + x*blocks[0] + y;
				n4 = 1.0 / sqrt(*p + *(p+1) + *(p+blocks[0]) + *(p+blocks[0]+1) + eps);

				t1 = 0;
				t2 = 0;
				t3 = 0;
				t4 = 0;

				// contrast-sensitive features
				src = hist + (x+1)*blocks[0] + (y+1);

				for (int o = 0; o < 18; o++)
				{
					h1 = min(*src * n1, 0.2);
					h2 = min(*src * n2, 0.2);
					h3 = min(*src * n3, 0.2);
					h4 = min(*src * n4, 0.2);
					*dst2 = 0.5 * (h1 + h2 + h3 + h4);
					t1 += h1;
					t2 += h2;
					t3 += h3;
					t4 += h4;
					dst2 += out[0]*out[1];
					src += blocks[0]*blocks[1];
				}

				// contrast-insensitive features
				src = hist + (x+1)*blocks[0] + (y+1);

				for (int o = 0; o < 9; o++)
				{
					sum = *src + *(src + 9*blocks[0]*blocks[1]);
					h1 = min(sum * n1, 0.2);
					h2 = min(sum * n2, 0.2);
					h3 = min(sum * n3, 0.2);
					h4 = min(sum * n4, 0.2);
					*dst2 = 0.5 * (h1 + h2 + h3 + h4);
					dst2 += out[0]*out[1];
					src += blocks[0]*blocks[1];
				}

				*dst2 = 0.2357 * t1;
				dst2 += out[0]*out[1];
				*dst2 = 0.2357 * t2;
				dst2 += out[0]*out[1];
				*dst2 = 0.2357 * t3;
				dst2 += out[0]*out[1];
				*dst2 = 0.2357 * t4;

				dst2 += out[0]*out[1];
				*dst2 = 0;
			}
		}

		delete[] hist;
		delete[] norm;
	}
	return mxfeat;
}


Matrix<double>* _FeatPyramid::padArray (Matrix<double> *mat, int* dimPad, double __alphainfo)
{
	int dims[3];

	if (mat == NULL)
	{
		for (int i = 0; i < 2; i++)
			dims[i] = dimPad[i] * 2;

		dims[2] = 32 + (dimPad[2] * 2);
	}
	else
	{
		for (int i = 0; i < 3; i++)
			dims[i] = mat->GetDims()[i] + (dimPad[i] * 2);
	}

	Matrix<double>* padMat = new Matrix<double>(3, dims);
	double* paddata = padMat->GetData();

	int psize = dims[0]*dims[1]*dims[2];
	for(int i = 0; i < psize; i++)
		paddata[i] = __alphainfo;

	if (mat != NULL)
	{
		double *matdata = mat->GetData();
		int md1 = padMat->GetDims()[0]*padMat->GetDims()[1], md2 = padMat->GetDims()[0];
		int md3 = mat->GetDims()[0]*mat->GetDims()[1], md4 = mat->GetDims()[0];
		for (int i = dimPad[0]; i < padMat->GetDims()[0] - dimPad[0]; i++)
			for (int j = dimPad[1]; j < padMat->GetDims()[1] - dimPad[1]; j++)
				for (int k = dimPad[2]; k < padMat->GetDims()[2] - dimPad[2]; k++)
					paddata[i+j*md2+k*md1] = matdata[(i-dimPad[0])+(j-dimPad[1])*md4+(k-dimPad[2])*md3];
	}

	return padMat;
}

void _FeatPyramid::_featPyramid (const Matrix<double> *im, int padX, int padY, double sbin, double interval)
{
	double sc;
	int imsize[2];
	double maxScale;
	Matrix<double>* imAux;
	int pad[3];

	if (padX == -1 && padY == -1)
	{
		padX = getPaddingX();
		padY = getPaddingY();
	}

	sc = pow (2, 1/(interval));
	imsize[0] = im->GetDims()[0];
	imsize[1] = im->GetDims()[1];
	maxScale = 1 + floor (log (min (imsize[0], imsize[1]) /
		(5*sbin)) / log(sc));

	if ((maxScale + interval) < (2*interval))
		setDim (int(2*interval));

	else
		setDim (int(maxScale + interval));

	_feat = new Matrix<double>*[getDim()];
	_scales = new double [getDim()];

	setImSize (imsize);
	for (int i = 0; i < interval; i++)
	{
		imAux = resize (im, (1/pow(sc, i)));

		setFeat(process(imAux, sbin/2), i);
		setScales(2/pow(sc, i), i);
		setFeat (process (imAux, sbin), (int)(i+interval));
		setScales (1/pow(sc, i), (int)(i+interval));

		// Remaining intervals
		Matrix<double>* imAux2; // mjmarin added
		for (int j = (int)(i+interval); j < (int)maxScale; j = j+(int)interval)
		{
			imAux2 = resize (imAux, 0.5);
			//imAux = imAux2;

			setFeat (process (imAux2, sbin), (int)(j+interval));
			setScales (0.5 * getScales()[j], (int)(j+interval));
			delete imAux2;
		}
		delete imAux;
	}
	for (int i = 0; i < getDim(); i++)
	{
		// Add 1 to padding because feature generation deletes a 1-cell
		// Wide border around the feature map
		pad[0] = padY + 1;
		pad[1] = padX + 1;
		pad[2] = 0;

		Matrix<double>* tmpfeat = getFeat()[i];
		Matrix<double>* tmppad = padArray (tmpfeat, pad, 0);
		delete _feat[i];
		setFeat (tmppad, i);

		// Write boundary occlusion feature
		Matrix<double> *tfeat = getFeat(i);
		int dm1 = tfeat->GetDims()[0]*tfeat->GetDims()[1]*31, dm2 = tfeat->GetDims()[0];

		double *tfd = tfeat->GetData();
		for(int j = 0; j < padY; j++)
			for(int k = 0; k < tfeat->GetDims()[1]; k++)
				tfd[j+k*dm2+dm1] = 1;

		for (int j = tfeat->GetDims()[0] - padY -1; j < tfeat->GetDims()[0]; j++)
			for (int k = 0; k < tfeat->GetDims()[1]; k++)
				tfd[j+k*dm2+dm1] = 1;

		for (int j = 0; j < tfeat->GetDims()[0]; j++)
			for (int k = 0; k <= padX; k++)
				tfd[j+k*dm2+dm1] = 1;

		for (int j = 0; j < tfeat->GetDims()[0]; j++)
			for (int k = tfeat->GetDims()[1] - padX - 1; k < tfeat->GetDims()[1]; k++)
				tfd[j+k*dm2+dm1] = 1;
	}

	setPadX (padX);
	setPadY (padY);
}

void _FeatPyramid::destroyFeatPyramid()
{
	if (_feat != NULL)
		for (int i = 0; i < getDim(); i++) delete _feat[i];
	if(_feat)
		delete []_feat;

	if (_scales != NULL)
	{
		delete[] _scales;
		_scales = NULL;
	}
	if(_imSize)
		delete []_imSize;
}

void _FeatPyramid::setImSize(int *imSize)
{
	_imSize = new int [2];

	_imSize[0] = imSize[0];
	_imSize[1] = imSize[1];
}