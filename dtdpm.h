#ifndef _DTDPM_H_
#include "_FeatPyramid.h"

typedef struct ObjBox {
	int x1;
	int x2;
	int y1;
	int y2;
	float score;
}ObjBox;

class DTDPM {
public:
	DTDPM(){};
	DTDPM( char *filename );
	~DTDPM();
public:
	void Save( char *filename );
	void Load( char *filename );
	void CascadeDetect( float thre, float overlap, _FeatPyramid *pyra, ObjBox **pick, int *pickNum );

	int* GetMaxSize(){ return maxsize; }
	int GetMaxSbin(){ return sbin; }
	int GetInterval(){ return interval; }
	
protected:
	_FeatPyramid* ProjPyra( _FeatPyramid *pyra );
	void fconv( Matrix<double> *features, double **rootfilter, int a, int filtersDim, Matrix<double>** result );
	void process( Matrix<double> *mxA, Matrix<double> *mxB, Matrix<double> *mxC );
	Matrix<double>* cascade( _FeatPyramid *fp1, _FeatPyramid *fp2, Matrix<double> ***_rootscores, int _numrootlocs, double *fscales, int _padx, int _pady, int s );
	double conv( int x, int y, const double *A, const unsigned int *A_dims, const double *B, const unsigned int *B_dims, int num_features );
	double rconv( int L, int filterind, int x, int y, int pca );
	double pconvdt( int L, int probex, int probey, int filterind, int defindex, int xstart, int xend,\
		int ystart, int yend, int pca, double defthresh ); 
	double partscore( int L, int defindex, int pfind, int x, int y, int pca, double defthresh );
	void  Init( _FeatPyramid *fp1, _FeatPyramid *fp2, int s );
	void Cleanup();

private:
	int    sbin;
	int    pcadim;
	int    numfeatures;
	double **rootfilters;
	unsigned int\
		   **rootfilterdims;
	double **partfilters[2];
	unsigned int\
		   **partfilterdims;
	double *offsets;
	int    numcomponents;
	int    *numparts;
	double thresh;
	int    **partorder;
	double **t;
	int    *tdim;
	double ***anchors;
	int    **anchorsdim;
	double **defs;
	int    defdim;

	int    **pfind;
	int    **defind;
	int    numpartfilters;
	int    numdefparams;
	int    numlevels;
	int    **featdims;
	int    *featdimsprod;
	int    interval;
	int    numrootlocs;
	Matrix<double> coeff;
	//_FeatPyramid\
		   *pyra;
	int    maxsize[2];

	double **rootfilterspca;
	unsigned int\
		   **rootfilterspcadims;
	int    *rootIndex;
	int    *offsetIndex;

	// pyramid level offsets for convolution
	int *LOFFCONV;
	// pyramid level offsets for distance transform
	int *LOFFDT;
	// convolution values
	double *PCONV[2];
	// distance transform values
	double *DT[2];
	// distance transform argmaxes, x dimension
	int *DXAM[2];
	// distance transform argmaxes, y dimension
	int *DYAM[2];
	// half-width of distance transform window
	const static int S = 5;
	// padding used in the feature pyramid
	int padx, pady;
	// precomputed deformation costs
	double **DXDEFCACHE;
	double **DYDEFCACHE;

	// feature pyramid levels
	// feat[0] holds non-PCA HOG features
	// feat[1] holds PCA of HOG features
	double **feat[2];

	//! tmp buffer
	void   *tmpMem;
};
#endif