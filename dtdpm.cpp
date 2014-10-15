#include "DTDPM.h"
#include <limits>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

static double  INFINITY = numeric_limits<double>::infinity();
static inline int square(int x) { return x*x; }

DTDPM::DTDPM(char *filename) {
	Load(filename);
}

DTDPM::~DTDPM() {
	for (int i = 0; i < numpartfilters; i++) {
		delete []partfilters[0][i];
		delete []partfilters[1][i];
	}
	delete []partfilterdims;
	delete []partfilters[0];
	delete []partfilters[1];

	for(int i = 0; i < numdefparams; i++)
		delete []defs[i];
    delete []defs;
	delete []offsets;
	delete []tdim;

	for (int i = 0; i < numcomponents; i++) {
		delete []rootfilterdims[i];
		delete []rootfilters[i];
		delete []rootfilterspcadims[i];
		delete []rootfilterspca[i];
		delete []pfind[i];
		delete []defind[i];
		delete []partorder[i];
		delete []anchorsdim[i];
		delete []t[i];
		for(int j = 0; j < numparts[i]; j++)
			delete []anchors[i][j];
		delete []anchors[i];
	}
	delete []rootfilterdims;
	delete []rootfilters;
	delete []rootfilterspcadims;
	delete []rootfilterspca;
	delete []pfind;
	delete []defind;
	delete []partorder;
	delete []anchorsdim;
	delete []anchors;
	delete []numparts;
	delete [] rootIndex;
	delete [] offsetIndex;
	delete []t;
}

void DTDPM::Save(char *filename) {
	fstream file;
	//! have to use binary mode
	file.open(filename, ios_base::out | ios::binary);

	file.write((const char*)&thresh, sizeof(int));
	file.write((const char*)&interval, sizeof(int));
	file.write((const char*)&numcomponents, sizeof(int));
	file.write((const char*)&sbin, sizeof(int));
	file.write((const char*)&numpartfilters, sizeof(int));
	file.write((const char*)&numdefparams, sizeof(int));
	file.write((const char*)numparts, sizeof(int)*numcomponents);
	for (int i = 0; i < numpartfilters; i++) {
		file.write((const char*)partfilterdims[i], sizeof(int)*3);
		int wdsize = partfilterdims[i][0]*partfilterdims[i][1]*partfilterdims[i][2];
		file.write((const char*)partfilters[0][i], sizeof(double)*wdsize);
		file.write((const char*)partfilters[1][i], sizeof(double)*wdsize);
	}

	file.write((const char*)&defdim, sizeof(int));
	for(int i = 0; i < numdefparams; i++) {
		file.write((const char*)defs[i], sizeof(double)*defdim);
	}
	file.write((const char*)offsets, sizeof(double)*numcomponents);
	file.write((const char*)tdim, sizeof(int)*numcomponents);
	for (int i = 0; i < numcomponents; i++) {
		file.write((const char*)rootfilterdims[i], sizeof(int)*3);
		int rfwsize = rootfilterdims[i][0]*rootfilterdims[i][1]*rootfilterdims[i][2];
		file.write((const char*)rootfilters[i], sizeof(double)*rfwsize);

		file.write((const char*)rootfilterspcadims[i], sizeof(int)*3);
		rfwsize = rootfilterspcadims[i][0]*rootfilterspcadims[i][1]*rootfilterspcadims[i][2];
		file.write((const char*)rootfilterspca[i], sizeof(double)*rfwsize);

		file.write((const char*)pfind[i], sizeof(int)*numparts[i]);
		file.write((const char*)defind[i], sizeof(int)*numparts[i]);
		file.write((const char*)partorder[i], sizeof(int)*(2*numparts[i]+2));
		file.write((const char*)anchorsdim[i], sizeof(int)*numparts[i]);
		for(int j = 0; j < numparts[i]; j++) {
			file.write((const char*)anchors[i][j], sizeof(double)*anchorsdim[i][j]);
		}
	}
	file.write((const char*)rootIndex, sizeof(int)*numcomponents);
	file.write((const char*)offsetIndex, sizeof(int)*numcomponents);

	int *coeffdims = coeff.GetDims();
	double *coeffdata = coeff.GetData();
	file.write((const char*)coeffdims, sizeof(int)*2);
	file.write((const char*)coeffdata, sizeof(double)*coeffdims[0]*coeffdims[1]);

	file.write((const char*)maxsize, sizeof(int)*2);

	file.write((const char*)tdim, sizeof(int)*numcomponents);
	for(int i = 0; i < numcomponents; i++)
		file.write((const char*)t[i], sizeof(double)*tdim[i]);

	file.close();
}

void DTDPM::Load(char *filename) {
	//! alloc buffers for model && load model from 'filename'
	fstream file;
	file.open(filename, ios_base::in | ios::binary);

	file.read((char*)&thresh, sizeof(int));
	file.read((char*)&interval, sizeof(int));
	file.read((char*)&numcomponents, sizeof(int));
	file.read((char*)&sbin, sizeof(int));
	file.read((char*)&numpartfilters, sizeof(int));
	file.read((char*)&numdefparams, sizeof(int));

	numparts        = new int[numcomponents]; 
	anchors         = new double**[numcomponents];
	defs            = new double*[numdefparams];
	rootfilters     = new double*[numcomponents];
	rootfilterspca  = new double*[numcomponents];
	partfilters[0]  = new double*[numpartfilters];
	partfilters[1]  = new double*[numpartfilters];
	rootfilterdims  = new unsigned int*[numcomponents];
	rootfilterspcadims \
		            = new unsigned int*[numcomponents];
	partfilterdims  = new unsigned int*[numpartfilters];
	pfind           = new int*[numcomponents];
	defind          = new int*[numcomponents];
	anchorsdim      = new int*[numcomponents];
	partorder       = new int*[numcomponents];

	file.read((char*)numparts, sizeof(int)*numcomponents);

	for (int i = 0; i < numpartfilters; i++) {
		partfilterdims[i] = new unsigned int[3];
		file.read((char*)partfilterdims[i], sizeof(int)*3);
		int wdsize = partfilterdims[i][0]*partfilterdims[i][1]*partfilterdims[i][2];
		partfilters[0][i] = new double[wdsize];
		partfilters[1][i] = new double[wdsize];
		file.read((char*)partfilters[0][i], sizeof(double)*wdsize);
		file.read((char*)partfilters[1][i], sizeof(double)*wdsize);
	}

	file.read((char*)&defdim, sizeof(int));
	for(int i = 0; i < numdefparams; i++) {
		defs[i] = new double[defdim];
		file.read((char*)defs[i], sizeof(double)*defdim);
	}

	offsets = new double[numcomponents];
	file.read((char*)offsets, sizeof(double)*numcomponents);
	tdim = new int[numcomponents];
	file.read((char*)tdim, sizeof(int)*numcomponents);

	for (int i = 0; i < numcomponents; i++) {
		rootfilterdims[i] = new unsigned int[3];
		file.read((char*)rootfilterdims[i], sizeof(int)*3);
		int rfwsize = rootfilterdims[i][0]*rootfilterdims[i][1]*rootfilterdims[i][2];
		rootfilters[i] = new double[rfwsize];
		file.read((char*)rootfilters[i], sizeof(double)*rfwsize);

		rootfilterspcadims[i] = new unsigned int[3];
		file.read((char*)rootfilterspcadims[i], sizeof(int)*3);
		rfwsize = rootfilterspcadims[i][0]*rootfilterspcadims[i][1]*rootfilterspcadims[i][2];
		rootfilterspca[i] = new double[rfwsize];
		file.read((char*)rootfilterspca[i], sizeof(double)*rfwsize);

		pfind[i] = new int[numparts[i]];
		file.read((char*)pfind[i], sizeof(int)*numparts[i]);
		defind[i] = new int [numparts[i]];
		file.read((char*)defind[i], sizeof(int)*numparts[i]);
		partorder[i] = new int[2*numparts[i]+2];
		file.read((char*)partorder[i], sizeof(int)*(2*numparts[i]+2));
		anchorsdim[i] = new int[numparts[i]];
		file.read((char*)anchorsdim[i], sizeof(int)*numparts[i]);
		anchors[i] = new double*[numparts[i]];
		for(int j = 0; j < numparts[i]; j++) {
			anchors[i][j] = new double[anchorsdim[i][j]];
			file.read((char*)anchors[i][j], sizeof(double)*anchorsdim[i][j]);
		}
	}
	rootIndex = new int[numcomponents];
	file.read((char*)rootIndex, sizeof(int)*numcomponents);
	offsetIndex = new int[numcomponents];
	file.read((char*)offsetIndex, sizeof(int)*numcomponents);

	int dims[2];
	file.read((char*)dims, sizeof(int)*2);
	double *data = new double[dims[0]*dims[1]];
	file.read((char*)data, sizeof(double)*dims[0]*dims[1]);
	coeff.Set(2, dims, data);
	delete []data;

	file.read((char*)maxsize, sizeof(int)*2);

	tdim = new int[numcomponents];
	t = new double*[numcomponents];
	file.read((char*)tdim, sizeof(int)*numcomponents);
	for(int i = 0; i < numcomponents; i++) {
		t[i] = new double[tdim[i]];
		file.read((char*)t[i], sizeof(double)*tdim[i]);
	}

	file.close();
}

_FeatPyramid* DTDPM::ProjPyra(_FeatPyramid *pyra) {
	//! create proj pyra
	_FeatPyramid *projpyra = new _FeatPyramid;
	projpyra->setDim(pyra->getDim());
	int* size = new int[2];
	memcpy(size, pyra->getImSize(), 2*sizeof(int));
	projpyra->setImSize(size);

	projpyra->setPadX(pyra->getPadX());
	projpyra->setPadY(pyra->getPadY());

	Matrix<double> **tmpFeat = new Matrix<double>*[pyra->getDim()];
	projpyra->setFeat(tmpFeat);
	double *scale = new double[pyra->getDim()];
	memcpy(scale, pyra->getScales(), pyra->getDim()*sizeof(double));
	projpyra->setScales(scale);

	for(int i = 0; i < pyra->getDim(); i++) {
		int fdim = 3;
		int size[7];
		int *fdims = pyra->getFeat(i)->GetDims();
		int *cdims = coeff.GetDims();
		for(int j = 0; j < fdim; j++) size[j] = fdims[j];
		size[fdim] = size[0]*size[1];
		size[fdim+1] = size[2];
		size[fdim+2] = size[fdim];
		size[fdim+3] = cdims[1];
		size[2] = size[fdim+3];

		Matrix<double> src(*pyra->getFeat(i));
		src.ReShape(2, &size[fdim]);
		Matrix<double>* dst = new Matrix<double>(2, &size[fdim+2]);
		Matrix<double>::Mul(&src, &coeff, dst);
		dst->ReShape(3, size);
		projpyra->setFeat(dst, i);
	}
	return projpyra;
}

static bool cmp_boxs(const ObjBox & a, const ObjBox & b){ return a.score > b.score; }

void DTDPM::CascadeDetect(float thre, float overlap, _FeatPyramid *pyra, ObjBox **pick, int *pickNum) {
	int numrootfilters = numcomponents;

	_FeatPyramid *projpyra = ProjPyra(pyra);

	int numrootlocs = 0;
	int nlevels = pyra->getDim();
	Matrix<double> ***rootscores = new Matrix<double>**[numcomponents];
	for(int i = 0; i < numcomponents; i++)
		rootscores[i] = new Matrix<double>*[nlevels];
	for(int i = 0; i < numcomponents; i++)
		for(int j = 0; j < pyra->getDim(); j++)
			rootscores[i][j] = 0;
	int s = 0;

	Matrix<double> **tmpScores = new Matrix<double>*[pyra->getDim()];
	for(int i = 0; i < pyra->getDim(); i++) {
		Matrix<double> *tmpFeat = pyra->getFeat(i);
		s += tmpFeat->GetDims()[0]*tmpFeat->GetDims()[1];
		tmpScores[i] = 0;
		if(i >= interval) {
			/*Matrix<double>* scores = */
			fconv(projpyra->getFeat()[i], rootfilterspca, 1, numrootfilters, &tmpScores[i]);
			for(int c = 0; c < numcomponents; c++) {
				int u = rootIndex[c] - 1;
				int v = offsetIndex[c] - 1;
				tmpScores[i][u].Add(offsets[v]);
				rootscores[c][i] = &tmpScores[i][u];
				numrootlocs += tmpScores[i][u].GetDims()[0]*tmpScores[i][u].GetDims()[1];
			}
		}
	}
	s *= numpartfilters;
	thresh = thre;
	Matrix<double> *det = cascade(pyra, projpyra, rootscores, numrootlocs, (double*)pyra->getScales(),\
		pyra->getPadX(), pyra->getPadY(), s);

	if(0 == det) {
		*pick = 0;
		*pickNum = 0;
		for(int i = 0; i < pyra->getDim(); i++)
		{
			if(tmpScores[i])
				delete []tmpScores[i];
		}
		delete tmpScores;
		delete projpyra;
		for(int i = 0; i < numcomponents; i++)
			delete []rootscores[i];
		delete []rootscores;
		return;
	}
	vector<ObjBox> boxs;
	int objNum = det->GetDataSize()/5;
	double *objData = det->GetData();
	int *areas = new int[objNum];
	int *mark = new int[objNum];
	for(int i = 0; i < objNum; i++) {
		ObjBox box;
		box.x1 = objData[i*5+0];
		box.y1 = objData[i*5+1];
		box.x2 = objData[i*5+2];
		box.y2 = objData[i*5+3];
		box.score = objData[i*5+4];
		boxs.push_back(box);
		areas[i] = (boxs[i].x2 - boxs[i].x1 + 1)*(boxs[i].y2 - boxs[i].y1 + 1);
		mark[i] = 0;
	}
	sort(boxs.begin(), boxs.end(), cmp_boxs);
	int pickDim = 0;
	for(int i = 0; i < objNum; i++) {
		if(mark[i] != 0)
			continue;
		mark[i] = 2;
		pickDim++;
		for(int j = objNum - 1; j >= 0; j--) {
			if(mark[j]) continue;
			int xx1 = boxs[i].x1 > boxs[j].x1 ? boxs[i].x1 : boxs[j].x1;
			int yy1 = boxs[i].y1 > boxs[j].y1 ? boxs[i].y1 : boxs[j].y1;
			int xx2 = boxs[i].x2 < boxs[j].x2 ? boxs[i].x2 : boxs[j].x2;
			int yy2 = boxs[i].y2 < boxs[j].y2 ? boxs[i].y2 : boxs[j].y2;
			float w = xx2-xx1+1;
			float h = yy2-yy1+1;
			if (w > 0 && h > 0) {
				float o = w * h / areas[i];
				if (o > overlap)
					mark[j] = 1;
			}
		}
	}
	ObjBox *rboxs = new ObjBox[pickDim];
	int tmp = 0;
	for(int i = 0; i < objNum; i++)
		if(mark[i] == 2)
			rboxs[tmp++] = boxs[i];
	*pick = rboxs;
	*pickNum = tmp;
	delete det;
	delete mark;
	delete areas;

	for(int i = 0; i < pyra->getDim(); i++) {
		if(tmpScores[i])
			delete []tmpScores[i];
	}
	delete tmpScores;
	delete projpyra;
	for(int i = 0; i < numcomponents; i++)
		delete []rootscores[i];
	delete []rootscores;
}

void DTDPM::fconv(Matrix<double> *features, double **rootfilter, int a, int filtersDim, Matrix<double>** result) {
	int num_bs = filtersDim;  

	int start = a - 1;
	int end = filtersDim - 1;

	int len = end - start + 1;
	*result = new Matrix<double>[len];

	Matrix<double> mxB;
	int height;
	int width;

	for (int i = 0; i < len; i++)  {
		mxB.Set(3, (const int*)rootfilterspcadims[i+start], rootfilter[i+start]);
		height = features->GetDims()[0] - mxB.GetDims()[0] + 1;
		width = features->GetDims()[1] - mxB.GetDims()[1] + 1;

		int dims[2];
		dims[0] = height;
		dims[1] = width;
        (*result)[i].Set(2, dims, 0);

		process(features, &mxB, &(*result)[i]);
	}
}

void DTDPM::process(Matrix<double> *mxA, Matrix<double> *mxB, Matrix<double> *mxC) {
	double *A = mxA->GetData();
	double *B = mxB->GetData();
	double *C = mxC->GetData();
	const int *ADims = mxA->GetDims();
	const int *BDims = mxB->GetDims();
	const int *CDims = mxC->GetDims();
	int num_features = ADims[2];

	double *dst;
	double *A_src;
	double *B_src;
	double val = 0;
	double *A_off;
	double *B_off;

	for (int f = 0; f < num_features; f++) {
		dst = C;
		A_src = A + f*ADims[0]*ADims[1];    
		B_src = B + f*BDims[0]*BDims[1];

		for (int x = 0; x < CDims[1]; x++) {
			for (int y = 0; y < CDims[0]; y++) {
				val = 0;

				for (int xp = 0; xp < BDims[1]; xp++) {
					A_off = A_src + (x+xp)*ADims[0] + y;
					B_off = B_src + xp*BDims[0];

					switch (BDims[0]) {
					case 20: val += A_off[19] * B_off[19];
					case 19: val += A_off[18] * B_off[18];
					case 18: val += A_off[17] * B_off[17];
					case 17: val += A_off[16] * B_off[16];
					case 16: val += A_off[15] * B_off[15];
					case 15: val += A_off[14] * B_off[14];
					case 14: val += A_off[13] * B_off[13];
					case 13: val += A_off[12] * B_off[12];
					case 12: val += A_off[11] * B_off[11];
					case 11: val += A_off[10] * B_off[10];
					case 10: val += A_off[9] * B_off[9];
					case 9: val += A_off[8] * B_off[8];
					case 8: val += A_off[7] * B_off[7];
					case 7: val += A_off[6] * B_off[6];
					case 6: val += A_off[5] * B_off[5];
					case 5: val += A_off[4] * B_off[4];
					case 4: val += A_off[3] * B_off[3];
					case 3: val += A_off[2] * B_off[2];
					case 2: val += A_off[1] * B_off[1];
					case 1: val += A_off[0] * B_off[0];
						break;
					default:        
						for (int yp = 0; yp < BDims[0]; yp++) 
							val += *(A_off++) * *(B_off++);
					}
				}
				*(dst++) += val;
			}
		}
	}
}

Matrix<double>* DTDPM::cascade(_FeatPyramid *fp1, _FeatPyramid *fp2, Matrix<double> ***_rootscores, int _numrootlocs, double *fscales, int _padx, int _pady, int s) {
	Init(fp1, fp2, s);

	Matrix<double> ***rootscores = _rootscores;
	int numrootlocs           = _numrootlocs;
	padx                      = _padx;
	pady                      = _pady;

	double **pcascore = new double*[numcomponents];
	for (int c = 0; c < numcomponents; c++)
		pcascore[c] = new double[numparts[c]+1];

	vector<double> coords;

	int nlevels = numlevels - interval;
	unsigned int dim[2];

	for (int comp = 0; comp < numcomponents; comp++) {
		for (int plevel = 0; plevel < nlevels; plevel++) {
			// root filter pyramid level
			int rlevel = plevel + interval;
			const Matrix<double> *mxA = rootscores[comp][rlevel];
			double *rtscore;
			if(mxA == 0) dim[0] = 0, dim[1] = 0;
			else
			{
				dim[0] = mxA->GetDims()[0],dim[1] = mxA->GetDims()[1];
				rtscore = mxA->GetData();
			}
			// process each location in the current pyramid level
			for (int rx = ceil(padx/2.0); rx < dim[1] - ceil(padx/2.0); rx++) {
				for (int ry = ceil(pady/2.0); ry < dim[0] - ceil(pady/2.0); ry++) {
					// get stage 0 score (PCA root + component offset)
					double score = *(rtscore + rx*dim[0] + ry);

					// record score of PCA filter ('score' has the component offset added 
					// to it, so we subtract it here to get just the PCA filter score)
					pcascore[comp][0] = score - offsets[comp];

					// cascade stages 1 through 2*numparts+2
					int stage = 1;
					int numstages = 2*numparts[comp]+2; 
					for (; stage < numstages; stage++) {
						// check for hypothesis pruning
						if (score < t[comp][2*stage-1])
							break;

						// pca == 1 if we're placing pca filters
						// pca == 0 if we're placing "full"/non-pca filters
						int pca = (stage < numparts[comp]+1 ? 1 : 0);
						// get the part# used in this stage
						// root parts have index -1, non-root parts are indexed 0:numparts
						int part = partorder[comp][stage];

						if (part == -1) {
							// we just finished placing all PCA filters, now replace the PCA root
							// filter with the non-PCA root filter
							double rscore = rconv(rlevel, comp, rx, ry, pca);
							score += rscore - pcascore[comp][0];
						} else {
							// place a non-root filter (either PCA or non-PCA)
							int px = 2*rx + (int)anchors[comp][part][0];
							int py = 2*ry + (int)anchors[comp][part][1];
							// lookup the filter and deformation model used by this part
							int filterind = pfind[comp][part];
							int tdefind = defind[comp][part];
							double defthresh = t[comp][2*stage] - score;
							double ps = partscore(plevel, tdefind, filterind, 
								px, py, pca, defthresh);
							if (pca == 1) {
								// record PCA filter score and update hypothesis score with ps
								pcascore[comp][part+1] = ps;
								score += ps;
							} else {
								// update hypothesis score by replacing the PCA filter score with ps
								score += ps - pcascore[comp][part+1];
							}
						}
					}
					// check if the hypothesis passed all stages with a final score over 
					// the global threshold
					//! debug
					if(stage == numstages)
						stage = stage;
					if(score >= thresh)
						stage = stage;
					if (stage == numstages && score >= thresh) {
						// compute and record image coordinates of the detection window
						double scale = sbin/fp1->getScales()[rlevel];
						double x1 = (rx-padx)*scale;
						double y1 = (ry-pady)*scale;
						double x2 = x1 + rootfilterdims[comp][1]*scale - 1;
						double y2 = y1 + rootfilterdims[comp][0]*scale - 1;
						// add 1 for matlab 1-based indexes
						coords.push_back(x1+1);
						coords.push_back(y1+1);
						coords.push_back(x2+1);
						coords.push_back(y2+1);
						// compute and record image coordinates of the part filters
						//! you can open it ,if you need this info
						/*scale = sbin/fp1->getScales()[plevel];
						for (int P = 0; P < numparts[comp]; P++) {
							int probex = 2*rx + (int)anchors[comp][P][0];
							int probey = 2*ry + (int)anchors[comp][P][1];
							int dind = defind[comp][P];
							int offset = LOFFDT[plevel] + dind*featdimsprod[plevel] 
							+ (probex-padx)*featdims[plevel][0] + (probey-pady);
							int px = *(DXAM[0] + offset) + padx;
							int py = *(DYAM[0] + offset) + pady;
							double x1 = (px-2*padx)*scale;
							double y1 = (py-2*pady)*scale;
							double x2 = x1 + partfilterdims[P][1]*scale - 1;
							double y2 = y1 + partfilterdims[P][0]*scale - 1;
							coords.push_back(x1+1);
							coords.push_back(y1+1);
							coords.push_back(x2+1);
							coords.push_back(y2+1);
						}
						// record component number and score
						coords.push_back(comp+1);*/
						coords.push_back(score);
					}
				} // end loop over root y
			}   // end loop over root x
		}     // end loop over pyramid levels
	}         // end loop over components

	Matrix<double> *result = 0;
	int cdims = coords.size();
	if(cdims > 0) {
		result = new Matrix<double>();
		result->Set(1, &cdims, &coords[0]);
	}

	for (int i = 0; i < numcomponents; i++)
		delete [] pcascore[i];
	delete [] pcascore;
	Cleanup();

	return result;
}

double DTDPM::conv(int x, int y, const double *A, const unsigned int *A_dims, const double *B, \
						const unsigned int *B_dims, int num_features) {
	double val = 0;
	const double *A_src = A + x*A_dims[0] + y;
	const double *B_off = B;
	int A_inc = A_dims[0]*A_dims[1];
	for (int f = 0; f < num_features; f++) {
		const double *A_off = A_src;
		for (int xp = 0; xp < B_dims[1]; xp++) {
			// valid only for filters with <= 20 rows
			switch(B_dims[0]) {
				case 20: val += A_off[19] * B_off[19];
				case 19: val += A_off[18] * B_off[18];
				case 18: val += A_off[17] * B_off[17];
				case 17: val += A_off[16] * B_off[16];
				case 16: val += A_off[15] * B_off[15];
				case 15: val += A_off[14] * B_off[14];
				case 14: val += A_off[13] * B_off[13];
				case 13: val += A_off[12] * B_off[12];
				case 12: val += A_off[11] * B_off[11];
				case 11: val += A_off[10] * B_off[10];
				case 10: val += A_off[9] * B_off[9];
				case 9: val += A_off[8] * B_off[8];
				case 8: val += A_off[7] * B_off[7];
				case 7: val += A_off[6] * B_off[6];
				case 6: val += A_off[5] * B_off[5];
				case 5: val += A_off[4] * B_off[4];
				case 4: val += A_off[3] * B_off[3];
				case 3: val += A_off[2] * B_off[2];
				case 2: val += A_off[1] * B_off[1];
				case 1: val += A_off[0] * B_off[0];
			}
			A_off += A_dims[0];
			B_off += B_dims[0];
		}
		A_src += A_inc;
	}
	return val;
}

double DTDPM::rconv(int L, int filterind, int x, int y, int pca)
{
	const unsigned int *A_dims = (const unsigned int*)featdims[L];
	const double *A = feat[pca][L];
	const unsigned int *B_dims = rootfilterdims[filterind];
	const double *B = rootfilters[filterind];
	int num_features = numfeatures;
	// compute convolution
	return conv(x, y, A, A_dims, B, B_dims, num_features);
}

double DTDPM::pconvdt(int L, int probex, int probey, int filterind, int defindex, int xstart, int xend, int ystart, int yend, int pca, double defthresh) {
	const unsigned int *A_dims = (const unsigned int*)featdims[L];
	const double *A = feat[pca][L];
	const unsigned int *B_dims = partfilterdims[filterind];
	const double *B = partfilters[pca][filterind];
	int num_features = (pca == 1 ? pcadim : numfeatures);

	double *ptrbase = PCONV[pca] + LOFFCONV[L] + filterind*featdimsprod[L];
	for (int x = xstart; x <= xend; x++) {
		double *ptr = ptrbase + x*featdims[L][0] + ystart-1;
		for (int y = ystart; y <= yend; y++) {
			ptr++;
			// skip if already computed
			if (*ptr > -INFINITY)
				continue;
			// check for deformation pruning
			double defcost = DXDEFCACHE[defindex][probex-x+S] + DYDEFCACHE[defindex][probey-y+S];
			if (defcost < defthresh)
				continue;
			// compute convolution
			*ptr = conv(x, y, A, A_dims, B, B_dims, num_features);
		}
	}

	// do distance transform over the region.
	// the region is small enough that brute force DT
	// is the fastest method.
	double max = -INFINITY;
	int xargmax = 0;
	int yargmax = 0;

	for (int x = xstart; x <= xend; x++) {
		double *ptr = ptrbase + x*featdims[L][0] + ystart-1;
		for (int y = ystart; y <= yend; y++) {
			ptr++;
			double val = *ptr + DXDEFCACHE[defindex][probex-x+S]
			+ DYDEFCACHE[defindex][probey-y+S];
			if (val > max) {
				max = val;
				xargmax = x;
				yargmax = y;
			}
		}
	}
	int offset = defindex*featdimsprod[L]
	+ probex*featdims[L][0] 
	+ probey;

	// record max and argmax for DT
	*(DXAM[pca] + LOFFDT[L] + offset) = xargmax;
	*(DYAM[pca] + LOFFDT[L] + offset) = yargmax;
	*(DT[pca] + LOFFDT[L] + offset) = max;
	return max;
}

double DTDPM::partscore(int L, int defindex, int pfind, int x, int y, int pca, double defthresh) {
	// remove virtual padding
	x -= padx;
	y -= pady;
	// check if already computed...
	int offset = defindex*featdimsprod[L]
	+ x*featdims[L][0] 
	+ y;
	double *ptr = DT[pca] + LOFFDT[L] + offset;
	if (*ptr > -INFINITY)
		return *ptr;

	// ...nope, define the bounds of the convolution and 
	// distance transform region
	int xstart = x-S;
	xstart = (xstart < 0 ? 0 : xstart);
	int xend = x+S;

	int ystart = y-S;
	ystart = (ystart < 0 ? 0 : ystart);
	int yend = y+S;

	const unsigned int *A_dims = (const unsigned int*)featdims[L];
	const unsigned int *B_dims = partfilterdims[pfind];
	yend = (B_dims[0] + yend > A_dims[0])
		? yend = A_dims[0] - B_dims[0]
	: yend;
	xend = (B_dims[1] + xend > A_dims[1])
		? xend = A_dims[1] - B_dims[1]
	: xend;

	// do convolution and distance transform in region
	// [xstart, xend] x [ystart, yend]
	return pconvdt(L, x, y, pfind, defindex, xstart, xend, ystart, yend, pca, defthresh);
}

void DTDPM::Init(_FeatPyramid *fp1, _FeatPyramid *fp2, int s) {
	numlevels    = fp1->getDim();
	featdims     = new int*[numlevels];
	featdimsprod = new int[numlevels];
	feat[0]      = new double*[numlevels];
	feat[1]      = new double*[numlevels];
	for (int l = 0; l < numlevels; l++) {
		featdims[l] = new int[3];
		Matrix<double> *tmp = fp1->getFeat(l);
		featdims[l][0] = tmp->GetDims()[0];
		featdims[l][1] = tmp->GetDims()[1];
		featdims[l][2] = tmp->GetDims()[2];
		featdimsprod[l]     = featdims[l][0]*featdims[l][1];

		feat[0][l] = tmp->GetData();
		tmp = fp2->getFeat(l);
		feat[1][l] = tmp->GetData();
	}
	numfeatures = fp1->getFeat(0)->GetDims()[2];
	pcadim = fp2->getFeat(0)->GetDims()[2];

	// allocate memory for storing convolution and
	// distance transform data pyramids
	int N    = s;
	PCONV[0] = new double[N];
	PCONV[1] = new double[N];
	DT[0]    = new double[N];
	DT[1]    = new double[N];
	for (int i = 0; i < N; i++) {
		PCONV[0][i] = -INFINITY;
		PCONV[1][i] = -INFINITY;
		DT[0][i]    = -INFINITY;
		DT[1][i]    = -INFINITY;
	}

	// each data pyramid (convolution and distance transform)
	// is stored in a 1D array.  since pyramid levels have
	// different sizes, we build an array of offset values
	// in order to index by level.  the last offset is the
	// total length of the pyramid storage array.
	LOFFCONV = new int[numlevels+1];
	LOFFDT = new int[numlevels+1];
	LOFFCONV[0] = 0;
	LOFFDT[0] = 0;
	for (int i = 1; i < numlevels+1; i++) {
		LOFFCONV[i] = LOFFCONV[i-1] + numpartfilters*featdimsprod[i-1];
		LOFFDT[i]   = LOFFDT[i-1]   + numdefparams*featdimsprod[i-1];
	}

	// cache of precomputed deformation costs
	DXDEFCACHE = new double*[numdefparams];
	DYDEFCACHE = new double*[numdefparams];
	for (int i = 0; i < numdefparams; i++) {
		const double *def = defs[i];
		DXDEFCACHE[i] = new double[2*S+1];
		DYDEFCACHE[i] = new double[2*S+1];
		for (int j = 0; j < 2*S+1; j++) {
			DXDEFCACHE[i][j] = -def[0]*square(j-S) - def[1]*(j-S);
			DYDEFCACHE[i][j] = -def[2]*square(j-S) - def[3]*(j-S);
		}
	}

	for (int p = 0; p < 2; p++) {
		// allocate memory (left uninitialized intentionally)
		DXAM[p] = new int[LOFFDT[numlevels]];
		DYAM[p] = new int[LOFFDT[numlevels]];
	}
}

void DTDPM::Cleanup() {
	delete []featdimsprod;
	delete feat[0];
	delete feat[1];
	for (int l = 0; l < numlevels; l++)
		delete []featdims[l];
	delete []featdims;

	delete []PCONV[0];
	delete []PCONV[1];
	delete []DT[0];
	delete []DT[1];

	delete []LOFFCONV;
	delete []LOFFDT;

	for (int i = 0; i < numdefparams; i++) {
		delete []DXDEFCACHE[i];
		delete []DYDEFCACHE[i];
	}
	delete []DXDEFCACHE;
	delete []DYDEFCACHE;
	for (int p = 0; p < 2; p++) {
		delete []DXAM[p];
		delete []DYAM[p];
	}
};