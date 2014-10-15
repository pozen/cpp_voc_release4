#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include "dtdpm.h"

struct Rect {
	int x;
	int y;
	int height;
	int width;
};

Matrix<double>* FitObjDetectionData(Rect rect, IplImage* img) {
	Matrix<double>* mat = new Matrix<double>();
	int dims[3];
	dims[0] = rect.height;
	dims[1] = rect.width;
	dims[2] = 1;
	mat->Set(3, dims, 0);
	double *tmpData = mat->GetData();

	int count = 0;

	for ( int i = rect.x; i < rect.x + rect.width; i++ ) {
		for ( int j = rect.y ; j < rect.y + rect.height; j++ ) {
			tmpData[count] = (unsigned char)img->imageData[j*img->widthStep+i];
			count++;
		}
	}
	return mat;
}

int main() {
	IplImage* src = cvLoadImage("16.jpg");
	IplImage* gray = cvCreateImage(cvGetSize(src), 8, 1);
	cvCvtColor(src, gray, CV_BGR2GRAY);

	DTDPM* dtdpm = new DTDPM("car.model");

	Rect testRect;
	testRect.x = testRect.y = 0;
	testRect.height = gray->height;
	testRect.width = gray->width;
	Matrix<double>* mat = FitObjDetectionData(testRect, gray);
	_FeatPyramid *pyr = new _FeatPyramid(mat, dtdpm->GetMaxSize(), -1, -1, dtdpm->GetMaxSbin(), dtdpm->GetInterval() );

	ObjBox *result;
	int pickNum;

	dtdpm->CascadeDetect(-0.5, 0.5, pyr, &result, &pickNum);

	for (int j = 0; j < pickNum; ++j){
		cvRectangle(src, cvPoint(result[j].x1, result[j].y1),
			cvPoint(result[j].x2, result[j].y2), \
			CV_RGB( 255, 255, 255),  1, 8,  0);
	}
	cvShowImage("test", src);
	cvWaitKey(0);
	delete mat;
	delete pyr;
	delete result;
	return 0;
}