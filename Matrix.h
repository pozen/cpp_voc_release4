#ifndef _MATRIX_H
#define _MATRIX_H

#include "stdio.h"

template<class T>
class Matrix {
public:
	Matrix() { _data = 0; _dims = 0; _numDims = 0; _memFlag = 0; }
	Matrix(int numDims, int *dims, void *mem = 0);
	Matrix(const Matrix& mat);
	~Matrix();
public:
	void Add(const T& val);
	static void Add(Matrix *mat, const T& val);
	void Mul(const T& val);
	static void Mul(Matrix *mat, const T& val);
	void Mul(Matrix *mat, Matrix *dst);
	static void Mul(Matrix *mat1, Matrix *mat2, Matrix *dst);
	void ReShape(const int& numDims, const int *dims);
	static void ReShape(Matrix *mat, const int& numDims, const int *dims);
public:
	//! note! all the stupid methods , just for speed reason here!!!
	T* GetData() const { return _data; }
	int* GetDims() const { return _dims; }
	int  GetDataSize() const { return _dsize; }
	int  GetDimsNum() const { return _numDims; }
	void SetData(const int& index, const T& val) { _data[index] = val; }
	void SetDim(const int& index, const int& size) { _dims[index] = size; }
	void Set(const int &numDims, const int *dims, const T *data, void* mem = 0);
private:
	T *_data;
	int *_dims;
	int _numDims;
	int _memFlag;
	int _dsize;
};

template<class T>
Matrix<T>::Matrix(int numDims, int *dims, void *mem) {
	_numDims = numDims;
	if(mem) {
		_memFlag = 1;
		_dims = (int*)mem;
		mem = (void*)(_dims + numDims);
		_data = (T*)mem;
		_dsize = 1;
		for(int i = 0; i < numDims; i++)
			_dims[i] = dims[i], _dsize *= dims[i];
	} else {
		_memFlag = 0;
		_dims = new int[numDims];
		_dsize = 1;
		for(int i = 0; i < numDims; i++)
			_dims[i] = dims[i], _dsize *= dims[i];
		_data = new T[_dsize];
	}
	for(int i = 0; i < _dsize; i++) _data[i] = 0;
}

template<class T>
Matrix<T>::Matrix(const Matrix<T> &mat) {
	_dsize = mat.GetDataSize();
	_numDims = mat.GetDimsNum();
	_memFlag = 0;
	_dims = new int[_numDims];
	for(int i = 0; i < _numDims; i++) _dims[i] = mat.GetDims()[i];
	_data = new T[_dsize];
	memcpy(_data, mat.GetData(), _dsize*sizeof(T));
}

template<class T>
Matrix<T>::~Matrix() {
	if(_dims) delete []_dims;
	if(_data) delete []_data;
}

template<class T>
void Matrix<T>::Add(const T &val) {
	for(int i = 0; i < _dsize; i++)
		_data[i] += val;
}

template<class T>
void Matrix<T>::Add(Matrix<T> *mat, const T &val) {
	T *data = mat->GetData();
	int datasize = mat->GetDataSize();
	for(int i = 0; i < datasize; i++)
		data[i] += val;
}

template<class T>
void Matrix<T>::Mul(const T &val) {
	for(int i = 0; i < dsize; i++)
		_data[i] *= val;
}

template<class T>
void Matrix<T>::Mul(Matrix<T> *mat, const T &val) {
	T *data = mat->GetData();
	int datasize = mat->GetDataSize();
	for(int i = 0; i < datasize; i++)
		data[i] += val;
}

template<class T>
void Matrix<T>::Mul(Matrix<T> *mat, Matrix<T> *dst) {
	assert(_numDims == 2);
	assert(mat->GetDimsNum() == 2);
	assert(dst->GetDimsNum() == 2);

	int *matdims = mat->GetDims();
	int *dstdims = mat->GetDims();

	assert(_dims[1] == matdims[0]);
	assert(_dims[0] == dstdims[0]);
	assert(matdims[1] == dstdims[1]);
	T *matdata = mat->GetData();
	T *dstdata = dst->GetData();

	for(int i = 0; i < _dims[0]; i++) {
		for(int j = 0; j < matdims[1]; j++) {
			int index = j*dstdims[0] + i;
			dstdata[index] = 0;
			for(int k = 0; k < _dims[1]; k++)
				dstdata[index] += _data[i+k*_dims[0]]*matdata[k+j*matdims[0]];
		}
	}
}

template<class T>
void Matrix<T>::Mul(Matrix<T> *mat1, Matrix<T> *mat2, Matrix<T> *dst) {
	int *mat1dims = mat1->GetDims();
	int *mat2dims = mat2->GetDims();
	int *dstdims = dst->GetDims();

	T *mat1data = mat1->GetData();
	T *mat2data = mat2->GetData();
	T *dstdata = dst->GetData();

	for(int i = 0; i < mat1dims[0]; i++) {
		for(int j = 0; j < mat2dims[1]; j++) {
			int index = j*dstdims[0] + i;
			dstdata[index] = 0;
			for(int k = 0; k < mat1dims[1]; k++)
				dstdata[index] += mat1data[i+k*mat1dims[0]]*mat2data[k+j*mat2dims[0]];
		}
	}
}

template<class T>
void Matrix<T>::ReShape(const int& numDims, const int *dims) {
	int tmpSize = 1;
	for(int i = 0; i < numDims; i++)
		tmpSize *= dims[i];

	_numDims = numDims;
	if(_dims) delete []_dims;
	_dims = new int[_numDims];
	for(int i = 0; i < numDims; i++)
		_dims[i] = dims[i];
}

template<class T>
void Matrix<T>::Set(const int &numDims, const int *dims, const T *data, void* mem) {
	_numDims = numDims;
	if(_dims && !_memFlag) delete []_dims;

	_memFlag = 0;
	if(mem) {
		_dims = (int*)mem;
		mem = (void*)(_dims + numDims);
		_memFlag = 1;
	} else _dims = new int[numDims];

	_dsize = 1;
	for(int i = 0; i < numDims; i++)
		_dims[i] = dims[i], _dsize *= dims[i];
	if(_data && !_memFlag) delete []_data;

	_memFlag = 0;
	if(mem) {
		_data = (T*)mem;
		_memFlag = 1;
	} else _data = new T[_dsize];
	if(data) memcpy(_data, data, _dsize*sizeof(T));
	else {
		for(int i = 0; i < _dsize; i++)
			_data[i] = 0.0;
	}
} 
#endif