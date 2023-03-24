#pragma once
#include <vector>

using std::vector;

template <typename T>
class Matrix {
private:
	int rows;
	int cols;
	vector<vector<T>> data;
public:
	Matrix();
	Matrix(int size);
	Matrix(int rows, int cols);
	Matrix(int rows, int cols, T defaultVal);
	Matrix(vector<vector<T>> data, int rows, int cols);
	Matrix(T**& data, int rows, int cols);
	Matrix(const Matrix& mat);
	~Matrix();

	const int getRows();
	const int getCols();

	const T get(int row, int col);
	void set(int row, int col, T value);

	vector<T>& operator[](int index);

	Matrix<T>& operator=(const Matrix<T>& m);
	Matrix<T>& operator+=(const Matrix<T>& m);
	Matrix<T>& operator-=(const Matrix<T>& m);
	Matrix<T>& operator*=(const Matrix<T>& m);
	Matrix<T>& operator*=(T val);
	Matrix<T>& operator/=(T val);

	template <typename V>
	friend Matrix<T> operator+(const Matrix<T>& m1, const Matrix<T>& m2);
	template <typename V>
	friend Matrix<T> operator-(const Matrix<T>& m1, const Matrix<T>& m2);
	template <typename V>
	friend Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2);
	template <typename V>
	friend Matrix<T> operator*(const Matrix<T>& m, T val);
	template <typename V>
	friend Matrix<T> operator*(T val, const Matrix<T>& m);
	template <typename V>
	friend Matrix<T> operator/(const Matrix<T>& m, T val);
};

template<typename T>
Matrix<T>::Matrix() : rows(0), cols(0) {
}

template<typename T>
Matrix<T>::Matrix(int size) : Matrix(size, size) {
}

template<typename T>
Matrix<T>::Matrix(int rows, int cols) : rows(rows), cols(cols) {
	data.resize(rows, vector<T>(cols));
}

template<typename T>
Matrix<T>::Matrix(int rows, int cols, T defaultVal) : rows(rows), cols(cols) {
	data.resize(rows, vector<T>(cols));
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			data[i][j] = defaultVal;
		}
	}
}

template<typename T>
Matrix<T>::Matrix(vector<vector<T>> data, int rows, int cols) : rows(rows), cols(cols) {
	this->data.resize(rows, vector<T>(cols));
	this->data = data;
}

template<typename T>
Matrix<T>::Matrix(T**& data, int rows, int cols) : rows(rows), cols(cols) {
	this->data.resize(rows, vector<T>(cols));
	try {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				this->data[i][j] = data[i][j];
			}
		}
	}
	catch (const std::exception&) {
		this->data.resize(rows, vector<T>(cols));
	}
}

template<typename T>
Matrix<T>::Matrix(const Matrix& mat) : rows(mat.rows), cols(mat.cols) {
	this->data = mat.data;
}

template<typename T>
Matrix<T>::~Matrix() {
	this->data.clear();
}

template<typename T>
const int Matrix<T>::getRows() {
	return this->rows;
}

template<typename T>
const int Matrix<T>::getCols() {
	return this->cols;
}

template<typename T>
const T Matrix<T>::get(int row, int col) {
	return this->data[row][col];
}

template<typename T>
void Matrix<T>::set(int row, int col, T value) {
	this->data[row][col] = value;
}

template<typename T>
vector<T>& Matrix<T>::operator[](int index) {
	return this->data[index];
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m) {
	if (this == &m) {
		return *this;
	}
	this->cols = m.cols;
	this->rows = m.rows;
	this->data.resize(rows, vector<T>(cols));
	this->data = m.data;
	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m) {
	if (this->cols != m.cols || this->rows != m.rows) {
		//throw std::invalid_argument("Matrices must have the same dimensions.");
		return *this;
	}
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			data[i][j] += m.data[i][j];
		}
	}
	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m) {
	if (this->cols != m.cols || this->rows != m.rows) {
		//throw std::invalid_argument("Matrices must have the same dimensions.");
		return *this;
	}
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			this->data[i][j] -= m.data[i][j];
		}
	}
	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& m) {
	if (this->cols != m.rows) {
		//throw std::invalid_argument("The number of columns in the first matrix must match the number of rows in the second matrix.");
		return *this;
	}
	Matrix<T> temp(this->rows, m.cols);
	for (int i = 0; i < temp.rows; ++i) {
		for (int j = 0; j < temp.cols; ++j) {
			for (int k = 0; k < this->cols; ++k) {
				temp.data[i][j] += (this->data[i][k] * m.data[k][j]);
			}
		}
	}
	return (*this = temp);
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(T val) {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			this->data[i][j] *= val;
		}
	}
	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator/=(T val) {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			this->data[i][j] /= val;
		}
	}
	return *this;
}


template<typename T>
Matrix<T> operator+(const Matrix<T>& m1, const Matrix<T>& m2) {
	Matrix<T> temp(m1);
	return (temp += m2);
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& m1, const Matrix<T>& m2) {
	Matrix<T> temp(m1);
	return (temp -= m2);
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2) {
	Matrix<T> temp(m1);
	return (temp *= m2);
}
template
Matrix<double> operator*(const Matrix<double>& m1, const Matrix<double>& m2);

template<typename T>
Matrix<T> operator*(const Matrix<T>& m, T val) {
	Matrix<T> temp(m);
	return (temp *= val);
}

template<typename T>
Matrix<T> operator*(T val, const Matrix<T>& m) {
	return m * val;
}

template<typename T>
Matrix<T> operator/(const Matrix<T>& m, T val) {
	Matrix<T> temp(m);
	return (temp /= val);
}