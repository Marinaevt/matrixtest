#pragma once
// стандартна€ ленточна€ матрица
#define BANDMATRIX CBandMatrix<DBL>
#define DBL double
#define MPI 3.141592653589793
#define EPS 0.000000000000001
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#include <vector>
#include <iostream>
#include <fstream>
#include <Windows.h>
#include <tchar.h>

using namespace std;


template <typename SizeT, typename Cell>
class CCustomArray
{
protected:

	Cell* data_;
	SizeT size_;

public:

	typedef Cell Cell;
	typedef SizeT SizeT;


	CCustomArray()
	{
		data_ = nullptr;
		size_ = 0;
	}

	CCustomArray(SizeT size, Cell* data = nullptr)
	{
		if (data_ && size) {
			data_ = nullptr;
			resize(size);
		}
	}


	CCustomArray(const CCustomArray& array)
	{
		if (data_ != array.data_) {
			data_ = nullptr;
			resize(array.size_);
			memcpy(data_, array.data_, size_ * sizeof(Cell));
		}
	}

	bool resize(SizeT newSize)
	{
		if (newSize == size_) {
			return true;
		}

		free();

		if (newSize) {
			data_ = new Cell[newSize];
			//if (data_) {
			size_ = newSize;
			return true;
			//}
			//return false;
		}

		return false;
	}

	void free()
	{
		//if (data_) 
		delete[] data_;
		//data_ = nullptr;
		size_ = 0;
	}

	void swap(CCustomArray& array)
	{
		Cell* tempData = data_;
		SizeT tempSize = size_;

		data_ = array.data_;
		size_ = array.size_;

		array.data_ = tempData;
		array.size_ = tempSize;
	}

	SizeT size() const
	{
		return size_;
	}

	Cell* data()
	{
		return data_;
	}

	const Cell* data() const
	{
		return data_;
	}

	Cell& operator[] (SizeT i)
	{
		return data_[i];
	}

	const Cell& operator[] (SizeT i) const
	{
		return data_[i];
	}

	CCustomArray& operator= (const CCustomArray& array)
	{
		if (data_ != array.data_) {
			data_ = nullptr;
			resize(array.size_);
			memcpy(data_, array.data_, size_ * sizeof(Cell));
		}
		return *this;
	}

	// индексирование 2мерной матрицы, когда заранее известно
	// конкретное число столбцов
	template <SizeT Cols>
	Cell& cell2d(SizeT i, SizeT j)
	{
		return data_[i * Cols + j];
	}

	template <SizeT Cols>
	const Cell& cell2d(SizeT i, SizeT j) const
	{
		return data_[i * Cols + j];
	}

	// индексирование 2мерной матрицы
	Cell& cell2d(SizeT i, SizeT j, SizeT cols)
	{
		return data_[i * cols + j];
	}

	const Cell& cell2d(SizeT i, SizeT j, SizeT cols) const
	{
		return data_[i * cols + j];
	}

	~CCustomArray()
	{
		free();
	}
};


// symmetric band matrix class
// stores only upper-half of band (with diagonal, of course)
template <typename Cell>
class CBandMatrix : public CCustomArray<size_t, Cell>
{
	static Cell zero_; // остальные €чейки

protected:

	size_t rows_; // кол-во строк
	size_t band_; // размер полосы

public:

	typedef CCustomArray<size_t, Cell> Array;
	typedef Cell Cell;

	CBandMatrix()
	{
		rows_ = band_ = 0;
	}

	CBandMatrix(size_t rows, size_t band) : Array(rows * band),
		rows_(rows),
		band_(band)
	{}

	CBandMatrix(const CBandMatrix& matrix) : CCustomArray(matrix),
		rows_(matrix.rows_),
		band_(matrix.band_)
	{}

	bool resize(size_t rows = 0, size_t band = 0, Cell value = Cell(0))
	{
		if (rows && band) {
			if (Array::resize(rows * band))
			{
				rows_ = rows;
				band_ = band;
				for (size_t i = 0; i < size_; i++)
					data_[i] = value;
				return true;
			}
			return false;

		}
		else {
			free();
		}
		return true;
	}

	void reset(const Cell& c = Cell(0))
	{
		for (SizeT i = 0; i < size_; i++)
			data_[i] = c;
	}

	//обнуление строк и полосы
	void free()
	{
		Array::free();
		rows_ = band_ = 0;
	}

	void swap(CBandMatrix& matrix)
	{
		Array::swap(matrix);

		rows_ = matrix.rows_;
		band_ = matrix.band_;
	}

	//возвращает кол-во строк
	size_t rows() const
	{
		return rows_;
	}

	//!возвращает размер полосы (или размер матрицы?)
	size_t band() const
	{
		return band_;
	}

	//возвращает кол-во строк (???)
	size_t cols() const
	{
		return rows();
	}

	///! пр€мое индексирование
	Cell& direct_cell(size_t i, size_t j)
	{
		return cell2d(i, j, band_);
	}

	///! пр€мое индексирование
	const Cell& direct_cell(size_t i, size_t j) const
	{
		return cell2d(i, j, band_);
	}

	///! универсальное ("безопасное") индексирование
	Cell& cell(size_t i, size_t j)
	{
		if ((j >= i) && (j - i < band_))
			return cell2d(i, j - i, band_);

		if ((j < i) && (i - j < band_))
			return cell2d(j, i - j, band_);

		return zero_;
	}

	CBandMatrix& operator= (const CBandMatrix& matrix)
	{
		if (data_ != matrix.data_) {
			Array::operator=(matrix);
			rows_ = rows;
			band_ = band;
		}
		return *this;
	}

	//записывает ленточную (Ќ≈ полную) матрицу
	// true - выводить значение €чейки, false - выводить *, если не 0, иначе - 0
	void WriteToLog(bool val) {

		DLOG(CString(_T("Matrix log start ->")), log_info);
		for (SizeT i = 0; i < rows(); i++) {
			CString tmp = _T("");
			//for (SizeT j = 0; j < cols(); j++)
			for (SizeT j = i; j < min(cols(), i + band_); j++) {

				if (val) {
					tmp += AllToString(cell(i, j)) + CString(_T(","));
				}
				else {
					if (cell(i, j)) {
						tmp += AllToString((int)(min(1, fabs(cell(i, j)) + 1))); //_T("*");
					}
					else {
						tmp += _T("0");
					}
				}
			}

			CString count = _T("0");
			if (i>9) {
				count.Empty();	//если больше 9 строк, то перестаЄм вставл€ть 0: "01, 02,...09, _10"
			}

			DLOG(CString(_T("[")) + count + AllToString(i) + CString(_T("]:")) + tmp, log_info);
		}
		DLOG(CString(_T("Matrix log end   <-")), log_info);
	}

	//записывает квадратную (полную) матрицу ѕ–» ќ“Ћјƒ ≈
	// true - выводить значение €чейки, false - выводить *, если не 0, иначе - 0
	void WriteToLogFullMatrix(bool val) {

		LOG(CString(_T("Matrix log start ->")), log_info);
		size_t row_size = rows();
		size_t col_size = cols();

		for (SizeT i = 0; i < row_size; i++) {
			CString tmp = _T("");
			for (SizeT j = 0; j < col_size; j++) {

				if (val) {
					tmp += AllToString(cell(i, j)) + CString(_T(","));
				}
				else {
					if (cell(i, j)) {
						tmp += AllToString((int)(min(1, fabs(cell(i, j)) + 1))); //_T("*");
					}
					else {
						tmp += _T("0");
					}
				}
			}

			CString count = _T("0");
			if (i>9) {
				count.Empty();	//если больше 9 строк, то перестаЄм вставл€ть 0: "01, 02,...09, _10"
			}

			LOG(CString(_T("[")) + count + AllToString(i) + CString(_T("]:")) + tmp, log_info);
		}
		LOG(CString(_T("Matrix log end   <-")), log_info);
	}

	//записывает квадратную (полную) матрицу ѕ–» ќ“Ћјƒ ≈
	// true - выводить значение €чейки, false - выводить *, если не 0, иначе - 0
	void WriteToLogFullMatrixOnDebug(bool val) {

		DLOG(CString(_T("Matrix log start ->")), log_info);
		size_t row_size = rows();
		size_t col_size = cols();

		for (SizeT i = 0; i < row_size; i++) {
			CString tmp = _T("");
			for (SizeT j = 0; j < col_size; j++) {

				if (val) {
					tmp += AllToString(cell(i, j)) + CString(_T(","));
				}
				else {
					if (cell(i, j)) {
						tmp += AllToString((int)(min(1, fabs(cell(i, j)) + 1))); //_T("*");
					}
					else {
						tmp += _T("0");
					}
				}
			}

			CString count = _T("0");
			if (i>9) {
				count.Empty();	//если больше 9 строк, то перестаЄм вставл€ть 0: "01, 02,...09, _10"
			}

			DLOG(CString(_T("[")) + count + AllToString(i) + CString(_T("]:")) + tmp, log_info);
		}
		DLOG(CString(_T("Matrix log end   <-")), log_info);
	}
	/*
	void WriteToLog()
	{
	DLOG(CString(_T("BandMatrix log start ->")), log_info);
	for (int i = 0; i < rows(); i++)
	{
	CString tmp = _T("");
	for (int j = 0; j < rows(); j++)
	tmp += AllToString(cell(i, j)) + _T(",");

	DLOG(CString(_T("[")) + AllToString(i) + CString(_T("]:")) + tmp, log_info);
	}
	DLOG(CString(_T("BandMatrix log end   <-")), log_info);
	}
	*/
};

template <typename Cell>
Cell CBandMatrix<Cell>::zero_ = Cell(0);



class CSLAE
{
public:
	BANDMATRIX			m_matr; // матрица CBandMatrix DBL
	std::vector<DBL>	m_rp;   // права€ часть Ax=b
	std::vector<DBL>	m_sol;  // решение Ax=b
	std::vector<bool>	m_flg;  // флажки дл€ контрол€ граничных условий	//Ќ≈ »—ѕќЋ№«”ё“—я
	bool			m_is_null;	// true, если только создана и не заполнена

	CSLAE();

	bool Init(size_t eqns, size_t band);
	void ZeroAll();

	// вращает локальную систему координат, заданную переменными k (x_k) и k+1 (y_k), на угол alpha
	// дл€ матрицы
	void RotateMatrixLCS(size_t k, DBL alpha);
	// дл€ правой части
	void RotateRPLCS(size_t k, DBL alpha);
	void Gauss();
	//~CSLAE();
};


