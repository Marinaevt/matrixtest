#include "SLAE.h"

CSLAE::CSLAE() {
	m_is_null = true;
}
void CSLAE::ZeroAll()
{
	m_matr.reset(0.0);
	for (size_t i = 0; i < m_rp.size(); i++) {
		m_rp[i] = m_sol[i] = 0;
	}
}
bool CSLAE::Init(size_t eqns, size_t band)
{
	if (!m_matr.resize(eqns, band, 0))
		return false;

//	try {
		m_rp.resize(eqns);
		m_sol.resize(eqns);
		m_flg.resize(eqns);
//	}
//	catch (CException* pEx) {
//		CDlgShowError cError(ID_ERROR_SLAE_INIT);	//_T("SLAE Init error"));
//		pEx->Delete();
//		return false;
//	}
	//TRY & END_CATCH_ALL	//избавляемся от MFC макросов и вообще от MFC по максимуму

	ZeroAll();

	return true;
}

void CSLAE::Gauss(/*int nxy2, int isl, bool bZZ*/)
{
	// nxy2 - кол-во уравнений
	// isl - ширина ленты (половины)

	size_t r, s, m, n, j;	// индексы
	double zn, anul;

	size_t nxy2 = m_matr.rows();
	size_t isl = m_matr.band();	//половина, так половина

	for (r = 0; r < nxy2; r++)
	{
		m_rp[r] /= m_matr.cell(r, 0);

		if (r == nxy2 - 1)		//?????
			break;

		zn = m_matr.cell(r, 0);	// ERROR_&_CRASH
		if (fabs(zn) < EPS) return;		// IF CRASH

		for (s = 1; s < isl; s++)
		{
			m_sol[s] = m_matr.cell(r, s);

			if (fabs(m_sol[s]) < EPS)
				continue;

			m_matr.cell(r, s) = m_sol[s] / zn;
		}

		for (m = 1; m < isl; m++)
		{
			zn = m_sol[m];

			if (fabs(zn) < EPS)
				continue;

			n = r + m;

			if (n > nxy2 - 1)
				continue;

			j = 0;
			for (s = m; s < isl; s++)
			{
				anul = m_matr.cell(n, j);
				m_matr.cell(n, j) -= zn * m_matr.cell(r, s);
				if (fabs(m_matr.cell(n, 0)) < EPS) {
					m_matr.cell(n, 0) = EPS; //10^(-18)
				}
				j++;
			}
			m_rp[n] -= zn * m_rp[r];
		}
	}

	// цикл вычисления решения
	for (r = nxy2 - 2; r >= 0; r--) {
		for (s = 1; s < isl; s++) {
			m = r + s;

			if (m > nxy2 - 1)
				continue;

			m_rp[r] -= m_matr.cell(r, s) * m_rp[m];
		}
	}

	// сохраняем решение в m_sol
	for (r = 0; r < nxy2; r++) {
		//Math::swap(m_sol[r], m_rp[r]);
		m_sol[r] = m_rp[r];
	}
}
void CSLAE::RotateRPLCS(size_t k, DBL alpha)
{
	DBL sina = sin(alpha),
		cosa = cos(alpha);

	DBL rp1 = m_rp[k] * cosa + m_rp[k + 1] * sina,
		rp2 = m_rp[k] * (-sina) + m_rp[k + 1] * cosa;

	m_rp[k] = rp1;
	m_rp[k + 1] = rp2;
}

void CSLAE::RotateMatrixLCS(size_t k, DBL alpha) //k - номер узла
{
	DBL sina = sin(alpha), cosa = cos(alpha);
	for (size_t r = (k >= m_matr.band()? k - m_matr.band() :0); r < k; r++) // перебор элементов столбца k
	{
		size_t c = k - r; // индекс столбца k
		DBL ai1, ai2;
		ai1 = m_matr.direct_cell(r, c) * cosa + m_matr.direct_cell(r, c + 1) * sina;
		ai2 = -m_matr.direct_cell(r, c) * sina + m_matr.direct_cell(r, c + 1) * cosa;

		m_matr.direct_cell(r, c) = ai1;
		m_matr.direct_cell(r, c + 1) = ai2;
	}
	for (size_t c = 2; c < m_matr.band(); c++)
	{
		DBL a1j, a2j;

		a1j = m_matr.direct_cell(k, c) * cosa + m_matr.direct_cell(k + 1, c - 1) * sina;
		a2j = -m_matr.direct_cell(k, c) * sina + m_matr.direct_cell(k + 1, c - 1) * cosa;

		m_matr.direct_cell(k, c) = a1j;
		m_matr.direct_cell(k + 1, c - 1) = a2j;
	}
	DBL akk = m_matr.cell(k, k) * cosa * cosa
		+ 2 * m_matr.cell(k, k + 1) * sina * cosa + m_matr.cell(k + 1, k + 1)* sina * sina;
	DBL akk1 = (m_matr.cell(k + 1, k + 1) - m_matr.cell(k, k))*sina * cosa
		+ m_matr.cell(k, k + 1) * (cosa * cosa - sina * sina);
	DBL ak1k1 = m_matr.cell(k + 1, k + 1) * cosa * cosa
		- 2 * m_matr.cell(k, k + 1) * sina * cosa + m_matr.cell(k, k)* sina * sina;
	m_matr.cell(k, k) = akk;
	m_matr.cell(k, k + 1) = akk1;
	m_matr.cell(k + 1, k + 1) = ak1k1;
}

int main() {
	// Дескриптор DLL-библиотеки
	HMODULE hDll;
	// Указатель на функцию
	int(*dllgauss) (double*, double*, double*, double*, int, int);

	// Загружаем динамически подключаемую библиотеку
	hDll = LoadLibraryEx(_T("ImprovedSystem.dll"), 0, DONT_RESOLVE_DLL_REFERENCES);

	if (!hDll)
	{
		cout << _T("Динамическая библиотека не загружена") << endl;
		return GetLastError();
	}
	dllgauss = (int(*)(double*, double*, double*, double*, int, int))GetProcAddress(hDll, "ImprovedGaussSystem");
	if (!dllgauss)
	{
		cout << _T("Ошибка получения адреса функции") << endl;
		return GetLastError();
	}

	CSLAE test;
	
	int n=533;
	double* tt1 = new double [n];
	double* tt2 = new double[n];
	test.Init(n, n);

	ifstream Kmat("Kbig.txt");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Kmat >> test.m_matr.cell(i, j);
		}
	}
	Kmat.close();
	ifstream Fv("Fbig.txt");
	for (int i = 0; i < n; i++) {
		Fv >> test.m_rp[i];
	}
	/*
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << test.m_matr.cell(i, j) << '\t';
		}
		cout << endl;
	}
	*/
	Fv.close();
	/*
	for (int i = 0; i < n; i++) {
		cout << test.m_rp[i] << endl;
	}
	test.RotateMatrixLCS(4, 45*MPI/180);
	
	ifstream Kmatr("Krot.txt");

	double z;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Kmatr >> z;
			cout << z - test.m_matr.cell(i, j) << '\t';
		}
		cout << endl;
	}
	Kmatr.close();
	*/
	//test.Gauss();
	dllgauss(&test.m_matr[0], &test.m_rp[0], tt2, tt1, n, n/2);
	ifstream Uv("Ubig.txt");
	double k;
	cout << "solve" << endl;
	for (int i = 0; i < n; i++) {
		Uv >> k;
		cout << k << "\t" << test.m_rp[i] << endl;
	}
	Uv.close();
	/*
	test.RotateRPLCS(4,-45 * MPI / 180);

	cout << "Rsolve" << endl;
	ifstream Urotv("Urot.txt");
	for (int i = 0; i < n; i++) {
		Urotv >> k;
		cout << k << '\t' << test.m_sol[i] << endl;
	}
	Urotv.close();
	*/
	// Отключаем библиотеку
	if (!FreeLibrary(hDll))
	{
		cout << _T("Ошибка выгрузки библиотеки из памяти") << endl;
		return GetLastError();
	}

}