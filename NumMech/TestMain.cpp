#include "Matrix.h"

using namespace std;

template <typename T>
void print(vector<T> vec);

int main()
{
	//cin >> input;

	vector<int> v1{ 1,2,3 };
	vector<int> v2{ 4,5,6 };
	vector<int> v3{ 7,8,9 };
	vector<vector<int>> M{ v1, v2, v3 };

	vector<vector<int>> M2{ {1,1,1}, {1,1,1}, {1,1,1} };

	function<double(double)> f{ [](double x) {return x * x * x; } };

	Matrix<int> mat(M);
	Matrix<int> mat2(M2);

	//mat = mat2;
	mat2 = mat.Transpose();

	mat.col(2) = vector<int>{ 1,1,1 };
	//mat.transpose();

	std::cout << '\n' << std::endl;

	//mat.funcToRow(f, 2);

	//print(init[0]);

	Matrix<int> mat3 = mat - mat2;


	mat.print();

	//cin >> input;
}
template <typename T>
void print(vector<T> vec)
{
	for (T el : vec)
		cout << el << ' ';
}