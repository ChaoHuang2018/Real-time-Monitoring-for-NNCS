#include "../flowstar-template/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	clock_t begin, end;
	

	Matrix<Interval> IM1(50, 50), IM2(50, 50);
	Interval I(-2,2);

	
	for(int i=0; i<50; ++i)
	{
		for(int j=0; i<50; ++i)
		{
				IM1[i][j] = I;
				IM2[i][j] = I;
		}
	}
	
	begin = clock();

	Matrix<Interval> IM3 = IM1 * IM2;


	end = clock();
	
	cout << "Computation Time: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

	return 0;
}

