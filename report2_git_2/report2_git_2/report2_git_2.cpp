#include "stdafx.h"
#include "GA.h"

#define MAX_GENERATION 30000
#define MAX_GENOM_LIST 50
#define VAR_NUM 2
#define __ENABLE_MUTATION__
//min:(-0.54719,-1.54719)=-1.9133

int main()
{
	double varMax[] = { 4.0,4.0 };
	double varMin[] = { -1.5,-3.0 };

	//配列をstd::vectorへ変換
	//std::vector<double> vMax(varMax, std::end(varMax));
	//std::vector<double> vMin(varMin, std::end(varMin));
	GA ga(MAX_GENOM_LIST, VAR_NUM, std::vector<double>(varMax, std::end(varMax)) , std::vector<double>(varMin, std::end(varMin)));//遺伝的アルゴリズム諸関数をまとめたクラスの宣言

	for (int i = 0; i <= MAX_GENERATION; i++)//メインのループ
	{
		bool change = ga.selection();//選択

		ga.blxAlphaCrossover();//交叉
#ifdef __ENABLE_MUTATION__
		ga.mutation();//突然変異
#endif
		if (i % (MAX_GENERATION / 10) == 0 || change)
		{
			std::cout << "i=" << std::to_string(i) << std::endl;
			ga.calc(true, i % (MAX_GENERATION / 10) != 0);//評価関数の計算
		}
		else
		{
			ga.calc(false);//評価関数の計算
		}
	}

	while (1)
	{
		if (_kbhit() && _getch() == 27)
			break;
	}
	return 0;
}

