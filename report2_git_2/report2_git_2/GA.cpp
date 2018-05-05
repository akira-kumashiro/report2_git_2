#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _var_num, std::vector<double> _varMax, std::vector<double> _varMin) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//dataの初期化
	eliteData(_var_num)
{
	//もらった変数をクラス内変数に格納
	max_genom_list = _max_genom_list;
	var_num = _var_num;
	varMax = _varMax;
	varMin = _varMin;

	for (int i = 0; i < max_genom_list; i++)
	{
		for (int j = 0; j < var_num; j++)
		{
			data[i].x[j] = random(varMin[j], varMax[j]);//遺伝子の初期設定
		}
	}
	prev_data = data;
	calcResult();

	displayValues(false);
}

bool GA::selection()
{
	int max_num = 0;//最も評価の良い個体の番号
	bool ret = isChanged;//最も評価の良い個体の変化の監視(デバッグ用)
	isChanged = false;

	/*for (int i = 0; i < max_genom_list; i++)//ルーレット選択用に評価関数の合計と一番評価の良い番号を取得
	{
		if (data[i].result > data[max_num].result)
			max_num = i;
	}

	eliteData = data[max_num];*/
	eliteData = searchRank(0);//最も評価の良い個体を保持
	prev_data = data;
	for (int i = 0; i < max_genom_list; i++)
	{
		double selector = random(0.0, 1.0);//乱数を生成
		double needle = 0;//ルーレットの針を生成
		int j = 0;
		for (; ; j++)
		{
			needle += (prev_data[j].result / resultSumValue);//ルーレットの針を乱数の値まで進める
			if (needle > selector)
				break;
			if (j == (max_genom_list - 1))
				break;
		}
		data[i] = prev_data[j];
	}
	return ret;
}

bool GA::blxAlphaCrossover()
{
	prev_data = data;

	for (int i = 0; i < max_genom_list; i += 2)//2個ずつ交叉
	{
		for (int j = 0; j < var_num; j++)
		{
			double ave = (data[i].x[j] + data[i + 1].x[j]) / 2;
			double length = std::abs((data[i].x[j] - data[i + 1].x[j]));

			data[i].x[j] = random(ave - length * (1 + alpha * 2) / 2, ave + length * (1 + alpha * 2) / 2);
			data[i + 1].x[j] = random(ave - length * (1 + alpha * 2) / 2, ave + length * (1 + alpha * 2) / 2);
		}
	}
	return true;
}

bool GA::mutation()
{
	for (int i = 0; i < max_genom_list; i++)
	{
		if (random(0.0, 1.0) <= individualMutationRate)//個体突然変異率の計算
		{
			int pos = random(0, var_num - 1);
			data[i].x[pos] += random(varMin[pos] * random(0, genomMutationRate), varMax[pos] * random(0, genomMutationRate));

			//			for (int j = 0; j < var_num; j++)
			//			{
			//				data[i].x[j] = random(varMin[j], varMax[j]);
			//			}
		}
	}
	return true;
}

bool GA::calc(bool enableDisplay, bool enableOneLine)
{
	calcResult();
	for (int i = 0; i < max_genom_list; i++)//評価関数が最小の奴と最大のやつを検索
	{
		if (data[i].result < data[minNum].result)
		{
			minNum = i;
		}
		if (data[i].result > data[maxNum].result)
		{
			maxNum = i;
		}
	}
	if (data[maxNum].functionValue - eliteData.functionValue > 0.0001)
	{
		isChanged = true;
	}
	//評価関数が最もいいやつを保存
	data[minNum] = eliteData;

	calcResult();

	if (enableDisplay)
		displayValues(enableOneLine);

	return true;
}

bool GA::calcResult(bool enableSort)
{
	int maxNum = 0;

	for (int i = 0; i < max_genom_list; i++)
	{
		data[i].functionValue = std::sin(data[i].x[0] + data[i].x[1]) + std::pow((data[i].x[0] - data[i].x[1]), 2.0) - 1.5*data[i].x[0] + 2.5*data[i].x[1] + 1;//与えられた関数

		 //min:(-0.54719,-1.54719)=-1.9133//min:(-0.54719,-1.54719)=-1.9133
		if (data[maxNum].functionValue < data[i].functionValue)//座標の中で最も関数が大きいやつを検索
		{
			maxNum = i;
		}
	}
	double seg = data[maxNum].functionValue;//評価関数の切片を与えられた関数が最も大きいやつにセット
	double seg2 = searchRank(var_num - 2).functionValue - seg;
	resultSumValue = 0;
	double coefficient = 0.001 / var_num;//評価関数用の定数

	for (int i = 0; i < max_genom_list; i++)
	{
		bool flag = true;
		
		for (int j = 0; j < var_num; j++)
		{
			if (data[i].x[j] > varMax[j] || data[i].x[j] < varMin[j])//座標が場外にいるやつの処理
				flag = false;
		}
		data[i].result = std::pow((data[i].functionValue - seg)/seg2/coefficient, 2.0);//与えられた関数の値から切片で設定した値を引いて2乗する→与えられた関数の値が小さいやつが強くなる
		//data[i].result = std::abs(data[i].functionValue - seg);

		if (!flag)//場外に出たやつの処理
			data[i].result *= coefficient;
		resultSumValue += data[i].result;
	}
	if (enableSort)
		std::sort(data.begin(), data.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });

	return true;
}

int GA::random(int min, int max)
{
	//乱数の設定
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(engine);
}

double GA::random(int min, double max)
{
	return random((double)min, max);
}

double GA::random(double min, int max)
{
	return random(min, (double)max);
}

double GA::random(double min, double max)
{
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(engine);
}

bool GA::displayValues(bool enableOneLine)
{
	if (enableOneLine)
	{
		for (int j = 0; j < var_num; j++)
		{
			printf_s("%10.7lf,", eliteData.x[j]);//デバッグ用
		}
		printf_s(" \t f(x,y)=%10.7lf\t Result=%10.7lf\n", eliteData.functionValue, eliteData.result);

	}
	else
	{
		std::vector<Data> data_temp = data;
		std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });

		for (int i = 0; i < max_genom_list; i++)
		{
			for (int j = 0; j < var_num; j++)
			{
				printf_s("%10.7lf,", data_temp[i].x[j]);//デバッグ用
			}
			printf_s(" \t f(x,y)=%10.7lf\t Result=%10.7lf\n", data_temp[i].functionValue, data_temp[i].result);
		}
	}
	return true;
}

GA::Data GA::searchRank(int num)
{
	std::vector<Data> data_temp = data;
	std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue < y.functionValue; });
	return data_temp[num];
}



GA::~GA()
{

}
