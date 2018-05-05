#pragma once

#include<vector>
#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>

class GA
{
private:
	int max_genom_list;//個体数
	int var_num;//品物の個数
	double individualMutationRate = 0.05;//個体突然変異率
	double genomMutationRate = 0.1;
	int minNum = 0, maxNum = 0;
	double alpha = 0.3;
	std::vector<double> varMax, varMin;//変数の最小値・最大値
	bool isChanged = false;
public:
	double resultSumValue;//評価関数の合計

	class Data//データ格納用クラス
	{
	private:
		int var_num;//変数の数
	public:
		std::vector<double> x;//座標
		double functionValue;//与えられた関数の値
		double result;

		Data(int _var_num)//コンストラクタ
		{
			var_num = _var_num;

			x.resize(var_num);//isIncludedの配列の長さの設定
		}
	};

	std::vector<Data> data, prev_data;//操作前後で値を保持するために2個
	Data eliteData;
	GA(int _max_genom_list, int _var_num, std::vector<double> _varMax, std::vector<double> _varMin);	//コンストラクタ
	bool selection();//選択

	bool blxAlphaCrossover();
	bool mutation();//突然変異
	bool calc(bool enableDisplay, bool enableOneLine = false);//評価関数の計算
private:
	bool calcResult(bool enableSort = false);
	int random(int min, int max);
	double random(int min, double max);
	double random(double min, int max);
	double random(double min, double max);
	bool displayValues(bool enableOneLine);
	Data searchRank(int num);
public:
	~GA();//デコンストラクタ
};


