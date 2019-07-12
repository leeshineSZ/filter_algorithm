#include "DataProcess.h"


#ifdef FAST_SQRT_CAL_ENABLE
#include <stdint.h>
static int max_int_of_sqrt_within100(uint8_t num)
{
	if (num > 100)
	{
		return -1;
	}
	if (0 == num)
	{
		return 0;
	}
	int i = 10;
	while (i*i > num)
	{
		i--;
	}
	return i;
}
static int fastsqrt(uint32_t num) //最大不超过32位的数
{
	if (num <= 100)
	{
		return max_int_of_sqrt_within100(num);
	}
	uint8_t tem[6] = { 0 }; //2^32 = 42 94 96 72 96
	uint32_t temp1 = num;
	int i = 0;
	while (temp1)
	{
		tem[i] = (uint8_t)(temp1 % 100);
		temp1 = temp1 / 100;
		i++;
	}
	uint32_t ret_val = 0;
	uint32_t remd; //余数
	uint32_t divs;//除数
	ret_val = max_int_of_sqrt_within100(tem[i - 1]);
	remd = tem[i - 1] - ret_val*ret_val;
	i = i - 1;
	while (i > 0)
	{
		remd = remd * 100 + tem[i - 1];
		divs = ret_val * 20;
		uint32_t quots_tem = divs;
		uint32_t tem_quot = remd / divs;
		while (1)
		{
			divs = quots_tem + tem_quot;
			if (divs* tem_quot > remd)
			{
				tem_quot--;
			}
			else
			{
				break;
			}
		}
		ret_val = ret_val * 10 + tem_quot;
		remd = remd - tem_quot* divs;
		i = i - 1;
	}
	return ret_val;
}
float GetFastSqrtf(float data)
{
	int32_t u1 = (int32_t)(data);
	if (data < 0.f)
	{
		return -1;
	}
	else if (u1 < 0 )
	{
		uint32_t m_data = (uint32_t)(data / 4);
		return  (2 * fastsqrt(m_data));
	}
	else if (data > 400000)
	{
		return fastsqrt((uint32_t)data);
	}
	else if (data > 10000)
	{
		uint32_t m_data = (uint32_t)(data * 10000);
		return fastsqrt(m_data)/100.f;
	}
	else // (data < 10000)
	{
		uint32_t m_data = (uint32_t)(data * 360000);
		return fastsqrt(m_data) / 600.f;
	}
}
#endif//FAST_SQRT_CAL_ENABLE

int  average_value_filter_init(AVG_VL_PARAM_T* parg, int len)
{
	if (len > 1)
	{
		parg->data_count = 0;
		parg->agv_v = 0;
		parg->avg_len = len;
		return 0;
	}
	else
	{
		parg->data_count = 0;
		parg->agv_v = 0;
		parg->avg_len = 0;
		return -1;
	}
}
float average_value_filter(AVG_VL_PARAM_T* parg, float v)
{
	parg->data_count++;
	if (parg->data_count < parg->avg_len)
	{
		parg->agv_v = parg->agv_v + (v - parg->agv_v) / parg->data_count;
	}
	else
	{
		//parg->agv_v = (parg->agv_v * (parg->avg_len -1 ) + v)/ parg->avg_len;
		parg->agv_v = parg->agv_v + (v - parg->agv_v) / parg->avg_len;
	}
	return parg->agv_v;
}

int	  ButterScndOrderLP_HP_CoefCal(SCND_FILTER_COEF_T* parg, int  itype, float fs, float f_cut)
{
	if (FILTER_HIGH_PASS_TYPE != itype && FILTER_LOW_PASS_TYPE != itype)
	{
		return -1;
	}
	if (fs < 2 * f_cut)
	{
		return -1;
	}
	memset(parg, 0, sizeof(SCND_FILTER_COEF_T));
	parg->order = 2;
	float Q = 0.707f;
	float dBgain = 0;

	float pi = 2 * asinf(1);
	float A = 0;
	float w0 = 2 * pi* f_cut / fs;
	float cosw0 = cosf(w0);
	float sinw0 = sinf(w0);
	float alpha = sinw0 / (2 * Q);
	float gain = powf(10, dBgain / 20.0f);
	float b[3] = { 0 };
	float a[3] = { 0 };

	if (FILTER_LOW_PASS_TYPE == itype)
	{
		parg->Type = FILTER_HIGH_PASS_TYPE;
		//low pass filter coefficients

		b[0] = gain*((1 - cosw0) / 2);
		b[1] = gain*(1 - cosw0);
		b[2] = gain*((1 - cosw0) / 2);
		a[0] = 1 + alpha;
		a[1] = -2 * cosw0;
		a[2] = 1 - alpha;
	}
	else if (FILTER_HIGH_PASS_TYPE == itype)
	{
		//high pass filter coefficients
		b[0] = gain*((1 + cosw0) / 2);
		b[1] = gain*(-(1 + cosw0));
		b[2] = gain*((1 + cosw0) / 2);
		a[0] = 1 + alpha;
		a[1] = -2 * cosw0;
		a[2] = 1 - alpha;
	}
	//Normalize so that A0 = 1;
	parg->coef[0] = b[0] / a[0];
	parg->coef[1] = b[1] / a[0];
	parg->coef[2] = b[2] / a[0];
	parg->coef[3] = a[1] / a[0];
	parg->coef[4] = a[2] / a[0];
	return 0;
}

float ButterScndOrderFilterProcess(SCND_FILTER_COEF_T* parg, float samp)
{
	if (NULL == parg)
	{
		return 0;
	}
	float Yn = 0;
	{
		float *pcoef = parg->coef;
		int  i;
		Yn = pcoef[0] * samp;
		for (i = 1; i <= parg->order; i++)
		{
			Yn += pcoef[i] * parg->x[i - 1];
			Yn -= pcoef[i + parg->order] * parg->y[i - 1];
		}
		for (i = (parg->order - 1); i > 0; i--)
		{
			parg->x[i] = parg->x[i - 1];
			parg->y[i] = parg->y[i - 1];
		}
		parg->x[0] = samp;
		parg->y[0] = Yn;
	}
	return Yn;
}

float BandPassOrEliminateFilter(BPFILTER_DATA_PARAM_T* parg, float samp)//只针对二阶带通带阻滤波器
{
	float Yn = 0;
	const float *pcoef = parg->coefs;
	int  i;
	Yn = pcoef[0] * samp;
	for (i = 1; i <= 4; i++)
	{
		Yn += pcoef[i] * parg->Data_X[i - 1];
		Yn -= pcoef[i + 4] * parg->Data_Y[i - 1];
	}
	for (i = (4 - 1); i > 0; i--)
	{
		parg->Data_X[i] = parg->Data_X[i - 1];
		parg->Data_Y[i] = parg->Data_Y[i - 1];
	}
	parg->Data_X[0] = samp;
	parg->Data_Y[0] = Yn;
	return Yn;
}
float KalmanFilter(float data, float *pXk_1, float *pPk_1, float Q, float R)
{
	float Xk, Pk, kg;
	Xk = *pXk_1;
	Pk = *pPk_1 + Q;
	kg = Pk / (Pk + R);
	*pXk_1 = Xk + kg * (data - Xk);
	*pPk_1 = (1.0f - kg) * Pk;
	return *pXk_1;
}
float HighPassFilter(const float *pb, float *x, float *y, float samp)//一阶高通滤波
{
	const float *pcoff = pb;
	x[0] = samp;
	y[0] =  pcoff[0]*x[0] + pcoff[1]*x[1]  - pcoff[2] * y[1];

	y[1] = y[0];
	x[1] = x[0];
	return y[0];
}
float ButterFilter(const float *pb, float *x, float *y, float samp)//巴特沃斯二阶滤波
{
	const float *pcoff = pb;
	float ftem = 0;
	/* x[0]为上一个采样值，x[1]为上上个采样值，y[0]为上次滤波后的值，y[1]为上上次滤波后的值*/
	ftem = pcoff[0] * samp + pcoff[1] * x[0] + pcoff[2] * x[1] - (pcoff[3] * y[0] + pcoff[4] * y[1]);
	x[1] = x[0];
	x[0] = samp;

	y[1] = y[0];
	y[0] = ftem;
	return ftem;
}
float ChebyFilter(const float *pb, float *x, float *y, float samp) //切比雪夫三阶滤波器
{
	const float *pcoff = pb;
	float ftem = 0;
	/* x[0]为上一个采样值，x[1]为上上个采样值，x[2]为上上上次采样值，y[0]为上次滤波后的值，y[1]为上上次滤波后的值*/
	ftem = pcoff[0] * samp + pcoff[1] * x[0] + pcoff[2] * x[1] + pcoff[3] * x[2]
		- (pcoff[4] * y[0] + pcoff[5] * y[1] + pcoff[6] * y[2]);
	x[2] = x[1];
	x[1] = x[0];
	x[0] = samp;
	y[2] = y[1];
	y[1] = y[0];
	y[0] = ftem;
	return ftem;
}
float LowPassFilter(float data, float ratio, float* Yn_1)//低通滤波
{
	*Yn_1 = ratio * data + (1 - ratio) * (*Yn_1);
	return *Yn_1;
}
float HighOrLowPassFilter(FILTER_DATA_PARAM_T* parg, float samp)//针对二阶或者三四阶的高低通初始化并滤波
{
	float Yn = 0;
	if (parg->DataInitCnt == 0)
	{
		const float *pcoef = parg->coefs;
		int  i;
		Yn = pcoef[0] * samp;
		for (i = 1; i <= parg->FilterOrder; i++)
		{
			Yn += pcoef[i] * parg->Data_X[i - 1];
			Yn -= pcoef[i + parg->FilterOrder] * parg->Data_Y[i - 1];
		}
		for (i = (parg->FilterOrder - 1); i > 0; i--)
		{
			parg->Data_X[i] = parg->Data_X[i - 1];
			parg->Data_Y[i] = parg->Data_Y[i - 1];
		}
		parg->Data_X[0] = samp;
		parg->Data_Y[0] = Yn;
	}
	else
	{
		parg->DataInitCnt--;
		if (FILTER_LOW_PASS_TYPE == parg->FilterType)
		{
			parg->Data_X[parg->DataInitCnt] = samp;
			Yn = parg->Data_Y[parg->DataInitCnt] = samp;
			if (0 == parg->DataInitCnt)
			{
				float mean_val = 0;
				for (int i = 0; i < parg->FilterOrder;i++)
				{
					mean_val += parg->Data_Y[i];
				}
				mean_val = mean_val / parg->FilterOrder;
				for (int i = 0; i < parg->FilterOrder; i++)
				{
					parg->Data_Y[i] = mean_val;
				}
			}
		}
		else if (FILTER_HIGH_PASS_TYPE == parg->FilterType)
		{
			parg->Data_X[parg->DataInitCnt] = parg->Data_Y[parg->DataInitCnt] =  samp;
			Yn = 0;
			if (0 == parg->DataInitCnt)
			{
				float mean_val = 0;
				for (int i = 0; i < parg->FilterOrder; i++)
				{
					mean_val += parg->Data_Y[i];
				}
				mean_val = mean_val / parg->FilterOrder;
				for (int i = 0; i < parg->FilterOrder; i++)
				{
					parg->Data_Y[i] -= mean_val;
				}
			}
		}
	}
	return Yn;
}

#if 0
int Interpol(const double* x, const double *fx, int n, double *z, double *pz, int m)
{
	double term, *table, *coeff;
	int i, j, k;
	if ((table = (double*)malloc(sizeof(double)*n)) == NULL)
		return -1;
	if ((coeff = (double*)malloc(sizeof(double)*n)) == NULL)
	{
		free(table);
		return -1;
	}
	/*Initialize the coeffients*/
	memcpy(table, fx, sizeof(double)*n);
	/*Determine the coefficient of the interpolating polynimial*/
	coeff[0] = table[0];
	for (k = 1; k < n; k++)
	{
		for (i = 0; i < n - k; i++)
		{
			j = i + k;
			table[i] = (table[i + 1] - table[i]) / (x[j] - x[i]);
		}
		coeff[k] = table[0];
	}
	free(table);
	/*Evaluate the interpolating polynomial at the specified points */
	for (k = 0; k < m; k++)
	{
		pz[k] = coeff[0];
		for (j = 1; j < n; j++)
		{
			term = coeff[j];
			for (i = 0; i < j; i++)
			{
				term = term*(z[k] - x[i]);
			}
			pz[k] = pz[k] + term;
		}	 
	}
	free(coeff);
	return 0;
}
#endif
int FitCurveAndGetMax(const float x[4], const float fx[4], float *z, float *pz)//用fx = ax3+bx2+cx+d对数据点进行拟合，然后求最大值及最大值对应的位置
{
	float  coeff[4] = {0};
	coeff[0] = (fx[3] - 3 * fx[2] + 3 * fx[1] - fx[0]) / 6;
	coeff[1] = (fx[3] -2 * fx[2] + fx[1] - coeff[0]*6* x[2]) / 2;
	//coeff[2] = (fx[3] - fx[2]) - coeff[0]*(3*powf(x[2],2)+3* x[2]+1)- coeff[1]*(2* x[2]+1);
	coeff[2] = (fx[3] - fx[2]) - coeff[0] * (3 * x[2]* x[2] + 3 * x[2] + 1) - coeff[1] * (2 * x[2] + 1);
	//coeff[3] = fx[2]- coeff[0]* powf(x[2], 3) - coeff[1]* powf(x[2], 2) - coeff[2]* x[2];
	coeff[3] = fx[2] - coeff[0] * x[2]* x[2]* x[2] - coeff[1] * x[2]* x[2] - coeff[2] * x[2];

	double a, b, c, deta, x1, x2;
	a = 3 * coeff[0];
	b = 2 * coeff[1];
	c = coeff[2];
	deta = sqrtf(b*b - 4 *a*c);
	x1 = (-b + deta) / (2 * a);
	x2 = (-b - deta) / (2 * a);
	if (0 == coeff[0])
	{
		*z = x[2];
		*pz = fx[2];
		return 0;
	}
	
	if (x1 >= x[0] && x1 <= x[3])
	{
		*z = (float)x1;
	}
	if (x2 >= x[0] && x2 <= x[3])
	{
		*z = (float)x2;
	}
	//*pz = coeff[0] * powf(*z, 3) + coeff[1] * powf(*z, 2) + coeff[2] * (*z) + coeff[3];
	*pz = coeff[0] * (*z)*(*z)*(*z) + coeff[1] * (*z)*(*z) + coeff[2] * (*z) + coeff[3];
	return 0;
}
int FitCurveWithCubicCurve(const float x[4], const float fx[4], double coef[4])//用fx = ax3+bx2+cx+d对曲线进行拟合,求a,b,c,d系数，shl 2017-04-12
{
	double CK01, CK02, CK03;
	double CK11, CK12, CK13;
	double CK21, CK22, CK23;
	double CK31, CK32, CK33;
	//CK01 = powf(x[0], 3);
	CK01 = x[0]* x[0]* x[0];
	//CK02 = powf(x[0], 2);
	CK02 = x[0]* x[0];
	CK03 = x[0];
	//CK11 = powf(x[1], 3);
	//CK12 = powf(x[1], 2);
	CK11 = (x[1]* x[1]* x[1]);
	CK12 = (x[1]* x[1]);
	CK13 = x[1];
	//CK21 = powf(x[2], 3);
	//CK22 = powf(x[2], 2);
	CK21 = (x[2]* x[2]* x[2]);
	CK22 = (x[2]* x[2]);
	CK23 = x[2];
	//CK31 = powf(x[3], 3);
	//CK32 = powf(x[3], 2);
	CK31 = (x[3]* x[3]* x[3]);
	CK32 = (x[3]* x[3]);
	CK33 = x[3];
	double G1, G2, G3;
	G1 = CK13 - CK03;
	G2 = CK23 - CK13;
	G3 = CK33 - CK23;
	double P1, P2, P3, P4;
	P1 = G1*(CK21 - CK11) - G2*(CK11 - CK01);
	P2 = G1*(CK22 - CK12) - G2*(CK12 - CK02);
	P3 = G1*(CK31 - CK21) - G3*(CK11 - CK01);
	P4 = G1*(CK32 - CK22) - G3*(CK12 - CK02);

	coef[0] = (P2*(G1*(fx[3] - fx[2]) - G3*(fx[1] - fx[0])) - P4*(G1*(fx[2] - fx[1]) - G2*(fx[1] - fx[0]))) / (P3*P2 - P1*P4);
	coef[1] = ((G1*(fx[2] - fx[1]) - G2*(fx[1] - fx[0])) - P1*coef[0]) / P2;
	coef[2] = ((fx[1] - fx[0]) - coef[0] * (CK11 - CK01) - coef[1] * (CK12 - CK02)) / G1;
	coef[3] = fx[0] - coef[0] * CK01 - coef[1] * CK02 - coef[2] * CK03;
	return 1;
}
int FitCurveWithQuadraticCurve(const float x[3], const float fx[3], double coef[3])//用fx = ax2+bx+c对曲线进行拟合,求a,b,c系数，shl 2017-04-12
{
	float G1 = x[1] - x[0];
	float G2 = x[2] - x[0];
	//double G3 = powf(x[2], 2) - powf(x[0], 2);
	//double G4 = powf(x[1], 2) - powf(x[0], 2);
	double G3 = (x[2]* x[2]) -(x[0]* x[0]);
	double G4 = (x[1]*x[1]) - (x[0]* x[0]);
	coef[0] = (G1*(fx[2] - fx[0]) - G2*(fx[1] - fx[0])) / (G1*G3 - G2*G4);
	coef[1] = ((fx[1] - fx[0]) - coef[0] * G4) / G1;
	//coef[2] = fx[0] - coef[0] * powf(x[0], 2) - coef[1] * x[0];
	coef[2] = fx[0] - coef[0] * (x[0]* x[0]) - coef[1] * x[0];
	return 1;
}
int FitCurveAndGetMin(const float x[3], const float fx[3], float *z, float *pz)
{
	float coeff[3] = {0};
	coeff[0] = (fx[2] - 2 * fx[1] + fx[0]) / 2;
	coeff[1] = (fx[2] - fx[0] - 2 * x[1] * (fx[2] - 2 * fx[1] + fx[0])) / 2;
	coeff[2] = (float)(2 * fx[1] - (fx[2] - fx[0])* x[1] + (fx[2] - 2 * fx[1] + fx[0])*(x[1]* x[1])) / 2;
	*z = -coeff[1] / (2 * coeff[0]);
	*pz = (4 * coeff[0] * coeff[2] - (coeff[1]* coeff[1])) / (4 * coeff[0]);
	return 0;
}
int GetMaxPntWithCrsPntOfStrLine(const float x[5], const float fx[5], float *z, float *pz)
{
	if (fx[1] > fx[3])
	{
		*z = (fx[2]+ fx[2]* x[2] - fx[3]* x[2]+ fx[0]+ fx[1] * x[2] - 2* fx[1] - fx[0] * x[2]) / (fx[1] - fx[0] - fx[3] + fx[2]);
		*pz = (fx[1] - fx[0])*(*z) + fx[0] -(fx[1] - fx[0])*(x[2] - 2);
	}
	if (fx[1] < fx[3])
	{
		*z = (fx[3]- fx[2]+(fx[2] - fx[1])*x[2] - (fx[4] - fx[3])*(x[2]+1)) / (fx[2] - fx[1] - fx[4] + fx[3]);
		*pz = (fx[2] - fx[1])*(*z) +  fx[2]-(fx[2] - fx[1])*x[2];
	}
	return 0;
}
int GetIntersecOfStr_LineWithX_Axis(const float x[2], const float fx[2], float *z, float *pz)//获取直线与X轴的交点
{
	if (x[1] == x[0])
	{
		*z = x[1];
		*pz = 0;
		return  0;
	}
	if (fx[1] == fx[0])
	{
		return -1;
	}
	float a, b;
	a = (fx[1] - fx[0])/(x[1]-x[0]);
	b = fx[0] - a*x[0];
	*z = -b / a;
	*pz = 0;
	return 0;
}
void Lsqe(const float *x, const float *y, int n, float *b1, float* b0)
{
	float sumx, sumy, sumx2, sumxy;
	int i;
	/*Compute the required summations*/
	sumx = 0;
	sumy = 0;
	sumx2 = 0;
	sumxy = 0;
	for (i = 0; i < n; i++)
	{
		sumx = sumx + x[i];
		sumy = sumy + y[i];
		//sumx2 = sumx2 + powf(x[i], 2);
		sumx2 = sumx2 + (x[i]*x[i]);
		sumxy = sumxy + (x[i] * y[i]);
	}
	/*Compute the least_squares estimators*/
	//*b1 = (sumxy - ((sumx*sumy) / (float)n)) / (sumx2 - (powf(sumx, 2) / (float)n));
	*b1 = (sumxy - ((sumx*sumy) / (float)n)) / (sumx2 - ((sumx*sumx) / (float)n));
	*b0 = (sumy - ((*b1)*sumx)) / (float)n;
}

/*最小二乘估计，y= b1/x+b0,估计b1 和b0 的值去拟合,这样，y=b1/x+b0是由参数x和y所指定的点集的最佳拟合线*/
void LsqeAntyRateCurve(const float *x, const float *y, int n, float *b1, float* b0)
{
	float sumx, sumy, sumx2, sumxy;
	/*Compute the required summations*/
	sumx = 0;
	sumy = 0;
	sumx2 = 0;
	sumxy = 0;
	int i;
	for (i = 0; i < n; i++)
	{
		sumx += 1.0f / x[i];
		sumy += y[i];
		//sumx2 += 1.0f / powf(x[i], 2);
		sumx2 += 1.0f / (x[i]*x[i]);
		sumxy += (float)y[i] / x[i];
	}
	//*b1 = (sumxy - (sumy* sumx / (float)n)) / (sumx2 - (pow(sumx, 2) / (float)n));
	*b1 = (sumxy - (sumy* sumx / (float)n)) / (sumx2 - ((sumx*sumx) / (float)n));
	*b0 = (sumy - ((*b1)*sumx)) / (float)n;
}
int	  MeanFilterInit(MEAN_VAL_FILTER_PARAM_T* parg, float* data_cache, int cache_size)//均值滤波
{
	if (NULL == data_cache)
	{
		return -1;
	}	
	memset(parg, 0, sizeof(MEAN_VAL_FILTER_PARAM_T));
	parg->DataBuf = data_cache;
	parg->CacheSize = cache_size;
	return 0;
}
float MeanFilter(MEAN_VAL_FILTER_PARAM_T* parg, float samp) //
{
	if (NULL == parg->DataBuf)
	{
		return 0;
	}
	parg->DataSum += samp;
	if (parg->DataCnt < parg->CacheSize)
	{
		parg->DataBuf[parg->DataCnt++] = samp;
		return parg->DataSum / parg->DataCnt;
	}
	else
	{
		parg->DataSum -= parg->DataBuf[parg->DataIndex];
		parg->DataBuf[parg->DataIndex++] = samp;
		if (parg->DataIndex == parg->CacheSize)
		{
			parg->DataIndex = 0;
		}
		return (float)(parg->DataSum / parg->CacheSize);
	}
}
int   MedianFilterInit(MEAN_VAL_FILTER_PARAM_T* parg, float* data_cache, int cache_size)//中值滤波
{
	if (NULL == data_cache)
	{
		return -1;
	}	
	if (cache_size < 3)
	{
		return -1;
	}
	memset(parg, 0, sizeof(MEAN_VAL_FILTER_PARAM_T));
	parg->DataBuf = data_cache;
	parg->CacheSize = cache_size;
	return 0;
}
float MedianFilterPutValue(MEAN_VAL_FILTER_PARAM_T* parg, float samp)//均值滤波的数据个数小于len时，返回samp
{	
	if (NULL == parg->DataBuf)
	{
		return 0;
	}
	float ret = 0;
	parg->DataSum += samp;
	if (parg->DataCnt < parg->CacheSize)
	{
		parg->DataBuf[parg->DataCnt++] = samp;	
	}
	else
	{
		parg->DataSum -= parg->DataBuf[parg->DataIndex];
		parg->DataBuf[parg->DataIndex++] = samp;
		if (parg->DataIndex == parg->CacheSize)
		{
			parg->DataIndex = 0;
		}
	}
	if (parg->DataCnt >= 3)
	{
		float temMin = parg->DataBuf[0];
		float temMax = parg->DataBuf[0];
		for (int i = 1; i < parg->DataCnt; i++)
		{
			if (parg->DataBuf[i] > temMax)
			{
				temMax = parg->DataBuf[i];
			}
			if (parg->DataBuf[i] < temMin)
			{
				temMin = parg->DataBuf[i];
			}
		}
		ret = (parg->DataSum - temMax - temMin) / (parg->DataCnt - 2);
	}
	return ret;
}
int  MeanSquareParamInit(MEAN_SQUARE_PARAM_T *parg, float* data_cache, int cache_size) //求均方差初始化
{
	if (NULL == data_cache)
	{
		return -1;
	}
	if (cache_size < 3)
	{
		return -1;
	}
	memset(parg, 0, sizeof(MEAN_SQUARE_PARAM_T));
	parg->DataCache = data_cache;
	parg->CacheSize = cache_size;
	return 0;
}
float MeanSquareParamPutValue(MEAN_SQUARE_PARAM_T *parg, float ps) //求均方差
{
	if (NULL == parg->DataCache)
	{
		return 0;
	}
	parg->DataSum += ps;
	if (parg->DataCount < parg->CacheSize)
	{
		parg->DataCache[parg->DataCount++] = ps;
	}
	else
	{
		parg->DataSum -= parg->DataCache[parg->DataIndex];
		parg->DataCache[parg->DataIndex++] = ps;
		if (parg->DataIndex == parg->CacheSize)
		{
			parg->DataIndex = 0;
		}
	}
	parg->DataMSV = 0;
	parg->DataMeanValue = parg->DataSum / parg->DataCount;

	for (int i = 0; i < parg->DataCount; i++)
	{
		float diff_val = parg->DataCache[i] - parg->DataMeanValue;
		//parg->DataMSV += powf(parg->DataCache[i] - parg->DataMeanValue, 2);
		parg->DataMSV += diff_val*diff_val;
	}
#ifdef FAST_SQRT_CAL_ENABLE
	parg->DataMSV = GetFastSqrtf(parg->DataMSV / parg->DataCount);
#else
	parg->DataMSV = sqrtf(parg->DataMSV / parg->DataCount);
#endif	
	return parg->DataMSV;
}
int RCLowPassFilterInit(RC_FILTER_PARAM_T* parg, float fs, float cutoff_fs, float rs)
{
	memset(parg, 0, sizeof(RC_FILTER_PARAM_T));
	float pi = 2 * asinf(1);
	parg->a_default_coef = cutoff_fs * 2 * pi * 1.0f / fs;
	parg->RS_Threshold = rs;
	parg->Yn_1 = 0;
	parg->unequal_cnt = 0;
	return 0;
}
float RCLowPassFilter(RC_FILTER_PARAM_T* parg, float data)
{
	float diff_vl = data - parg->Yn_1;
	//float Rsq = powf((data - parg->Yn_1), 2);
	float Rsq = diff_vl*diff_vl;
	parg->a_curr_coef = Rsq < parg->RS_Threshold ? parg->a_default_coef:(parg->RS_Threshold / Rsq*parg->a_default_coef);
	if (parg->a_curr_coef == parg->a_default_coef)
	{
		parg->unequal_cnt = 0;
	}
	else
	{
		parg->unequal_cnt++;
	}
	if (parg->unequal_cnt > 5)
	{
		parg->Yn_1 = data;
		parg->a_curr_coef = parg->a_default_coef;
		parg->unequal_cnt = 0;
	}
	float Yn = parg->a_curr_coef* data + (1 - parg->a_curr_coef)* parg->Yn_1;
	parg->Yn_1 = Yn;

	return Yn;
}

int  FirFilterInit(FIR_FILTER_PARAM_T* parg, const float* coef, float* Xn, int Order) //系数关于中心对称，
{
	if (NULL == coef || NULL == Xn )
	{
		return -1;
	}
	parg->coefs = coef;
	parg->Xn = Xn;
	parg->Fir_Order = Order;
	return 0;
}
float FirFilterProcess(FIR_FILTER_PARAM_T* parg, float m_data)
{
#if 0
	int i = parg->Fir_Order;
	double sum = 0;
	while (i)
	{
		sum += parg->coefs[i] * parg->Xn[i - 1];
		if (i > 1)
		{
			parg->Xn[i - 1] = parg->Xn[i - 2];
		}
		i--;
	}
	sum += parg->coefs[0] * m_data;
	parg->Xn[0] = m_data;
	return (float)sum;

#endif
	int i = parg->Fir_Order; 
	int j = parg->Fir_Order / 2; //
	int ord = parg->Fir_Order;
	double sum = 0;
	const float *pcoef = parg->coefs;
	while (i)
	{
		if (i > j )
		{
			sum += pcoef[ord - i ] * parg->Xn[i - 1];
		}
		else
		{
			sum += pcoef[i] * parg->Xn[i - 1];
		}		
		if (i > 1)
		{
			parg->Xn[i - 1] = parg->Xn[i - 2];
		}
		i--;
	}
	sum += pcoef[0] * m_data;
	parg->Xn[0] = m_data;
	return (float)sum;
}
int fst_order_highpass_filter_param_init(FST_ODR_HP_FIL_PARAM_T* parg, float radio)
{
	if (radio < 0 || radio >1)
	{
		return -1;
	}
	memset(parg, 0, sizeof(FST_ODR_HP_FIL_PARAM_T));
	parg->radio = radio;
	return 0;
}
float fst_order_highpass_filter(FST_ODR_HP_FIL_PARAM_T* parg, float samp)
{
	float Yn = 0;
	if (0 == parg->is_use)
	{
		parg->is_use = 1;
		parg->pre_xn = samp;
	}
	Yn = parg->pre_yn * parg->radio + (samp - parg->pre_xn);
	parg->pre_xn = samp;
	parg->pre_yn = Yn;
	return Yn;
}

void BubbleSort(float a[], int n)// 冒泡排序
{
	int i, j;
	float temp;
	for (j = 0; j < n - 1; j++)
		for (i = 0; i < n - 1 - j; i++)
		{
			if (a[i] > a[i + 1])
			{
				temp = a[i];
				a[i] = a[i + 1];
				a[i + 1] = temp;
			}
		}
}
void ClearValueOfFtpt(FEATURE_POINT* pt1)//特征点的清零操作
{
	pt1->offSet = 0;
	pt1->value = 0;
}

void GetFitedMinValueInWave(int dots,float buf[], float *z, float *pz)
{
	float x1[3] = { 0 };
	float fx1[3] = { 0 };
	for (int i = dots - 1; i < dots + 2; i++)
	{
		int j = i - (dots - 1);
		x1[j] = i;
		fx1[j] = buf[i];
	}
	FitCurveAndGetMin(x1, fx1, z, pz);
}
void GetFitedMaxValueInWave(int dots, float buf[], float *z, float *pz)
{
	float x1[4] = { 0 };
	float fx1[4] = { 0 };
	for (int i = dots - 2; i < dots + 2; i++)
	{
		int j = i - (dots - 2);
		x1[j] = i;
		fx1[j] = buf[i];
	}
	FitCurveAndGetMax(x1, fx1, z, pz);
}
void GetZeroPointOfStrLineWithX_Axis(FEATURE_POINT* fp, float* z, float *pz)
{
	float x[2] = { 0 };
	float fx[2] = { 0 };
	for (int i = 0; i < 2; i++)
	{
		x[i] = fp[i].offSet;
		fx[i] = fp[i].value;
	}
	GetIntersecOfStr_LineWithX_Axis(x, fx, z, pz);
}
void GetMaxPntOfTwoStrLineCross(int dots, float buf[], float* z, float *pz)
{
	float x[5] = { 0 };
	float fx[5] = { 0 };
	for (int i = dots - 2; i <= dots + 2; i++)
	{
		int j = i - (dots - 2);
		x[j] = i;
		fx[j] = buf[i];
	}
	GetMaxPntWithCrsPntOfStrLine(x, fx, z, pz);
}
