#ifndef __DATA_PROCESS_H__
#define __DATA_PROCESS_H__       
/**********************************************************************************
***@func: 定义相应参数结构,********************************************************
***@Warning: FEATURE_POINT中的地址偏移最好以基点所在的点的横坐标为起点的偏移*******
**********************************************************************************/
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif  

#define FILTER_HIGH_PASS_TYPE     0X01
#define FILTER_LOW_PASS_TYPE	  0X10


typedef struct feature_point
{
	int				offSet;			 //偏移
	float			value;
}FEATURE_POINT;						//特征点
typedef struct mean_square_param_t //均方差
{
	float*	DataCache;
	float	DataSum; 
	float	DataMeanValue;
	float	DataMSV; //  DataMeanSquareValue
	int		DataCount;
	int		DataIndex;
	int		CacheSize;
}MEAN_SQUARE_PARAM_T;
typedef struct mean_val_filter_param_t //均值滤波
{
	double		DataSum;	//数值总和
	int			DataCnt;	//数据个数计数
	float*		DataBuf;
	int			DataIndex;
	int			CacheSize;
}MEAN_VAL_FILTER_PARAM_T;

typedef struct klm_filter_param_t
{
	float Xk_1;
	float Pk_1;

	float Q;
	float Q0;
	float R;
}KLM_FILTER_PARAM_T;
typedef struct filter_data_param_t
{
	float Data_X[4];
	float Data_Y[4];
	const float* coefs;
	short  FilterType;	//滤波器类型
	short  DataInitCnt; //需要初始化的数据个数
	short  FilterOrder; //滤波器阶数
}FILTER_DATA_PARAM_T;

typedef struct bpfilter_data_param_t //二阶带通或者带阻滤波器参数
{
	float Data_X[4];
	float Data_Y[4];
	const float* coefs;
}BPFILTER_DATA_PARAM_T;

typedef struct rc_filter_param_t
{
	float	a_curr_coef;
	float	a_default_coef;
	int		unequal_cnt;
	float	Yn_1;
	float	RS_Threshold;
}RC_FILTER_PARAM_T;
typedef struct fir_filter_param_t
{
	const float *coefs;
	float		*Xn;
	int			Fir_Order;
}FIR_FILTER_PARAM_T;

typedef struct _average_value_filter_param_t
{
	int 	avg_len; // 求多少个数据的平均值
	int 	data_count; //已经输入了多少个数
	float   agv_v;		//平均值
}AVG_VL_PARAM_T;

typedef struct _first_order_hp_filter_param_t
{
	int is_use;
	float pre_xn;
	float pre_yn;
	float radio; // 高通系数，取值小于1，越靠近1，滤波效果越差；
}FST_ODR_HP_FIL_PARAM_T;

typedef struct scnd_filter_coef_t
{
	int Type;	//	Lp or Hp or Bp
	int order;	//都会初始化为2
	float x[2];
	float y[2];
	float coef[5];//因为是二阶系数，因此，A0 = 1，因此，只有五个系数
}SCND_FILTER_COEF_T;

#ifdef FAST_SQRT_CAL_ENABLE
//均方差快速算法，数据大于400000.则结果误差在1以内，小于等于400000，则结果误差在0.01以内；
float GetFastSqrtf(float data);

#endif

int  average_value_filter_init(AVG_VL_PARAM_T* parg, int len);
float average_value_filter(AVG_VL_PARAM_T* parg, float v);
int	  ButterScndOrderLP_HP_CoefCal(SCND_FILTER_COEF_T* parg, int  itype, float fs, float f_cut); //注意，巴特沃斯二阶滤波器，在截止频率处，下降3db																							
float ButterScndOrderFilterProcess(SCND_FILTER_COEF_T* parg, float data);

float BandPassOrEliminateFilter(BPFILTER_DATA_PARAM_T* parg, float samp); //二阶带通或者带阻滤波器
float KalmanFilter(float data, float *pXk_1, float *pPk_1, float Q, float R);
float HighPassFilter(const float *pb, float *x, float *y, float samp);//一阶高通滤波
float ButterFilter(const float *pb, float *u, float *y, float samp);//巴特沃斯二阶滤波器
float ChebyFilter(const float *pb, float *x, float *y, float samp); //切比雪夫三阶滤波器
float LowPassFilter(float data, float ratio, float* Yn_1);//低通滤波
float HighOrLowPassFilter(FILTER_DATA_PARAM_T* parg, float samp);//针对二阶或者三四阶的高低通初始化并滤波
/*多项式插值运算，其中，x为已知的横坐标，fx为已知纵坐标，n为已知横坐标和纵坐标的点数，
z为需要插值的数，pz为计算出来的插值的值，m为待插值的个数。
插值成功返回0，失败返回-1；*/
int Interpol(const double* x, const double *fx, int n, double *z, double *pz, int m);
/*用fx = ax3+bx2+cx+d对数据点进行拟合，然后求最大值及最大值对应的位置*****
**warning: 其中，fx[2]为已知最大值,x[2]为已知最大值对应的点数，x数组内部存放四个连续的点对应的偏移，成功返回0，失败返回-1；*/
int FitCurveAndGetMax(const float x[4], const float fx[4], float *z, float *pz);
int FitCurveWithCubicCurve(const float x[4], const float fx[4], double coef[4]);//用fx = ax3+bx2+cx+d对曲线进行拟合,求a,b,c,d系数，shl 2017-04-12
int FitCurveWithQuadraticCurve(const float x[3], const float fx[3], double coef[3]);//用fx = ax2+bx+c对曲线进行拟合,求a,b,c系数，shl 2017-04-12
/*用fx = ax2+bx+c对数据点进行拟合，然后求最小值以及最小值对应的位置，
其中，fx[1]为已知最小点，x[1]为已知最小值对应的点数，成功返回0，失败返回-1；*/
int FitCurveAndGetMin(const float x[3], const float fx[3], float *z, float *pz);
/*用两直线相交的方法求取交点，fx[2]为最值点*********/
int GetMaxPntWithCrsPntOfStrLine(const float x[5], const float fx[5], float *z, float *pz);
int GetIntersecOfStr_LineWithX_Axis(const float x[2], const float fx[2], float *z, float *pz);//获取直线与X轴的交点
/*最小二乘估计，y= b1*x+b0,估计b1 和b0 的值去拟合,这样，y=b1*x+b0是由参数x和y所指定的点集的最佳拟合线*/
void Lsqe(const float *x, const float *y, int n, float *b1, float* b0);
/*最小二乘估计，y= b1/x+b0,估计b1 和b0 的值去拟合,这样，y=b1/x+b0是由参数x和y所指定的点集的最佳拟合线*/
void LsqeAntyRateCurve(const float *x, const float *y, int n, float *b1, float* b0);

int	  MeanFilterInit(MEAN_VAL_FILTER_PARAM_T* parg, float* data_cache, int cache_size);//均值滤波
float MeanFilter(MEAN_VAL_FILTER_PARAM_T* parg, float samp); //
int   MedianFilterInit(MEAN_VAL_FILTER_PARAM_T* parg, float* data_cache, int cache_size); //中值滤波
float MedianFilterPutValue(MEAN_VAL_FILTER_PARAM_T* parg, float samp);//均值滤波的数据个数小于len时，返回samp
int	  MeanSquareParamInit(MEAN_SQUARE_PARAM_T *parg, float* data_cache,int cache_size); //求均方差初始化
float MeanSquareParamPutValue(MEAN_SQUARE_PARAM_T *parg, float ps); //求均方差

int	  RCLowPassFilterInit(RC_FILTER_PARAM_T* parg, float fs, float cutoff_fs, float rs);
float RCLowPassFilter(RC_FILTER_PARAM_T* parg, float data);

int   FirFilterInit(FIR_FILTER_PARAM_T* parg, const float* coef, float* Xn, int Order);
float FirFilterProcess(FIR_FILTER_PARAM_T* parg, float m_data);

int	  fst_order_highpass_filter_param_init(FST_ODR_HP_FIL_PARAM_T* parg, float radio);
float fst_order_highpass_filter(FST_ODR_HP_FIL_PARAM_T* parg, float samp);

void  BubbleSort(float a[], int n);//冒泡排序
void  ClearValueOfFtpt(FEATURE_POINT* pt1); //特征点的清零操作

void  GetFitedMinValueInWave(int dots, float buf[], float *z, float *pz);//获取波形顶点
void  GetFitedMaxValueInWave(int dots, float buf[], float *z, float *pz);//获取波形底点
void  GetZeroPointOfStrLineWithX_Axis(FEATURE_POINT* fp, float* z, float *pz);//获取直线与X轴相交的点
void  GetMaxPntOfTwoStrLineCross(int dots, float buf[], float* z, float *pz);//获取两相交直线的交点
 
#ifdef __cplusplus
}
#endif

#endif


