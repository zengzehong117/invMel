#pragma once
 
#define      MAX_MATRIX_SIZE                   4194304             // 2048 * 2048
//#define      PI                                           3.141592653
#define      MAX_VECTOR_LENGTH              10000             //
 
typedef struct Complex
{
	float rl;
	float im;
}Complex;
 
class CFft1
{
public:
	CFft1(void){}
	~CFft1(void){}
 
public:
	bool fft(Complex const inVec[], int const vecLen, Complex outVec[]){
		char msg[256] = "11111 ";

		if ((vecLen <= 0) || (NULL == inVec) || (NULL == outVec))
			return false;
		if (!is_power_of_two(vecLen))
			return false;

		// create the weight array
		Complex         *pVec = new Complex[vecLen];
		Complex         *Weights = new Complex[vecLen];
		Complex         *X = new Complex[vecLen];
		int                   *pnInvBits = new int[vecLen];

		memcpy(pVec, inVec, vecLen*sizeof(Complex));

		// 计算权重序列
		float fixed_factor = (-2 * PI) / vecLen;
		for (int i = 0; i < vecLen / 2; i++) {
			float angle = i * fixed_factor;
			Weights[i].rl = cos(angle);
			Weights[i].im = sin(angle);
		}
		for (int i = vecLen / 2; i < vecLen; i++) {
			Weights[i].rl = -(Weights[i - vecLen / 2].rl);
			Weights[i].im = -(Weights[i - vecLen / 2].im);
		}

		int r = get_computation_layers(vecLen);

		// 计算倒序位码
		int index = 0;
		for (int i = 0; i < vecLen; i++) {
			index = 0;
			for (int m = r - 1; m >= 0; m--) {
				index += (1 && (i & (1 << m))) << (r - m - 1);
			}
			pnInvBits[i] = index;
			X[i].rl = pVec[pnInvBits[i]].rl;
			X[i].im = pVec[pnInvBits[i]].im;
		}

		// 计算快速傅里叶变换
		for (int L = 1; L <= r; L++) {
			int distance = 1 << (L - 1);
			int W = 1 << (r - L);

			int B = vecLen >> L;
			int N = vecLen / B;

			for (int b = 0; b < B; b++) {
				int mid = b*N;
				for (int n = 0; n < N / 2; n++) {
					int index = n + mid;
					int dist = index + distance;
					pVec[index].rl = X[index].rl + (Weights[n*W].rl * X[dist].rl - Weights[n*W].im * X[dist].im);                      // Fe + W*Fo
					pVec[index].im = X[index].im + Weights[n*W].im * X[dist].rl + Weights[n*W].rl * X[dist].im;
				}
				for (int n = N / 2; n < N; n++) {
					int index = n + mid;
					pVec[index].rl = X[index - distance].rl + Weights[n*W].rl * X[index].rl - Weights[n*W].im * X[index].im;        // Fe - W*Fo
					pVec[index].im = X[index - distance].im + Weights[n*W].im * X[index].rl + Weights[n*W].rl * X[index].im;
				}
			}

			memcpy(X, pVec, vecLen*sizeof(Complex));
		}

		memcpy(outVec, pVec, vecLen*sizeof(Complex));
		if (Weights)      delete[] Weights;
		if (X)                 delete[] X;
		if (pnInvBits)    delete[] pnInvBits;
		if (pVec)           delete[] pVec;
		return true;
	}            // 基于蝶形算法的快速傅里叶变换
    bool ifft(Complex const inVec[], int const len, Complex outVec[]);
 
	bool is_power_of_two(int num){
		int temp = num;
		int mod = 0;
		int result = 0;

		if (num < 2)
			return false;
		if (num == 2)
			return true;

		while (temp > 1)
		{
			result = temp / 2;
			mod = temp % 2;
			if (mod)
				return false;
			if (2 == result)
				return true;
			temp = result;
		}
		return false;
	}
	int    get_computation_layers(int num){
		int nLayers = 0;
		int len = num;
		if (len == 2)
			return 1;
		while (true) {
			len = len / 2;
			nLayers++;
			if (len == 2)
				return nLayers + 1;
			if (len < 1)
				return -1;
		}
	}         // calculate the layers of computation needed for FFT
};
