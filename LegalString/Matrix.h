#include <assert.h>
#include <map>
using namespace std;

const int MOD = 1000000007;
#define MAXN 100

/**
* 扩展欧几里德算法:求a模b的逆元
* @param a
* @param b
* @return
*/
__int64 extendGCD(__int64 a, __int64 b) {
    __int64 res = b;
    __int64 x1 = 1, x2 = 0;
    __int64 y1 = 0, y2 = 1;
    __int64 x = 0, y = 0;
    __int64 r = a % b;
    __int64 q = (a - r) / b;
    while (r > 0) {
        x = x1 - q * x2;
        y = y1 - q * y2;
        x1 = x2;
        y1 = y2;
        x2 = x;
        y2 = y;
        a = b;
        b = r;
        r = a % b;
        q = (a - r) / b;
    }

    if (x < 0) {
        res += x;
    } else {
        res = x;
    }
    return res;
}

__int64             g_C[100][100]           = {0};
map<int, __int64*>  map_Factor;
map<int, __int64*>  map_FactorOfKPow;
map<int, __int64*>  map_FactorOfPolynomial;
/**
* 计算排列组合(取模)
* 
* @param n
* @param m
* @return
*/
__int64 C(const __int64 n, const __int64 m) {
    assert(n >= m);
    assert(n <= 100);
    assert(m <= 100);
    if (0 != g_C[n][m])
    {
        return g_C[n][m];
    }
    __int64 i, a, b, p;
    if (n < m) {
        return 0; // n不能小于m
    }
    p = 1;
    i = n - m;
    a = i < m ? i : m;
    b = i > m ? i : m;
    __int64 middle = 0;
    for (i = 1; i <= a; i++) {
        middle = p * b;
        middle %= MOD;
        middle *= extendGCD(i, MOD);
        middle %= MOD;
        p += middle;
        p %= MOD;
    }
    g_C[n][m] = p;
    return p;
}

/**
* 计算(2x-1)^n的系数(幂次从0到n)
* 
* @param n
* return n+1 length array
*/
__int64* getFactor(const int n) {
    if (map_Factor.end() != map_Factor.find(n))
        return map_Factor[n];
    __int64 *fac = new __int64[n + 1];
    __int64 k = 1; // 2的幂(开始为2的0次幂)
    __int64 count = n; // -1的幂次
    int sign = 0;
    for (int i = 0; i <= n; ++i) {
        sign = count-- % 2 == 0 ? 1 : -1; // -1提供的系数
        fac[i] = (C(n , i) % MOD);// 组合数提供的系数
        fac[i] *= k; // 计算2的幂次
        fac[i] %= MOD;

        if(sign == -1){
            fac[i] = MOD - fac[i];
        }

        k *= 2;
        k %= MOD;
    }
    map_Factor[n] = fac;
    return fac;
}

/**
* 获取1到x的k次方幂的系数(结果数组一次代表: 幂次从1到(n + 1))
* @param k
* @return K+1 length array
*/
__int64* getFactorOfKPow(const int k){
    if (map_FactorOfKPow.end() != map_FactorOfKPow.find(k))
    {
        return map_FactorOfKPow[k];
    }
    __int64 *fac = new __int64[k + 1];
    fac[0] = extendGCD(k + 1 , MOD);

    __int64 mid = 0 , val = 0;
    for(int i = 1 ; i <= k ; ++ i){
        fac[i] = extendGCD(k - i + 1, MOD);
        mid = C(k , i);
        for(int j = 0 ; j < i ; ++ j){
            val = C(k - j + 1 , i - j + 1);
            val *= (MOD - fac[j]);
            val %= MOD;
            mid += val;
            mid %= MOD;
        }

        fac[i] *= mid;
        fac[i] = (fac[i] + MOD) % MOD;
    }

    //系数倒转,变成1次幂,2次幂...k次幂
    int len = k+1;
    __int64 *res = new __int64[len];
    for(int i = 0 ; i < len ; ++ i){
        res[len - 1 - i] = fac[i];
    }

    map_FactorOfKPow[k] = res;
    //delete fac;
    return res;
}

/**
* 获取1到(2 * x - 1)的k次幂和系数: 结果数组一次代表幂次从0到k+1
* @param k
* @return k+2 length array
*/
const __int64* getFactorOfPolynomial(const int k){
    if (map_FactorOfPolynomial.end() != map_FactorOfPolynomial.find(k))
    {
        return map_FactorOfPolynomial[k];
    }
    __int64 *res = new __int64[k + 2];

    for(int i = 0 ; i < k+2; ++i){ // 初始化置零
        res[i] = 0;
    }

    __int64 *powK = getFactorOfKPow(k);// 获取1到(k + 1)的系数

    __int64 *powF = NULL;
    __int64 mid = 0;
    for(int i = 0 ; i < k+1/*powK.length*/ ; ++ i){
        powF = getFactor(i+1);
        for(int j = 0 ; j < i+2/*powF.length*/ ; ++ j){
            mid = powF[j];
            mid *= powK[i];
            mid %= MOD;

            res[j] += mid;
            res[j] %= MOD;
        }
    }
    map_FactorOfPolynomial[k] = res;
    return res;
}

// 求M的K次方
int M_Power_K(int M, int K)
{
    if (0 == M)
        return 0;
    assert(K >= 0);
    __int64 nRet = 1;
    for (int i = 0; i < K; i++)
    {
        nRet *= M;
        nRet %= MOD;
    }
    return (int)nRet;
}

// M 个数，求对应的长度i为1到logM+1合法链的数量的表达式f(i)的系数，同时求出g(i)
__int64 f[100][100] = {0};
__int64 g[100]      = {0};
__int64 G[100]      = {0};
int get_g(int M)
{
    memset(f, 0, sizeof(f));
    memset(g, 0, sizeof(g));
    int nWidth = 1;
    f[1][0] = 1;
    g[1] = M;
    for (int idx = 2; true; idx++)
    {
        M = M/2;
        if (0 == M)
        {
            break;
        }
        nWidth++;
        // g[idx-1] 会影响 f[idx]的第1个系数（对应的是0阶）
        f[idx][0] = g[idx-1];
        for (int col = 0; col < idx-1; col++) // 循环f(idx-1)的系数个数次，确定f(idx)的系数，比f(idx-1)的系数多一个，多一个高阶
        {
            // f[idx] 的系数 = g[idx-1] - sigam(f(idx-1)的系数) 范围是(1, 2x-1);
            // col = 1, 第一个系数，对应的是0阶的系数
            
            // f[idx-1][0] = k0 会影响 f[idx]的第1个，第2个系数，怎么影响？-k0*sigma(1)[范围为1,2x-1]
            // 先求出sigma(1)[范围为1,2x-1], 其对应的是1阶多项式
            const __int64 *p = getFactorOfPolynomial(col);
            for (int i = 0; i < col+2; i++)
            {
                __int64 tmp = p[i] * f[idx-1][col];
                tmp %= MOD;
                f[idx][i] -= tmp;
                if (f[idx][i] < 0)
                {
                    f[idx][i] += MOD;
                }
            }
            //delete p;
        }
        // 系数算完之后，求g(idx) = sigma(f(idx))[区间范围是1到M](M已经除2了)
        for (int i = 0; i < idx; i++)
        {
            __int64* p = getFactorOfKPow(i);
            for (int j = 0; j < i+1; j++)
            {
                __int64 tmp = p[j]*M_Power_K(M, j+1);
                tmp %= MOD;
                assert(f[idx][i] < MOD);
                assert(f[idx][i] >= 0);
                tmp *= f[idx][i];
                tmp %= MOD;
                assert(tmp < MOD);
                assert(tmp >= 0);
                g[idx] += tmp;
                g[idx] %= MOD;
                assert(g[idx] < MOD);
                assert(g[idx] >= 0);
            }
        }
    }
    return nWidth;
}

int get_G(int M)
{
    int nWidth = get_g(M);
    __int64 g2[100] = {0};
    memcpy(g2, g, sizeof(g));
    memset(g, 0, sizeof(g));
    get_g(M/2);
    for (int i = 1; i <= nWidth; i++)
    {
        G[i] = g2[i] - g[i];
        if (G[i] < 0)
        {
            G[i] += MOD;
        }
    }
    return nWidth;
}

void ReleaseMemory()
{
    map<int, __int64*>::iterator itFac = map_Factor.begin();
    map<int, __int64*>::iterator itFacOfKP = map_FactorOfKPow.begin();
    map<int, __int64*>::iterator itFacOfPo = map_FactorOfPolynomial.begin();

    for ( ; itFac != map_Factor.end(); itFac++)
    {
        delete itFac->second;
    }
    map_Factor.clear();
    for ( ; itFacOfKP != map_FactorOfKPow.end(); itFacOfKP++)
    {
        delete itFacOfKP->second;
    }
    map_FactorOfKPow.clear();
    for ( ; itFacOfPo != map_FactorOfPolynomial.end(); itFacOfPo++)
    {
        delete itFacOfPo->second;
    }
    map_FactorOfPolynomial.clear();
}

// 矩阵快速幂乘
class Matrix { 
public: 

    __int64 m[MAXN][MAXN];
    int     m_nWidth;
    //二维数组存放矩阵 
    Matrix(){} 
    //对数组的初始化 
    void init(__int64  num[MAXN][MAXN], int nWidth){ 
        assert(nWidth <= MAXN);
        m_nWidth = nWidth;
        for(int i = 0 ; i < MAXN ; i++){ 
            for(int j = 0 ; j < MAXN ; j++){ 
                m[i][j] = num[i][j]; 
            } 
        } 
    } 

    //重载矩阵的乘法运算 
    friend Matrix operator*(Matrix &m1 ,Matrix &m2) { 
        int i, j, k; 
        Matrix temp; 
        for (i = 0; i < m1.m_nWidth; i++) { 
            for (j = 0; j < m1.m_nWidth; j++) { 
                temp.m[i][j] = 0; 
                for(k = 0 ; k < m1.m_nWidth ; k++) 
                    temp.m[i][j] += (m1.m[i][k] * m2.m[k][j])% MOD;
                    temp.m[i][j] %= MOD; 
                //注意每一步都进行取模 
            } 
        } 
        return temp; 
    } 
}; 

//矩阵的快速幂 
Matrix QuickPow(Matrix &M , __int64 n){ 
    Matrix tempans; 
    //初始化为单位矩阵 
    //初始化 
    for(int i = 0 ; i < M.m_nWidth ; i++){ 
        for(int j = 0 ; j < M.m_nWidth ; j++){ 
            if(i == j) 
                tempans.m[i][j] = 1; 
            else 
                tempans.m[i][j] = 0; 
        } 
    } 
    //快速幂（类似整数） 
    while(n){ 
        if(n & 1)
            tempans = tempans * M; 
        //已经重载了* 
        n = n >> 1; 
        M = M * M; 
    } 
    return tempans; 
} 

//求出V(0)到V(p-1)
void CalcV(__int64 V[], int nWidth)
{
    V[0] = 1;
    //V[1] = G[1];  
    for (int i = 1; i < nWidth; i++)
    {
        for (int j = 1; j <= i; j++)
        {
            V[i] += G[j]*V[i-j]%MOD;
            V[i] %= MOD;
        }
    }
}

// 
int get_validstring_count(int nNumCount, int nLength)
{
    int nWidth = get_G(nNumCount);
    __int64 V[MAXN] = {0};
    CalcV(V, nWidth);

    //构造矩阵
    __int64 matrix_arr[MAXN][MAXN] = {0};
    for (int i = 0; i < nWidth; i++)
    {
        matrix_arr[0][i] = G[i+1];
    }
    for (int i = 1; i < nWidth; i++)
    {
        matrix_arr[i][i-1] = 1;
    }

    Matrix matrix;
    matrix.init(matrix_arr, nWidth);

    Matrix newMatrix = QuickPow(matrix, nLength);



    return 0;

}