#include <assert.h>

const int MOD = 1000000007;

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

/**
* 计算排列组合(取模)
* 
* @param n
* @param m
* @return
*/
__int64 C(__int64 n, __int64 m) {
    assert(n >= m);
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
    return p;
}

/**
* 计算(2x-1)^n的系数(幂次从0到n)
* 
* @param n
*/
__int64* getFactor(int n) {
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
    return fac;
}

/**
* 获取1到x的k次方幂的系数(结果数组一次代表: 幂次从1到(n + 1))
* @param k
* @return
*/
__int64* getFactorOfKPow(int k){
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
    k ++;
    __int64 *res = new __int64[k];
    for(int i = 0 ; i < k ; ++ i){
        res[k - 1 - i] = fac[i];
    }

    delete fac;
    return res;
}

/**
* 获取1到(2 * x - 1)的k次幂和系数: 结果数组一次代表幂次从0到k+1
* @param k
* @return
*/
__int64* getFactorOfPolynomial(int k){
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
    return res;
}