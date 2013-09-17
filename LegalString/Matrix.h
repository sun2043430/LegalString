#include <assert.h>

const int MOD = 1000000007;

/**
* ��չŷ������㷨:��aģb����Ԫ
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
* �����������(ȡģ)
* 
* @param n
* @param m
* @return
*/
__int64 C(__int64 n, __int64 m) {
    assert(n >= m);
    __int64 i, a, b, p;
    if (n < m) {
        return 0; // n����С��m
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
* ����(2x-1)^n��ϵ��(�ݴδ�0��n)
* 
* @param n
*/
__int64* getFactor(int n) {
    __int64 *fac = new __int64[n + 1];
    __int64 k = 1; // 2����(��ʼΪ2��0����)
    __int64 count = n; // -1���ݴ�
    int sign = 0;
    for (int i = 0; i <= n; ++i) {
        sign = count-- % 2 == 0 ? 1 : -1; // -1�ṩ��ϵ��
        fac[i] = (C(n , i) % MOD);// ������ṩ��ϵ��
        fac[i] *= k; // ����2���ݴ�
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
* ��ȡ1��x��k�η��ݵ�ϵ��(�������һ�δ���: �ݴδ�1��(n + 1))
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

    //ϵ����ת,���1����,2����...k����
    k ++;
    __int64 *res = new __int64[k];
    for(int i = 0 ; i < k ; ++ i){
        res[k - 1 - i] = fac[i];
    }

    delete fac;
    return res;
}

/**
* ��ȡ1��(2 * x - 1)��k���ݺ�ϵ��: �������һ�δ����ݴδ�0��k+1
* @param k
* @return
*/
__int64* getFactorOfPolynomial(int k){
    __int64 *res = new __int64[k + 2];

    for(int i = 0 ; i < k+2; ++i){ // ��ʼ������
        res[i] = 0;
    }

    __int64 *powK = getFactorOfKPow(k);// ��ȡ1��(k + 1)��ϵ��

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