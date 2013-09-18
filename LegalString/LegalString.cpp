// LegalString.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <Windows.h>
#include <winnt.h>
#include "Matrix.h"
#include <conio.h>

//#define PRINT

// n is id, m is length
int validstring_1(int n,int m) 
{
    if (0 == n || 0 == m)
    {
        return 0;
    }
    int *pCurr = new int[n+1];
    int *pNext = new int[n+1];
    for (int i = 1; i < n+1; i++)
    {
        pCurr[i] = 1;
    }
    for (int i = 1; ; i++)
    {
#ifdef PRINT
        for (int ii = 1; ii < n; ii++)
        {
            printf("%d(%d) ", pCurr[ii], pCurr[ii+1]-pCurr[ii]);
        }
        printf("%d\n", pCurr[n]);
#endif
        if (i >= m) break;
        for (int i = 1; i < n+1; i++)
        {
            pNext[i] = 0;
        }
        for (int k = 1; k <= n; k++)
        {
            int nNum = pCurr[k];
            int j = (2*k > n) ? 1 : 2*k;
            for ( ; j <= n; j++)
            {
                pNext[j] += nNum;
                if (pNext[j] >= 1000000007) 
                    pNext[j] -= 1000000007;

            }
        }
        int *pTemp = pCurr;
        pCurr = pNext;
        pNext = pTemp;
    }
    int sum = 0;
    for (int i = n/2 + 1; i <= n; i++)
    {
        sum += pCurr[i];
        if (sum >= 1000000007) 
            sum -= 1000000007;
    }
    delete pCurr;
    delete pNext;
    return sum;
}


 int __stdcall validstring(int n,int m) 
{
    if (0 == n || 0 == m)
    {
        return 0;
    }
    int *pCurr = new int[n+1];
    int *pNext = new int[n+1];
    for (int i = 1; i < n+1; i++)
    {
        pCurr[i] = 1;
    }
    int nCurrHalfEnd = n - n/2;
    int nCurrFull = n;
    for (int i = 1; ; i++)
    {
#ifdef PRINT
        for (int ii = 1; ii < n; ii++)
        {
            printf("%4d ", pCurr[ii]);
        }
        printf("%4d\n", pCurr[n]);
#endif
        if (i > m) break;
        // set first value
        pNext[1] = nCurrHalfEnd;
        nCurrHalfEnd = 0;
        for (int i = 2; i < n+1; i++)
        {
            if ( i&1 )
                pNext[i] = pNext[i-1];
            else
            {
                pNext[i] = pNext[i-1] + pCurr[i/2];
                if (pNext[i] >= 1000000007) 
                    pNext[i] -= 1000000007;
            }
            if (2*i > n) 
            {
                nCurrHalfEnd += pNext[i];
                if (nCurrHalfEnd >= 1000000007) 
                    nCurrHalfEnd -= 1000000007;
            }
        }
        int *pTemp = pCurr;
        pCurr = pNext;
        pNext = pTemp;
    }
    int nRet = pCurr[1];
    delete pCurr;
    delete pNext;
    return nRet;
}

int GetNext(int *pCurr, int nNumCount, int nState)
{
    __int64 nRet = 0;
    int nLeft = (nState/nNumCount)+1;
    int nRight = (nState%nNumCount)+1;

    for (int i = 1; i <= nNumCount; i++) //左半边结尾的字符
    {
        int nLeftCount = pCurr[(nLeft-1)*nNumCount + (i-1)];
        int nRightCount = 0;
        int j = (2*i > nNumCount) ? 1 : 2*i; //右半边开头的字符
        for ( ; j <= nNumCount; j++)
        {
            nRightCount += pCurr[(j-1)*nNumCount + (nRight-1)];
            nRightCount %= 1000000007;
        }
        nRet += (__int64)nLeftCount * nRightCount;
        nRet %= 1000000007;
    }
    return (int)nRet;
}

int validstring_2(int nNumCount, int nLength)
{
    if (0 == nNumCount || 0 == nLength)
    {
        return 0;
    }
    if (1 == nLength)
    {
        return (1+nNumCount)/2;
    }
    int k = 1;
    int nCount = 0;
    while (k < nLength) 
    { 
        k = k<<1;
        nCount++;
    }
    k = k>>1;
    printf("%d\n", k);

    int *pBuff = new int[(nNumCount)*(nNumCount)*(nCount)];
    int *pCurr = pBuff;
    for (int i = 0; i < (nNumCount)*(nNumCount); i++)
    {
        int nLeft = (i/nNumCount)+1;
        int nRight = (i%nNumCount)+1;
        pCurr[i] = (nRight >= 2*nLeft || nLeft*2 > nNumCount) ? 1 : 0;
    }
    for (int nIdxLen = 2; nIdxLen <= k; nIdxLen = nIdxLen<<1)//nIdxLen 是指数 = 2 表示长度为4
    {
        int *pNext = &pCurr[nNumCount*nNumCount];
        for (int i = 0; i < (nNumCount)*(nNumCount); i++)
        {
            pNext[i] = GetNext(pCurr, nNumCount, i);
        }
        pCurr = pNext;
    }
    int nRet = 0;
    for (int i = 1; i < nNumCount*nNumCount; i++)
    {
        int nLeft = (i/nNumCount)+1;
        int nRight = (i%nNumCount)+1;
        if (nRight*2 > nNumCount)
        {
            nRet += pCurr[i];
            nRet %= 1000000007;
        }
    }
    delete pBuff;
    return nRet;
}

int _tmain(int argc, _TCHAR* argv[])
{
    get_G(1000000000);

    assert(M_Power_K(1,3) == 1);
    assert(M_Power_K(2,3) == 8);
    assert(M_Power_K(2,0) == 1);
    assert(M_Power_K(2,10) == 1024);
    assert(M_Power_K(10,4) == 10000);
    assert(M_Power_K(0,0) == 0);
    assert(M_Power_K(100,0) == 1);
//    printf("%d\n", validstring_1(3, 3));
//    printf("%d\n", validstring(3, 3));
//    printf("%d\n", validstring(5, 5));
//    printf("%d\n", validstring(4, 2));
//    printf("%d\n", validstring(4, 3));
//    printf("%d\n", validstring(4, 4));
//    printf("%d\n", validstring(4, 5));
//    printf("%d\n", validstring(4, 6));
//    printf("%d\n", validstring(4, 7));
//    printf("%d\n", validstring(4, 8));
//    printf("%d\n", validstring(7, 7));
    printf("%d\n", validstring_2(10, 1024*1024*128));
    printf("%d\n", validstring(10, 2));
    printf("%d\n", validstring(10, 4));
    printf("%d\n", validstring(10, 8));
    printf("%d\n", validstring(10, 16));
    printf("%d\n", validstring(10, 32));
    printf("%d\n", validstring(10, 64));
    printf("%d\n", validstring(10, 128));
    printf("%d\n", validstring(10, 256));
    printf("%d\n", validstring(10, 256));
    printf("%d\n", validstring(10, 256));

    printf("%d\n", validstring(4, 1));
    printf("%d\n", validstring(4, 2));
    printf("%d\n", validstring(4, 3));
    printf("%d\n", validstring(4, 4));
//    printf("%d\n", validstring_1(10, 536870912));
//    printf("%d\n", validstring(10, 10));
//    printf("%d\n", validstring(6, 6));
//    printf("%d\n", validstring(66, 634));
    _getch();

  	return 0;
}

