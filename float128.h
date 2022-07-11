#pragma once

#include <exception>
#include <string>
#include <limits>

using std::exception;

#ifndef ONE_SHIFT
#define ONE_SHIFT(x)        ( 1ull << (x) )
#endif

#ifndef MAX_BITS_VALUE
#define MAX_BITS_VALUE(x)   ( ONE_SHIFT(x) - 1 )
#endif

#ifndef GET_BIT
#define GET_BIT(x, n)        ( ((x) >> n) & 1)
#endif

#ifndef GET_BITS
#define GET_BITS(x, b, count)        ( ((x) >> (b)) & MAX_BITS_VALUE(count))
#endif

typedef union _float32
{
    struct
    {
        unsigned int f : 23;
        unsigned int e : 8;
        unsigned int s : 1;
    };
    float native;
    unsigned int integer;
} float32_t;

typedef union _float64
{
    struct
    {
        unsigned long long f : 52;
        unsigned int e : 11;
        unsigned int s : 1;
    };
    double native;
    unsigned __int64 integer;
} float64_t;

class not_implemented : public exception
{
public:
    not_implemented(const char* msg=NULL) : exception(msg) {
    }
};

class float128
{
public:
    // ctors
    float128() {
        e = 0; s = 0; f[0] = 0; f[1] = 0; // set to zero
    }
    float128(double val) {
    #ifdef DEBUG
        float64_t* v = (float64_t*)(&val); // for debug
    #endif
        unsigned __int64 i = *((unsigned __int64*)(&val));
        s = GET_BIT(i, 63);
        e = GET_BITS(i, 52, 11) - 1023; 
        f[0] = 0ull;
        f[1] = (i & MAX_BITS_VALUE(52)) << (64 - 52);
    }
    float128(float val) {
        unsigned __int64 i = *((unsigned*)(&val));
        s = GET_BIT(i, 31);
        e = GET_BITS(i, 23, 8) - 127;
        f[0] = 0ull;
        f[1] = (i & MAX_BITS_VALUE(52)) << (64 - 23);
    }
    float128(int val) { //TODO:
        throw not_implemented();
    }
    float128(unsigned int val) { //TODO:
        throw not_implemented();
    }
    float128(long long int val) { //TODO:
        throw not_implemented();
    }
    float128(unsigned long long int val) { //TODO:
        throw not_implemented();
    }
    float128(const char* val) { //TODO:
        throw not_implemented();
    }
    //dtor
    ~float128() {}
    
    //conversions
    operator double() { 
        if (e > 1023)
            return INFINITY;
        else if (e < -1022)
            return 0.0;

        // TODO: handle NaN
        double res = 0;
        unsigned __int64* p = (unsigned __int64*)(&res);
        *p = ((s & 1ull) << 63) | ((e + 1023ull) << 52) | (f[1] >> (64 - 52));

        return res;
    }
    operator float() {
        if (e > 127)
            return INFINITY;
        else if (e < -126)
            return 0.0;
        
        float res = 0.0f;
        unsigned int* p = (unsigned int*)(&res);
        *p = ((s & 1ul) << 31) | ((e + 127ul) << 23) | (unsigned)((f[1] >> (64 - 23)));
        return res;
    }
    operator std::string() { //TODO:
        throw not_implemented();
    }
private:
    unsigned int s; //sign
    int e; //exponent
    unsigned long long f[2]; // MSB is 0 and LSB is 1
};