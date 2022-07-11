#include <exception>

#pragma once
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

template<int I>
class fixed_point128
{
private:
    // members
    union
    {
        char bytes[16];
        struct
        {
            unsigned __int64 low;
            unsigned __int64 high;
            int sign; // 0 = positive, 1 negative
        };
    };
    static constexpr int int_bits = I;
    static constexpr int frac_bits = 128 - I;
    static constexpr int upper_frac_bits = (frac_bits <= 64) ? 0 : frac_bits - 64;
    static constexpr double upper_unity = 1.0 / (double)(1ull << (upper_frac_bits));
    static constexpr double lower_unity = upper_unity / (double)(1ull << (32)) / (double)(1ull << (32));
public:
    fixed_point128() { low = high = 0ull; sign = 0; }
    fixed_point128(double val) {
        unsigned __int64 i = *((unsigned __int64*)(&val));
        sign = GET_BIT(i, 63);
        int e = GET_BITS(i, 52, 11) - 1023;
        unsigned __int64 f = (i & MAX_BITS_VALUE(52));
        if (e >= int_bits) {
            throw std::exception("overflow!");
        }
        // normal number
        if (e > -1023) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= ONE_SHIFT(52);
            int bits_to_shift = 64 - int_bits - 52 + e;

            // f fits in high QWORD
            if (bits_to_shift >= 0) {
                high = f << bits_to_shift;
                low = 0;
            }
            // shift right
            else {
                bits_to_shift = -bits_to_shift;
                // f has some bits in high QWORD
                if (bits_to_shift <= 53) {
                    high = f >> bits_to_shift;
                    low = GET_BITS(f, 0, bits_to_shift);
                    low <<= (64 - bits_to_shift);
                }
                // shift f into low QWORD
                else {
                    high = 0;
                    bits_to_shift -= 52;
                    f <<= (63 - 52);
                    low = f >> bits_to_shift;
                }
            }
        }
    }

    // operators
    inline operator unsigned int() { return (int)(high >> (64 - int_bits)) & ((unsigned)(-1)); }
    inline operator int() { return (int)((__int64)high >> (64 - int_bits)) & ((unsigned)(-1)); }
    inline operator double() {
        double res;
        // common case
    #if int_bits <= 64
        if (int_bits <= 64) {
            res = (double)(high * upper_unity); // high 64 bit part - bits 64-127
            res += (double)(low * lower_unity);
            // do the rest with 32 bit
            //res += lower_unity1 * (double)(low >> 32);            // bits 32-63
            //res += lower_unity2 * (double)(low & (unsigned)(-1)); // bits 0-31
        }
    #else
        else if (frac_bits < 64) {

            res = ((double)high / (double)(1ull << frac_bits));
            //res += 
        }
    #endif
        return res;
    }
    //
};
