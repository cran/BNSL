#ifndef ASTER_HPP
#define ASTER_HPP

#include <sys/time.h>
#include <unistd.h>

//typedef long long int int64_t;
//typedef unsinged long long int uint64_t;
typedef uint64_t vset;
typedef double score_t;

//const score_t INF = 99999999.0;

inline double sec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

inline int CountBit(uint64_t x)
{
    // only for gcc
    return __builtin_popcountll(x);
}

#endif // ASTER_HPP
