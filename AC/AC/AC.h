#ifndef _AC_H_
#define _AC_H_
#include "def.h"
#include "IO_AC.h"
#include "Str_AC.h"


#define REQUIRE_INPUTS(N) \
if (argc < N + 1){\
	printf("Need %d input files.\n", N);\
	return -1;\
}\

#ifdef PRINT_NODE_TRUE 
#define PRINT_NODE(a) printf(a);
#else
#define PRINT_NODE(a)
#endif

#define ForInt(IDX, RangeBegin, RangeEnd) for(int IDX = RangeBegin; IDX < RangeEnd; ++IDX)
#define For0Int(IDX, RangeEnd) ForInt(IDX, 0, RangeEnd)
#endif
