// Fortran Interface to STL vector class

// This class stores Fortran pointers converted to their integer representation

#include <iostream>
#include <vector>

#include "common.h"

#define VECTOR vector<void *>
#define VEC_ITERATOR VECTOR::iterator

typedef struct {
	void *ptr;
} POINTER_TO_POINTER;

using namespace std;

extern "C" PTR_TYPE CREATE_VECTOR();
extern "C" void DELETE_VECTOR(PTR_TYPE *vecint);
extern "C" void ADD_TO_VECTOR(PTR_TYPE *vecint, void *val);
extern "C" PTR_TYPE GET_VECTOR_BEGIN_ITERATOR(PTR_TYPE *mapint);
extern "C" PTR_TYPE GET_VECTOR_END_ITERATOR(PTR_TYPE *mapint);
extern "C" void VECTOR_ITERATE_NEXT(PTR_TYPE *iterint);
extern "C" void GET_VALUE_FOR_VECTOR_ITERATOR(PTR_TYPE *iterint, POINTER_TO_POINTER *valptr);
extern "C" bool COMPARE_VECTOR_ITERATORS(PTR_TYPE *iter1int, PTR_TYPE *iter2int);
extern "C" void DELETE_VECTOR_ITERATOR(PTR_TYPE *iterint);

PTR_TYPE CREATE_VECTOR()
{
	return reinterpret_cast<PTR_TYPE>(new VECTOR());
}

void DELETE_VECTOR(PTR_TYPE *vecint)
{
	VECTOR *vecptr = reinterpret_cast<VECTOR*>(*vecint);
	delete vecptr;
}

void ADD_TO_VECTOR(PTR_TYPE *vecint, void *val)
{
	VECTOR *vecptr = reinterpret_cast<VECTOR*>(*vecint);
	vecptr->push_back(val);
}

PTR_TYPE GET_VECTOR_BEGIN_ITERATOR(PTR_TYPE *vecint)
{
	VECTOR *vecptr = reinterpret_cast<VECTOR*>(*vecint);
	VEC_ITERATOR *begin = new VEC_ITERATOR(vecptr->begin());
	return reinterpret_cast<PTR_TYPE>(begin);
}

PTR_TYPE GET_VECTOR_END_ITERATOR(PTR_TYPE *vecint)
{
	VECTOR *vecptr = reinterpret_cast<VECTOR*>(*vecint);
	VEC_ITERATOR *end = new VEC_ITERATOR(vecptr->end());
	return reinterpret_cast<PTR_TYPE>(end);
}

void VECTOR_ITERATE_NEXT(PTR_TYPE *iterint)
{
	VEC_ITERATOR *iterptr = reinterpret_cast<VEC_ITERATOR*>(*iterint);
	++(*iterptr);
}

void GET_VALUE_FOR_VECTOR_ITERATOR(PTR_TYPE *iterint, POINTER_TO_POINTER *valptr)
{
	VEC_ITERATOR *iterptr = reinterpret_cast<VEC_ITERATOR*>(*iterint);
	valptr->ptr = **iterptr;
}

bool COMPARE_VECTOR_ITERATORS(PTR_TYPE *iter1int, PTR_TYPE *iter2int)
{
	VEC_ITERATOR *iter1ptr = reinterpret_cast<VEC_ITERATOR*>(*iter1int);
	VEC_ITERATOR *iter2ptr = reinterpret_cast<VEC_ITERATOR*>(*iter2int);
	return (*iter1ptr == *iter2ptr);
}

void DELETE_VECTOR_ITERATOR(PTR_TYPE *iterint)
{
	VEC_ITERATOR *iterptr = reinterpret_cast<VEC_ITERATOR*>(*iterint);
	delete iterptr;
}