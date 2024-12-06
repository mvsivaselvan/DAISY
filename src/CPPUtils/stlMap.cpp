// Fortran Interface to STL map class

// This class is an associative container for Fortran pointers. These pointers are converted to
// their integer representation and stored as <ID,pointer> pairs in the map

#include <iostream>
#include <map>

#include "common.h"

#define MAP map<ID_TYPE,void *>
#define ITERATOR MAP::iterator

typedef struct {
	void *ptr;
} POINTER_TO_POINTER;

using namespace std;

extern "C" PTR_TYPE CREATE_MAP();
extern "C" void DELETE_MAP(PTR_TYPE *mapint);
extern "C" void ADD_TO_MAP(PTR_TYPE *mapint, ID_TYPE *key, void *val);
extern "C" void REMOVE_FROM_MAP(PTR_TYPE *mapint, ID_TYPE *key);
extern "C" void GET_VAL_FOR_KEY(PTR_TYPE *mapint, ID_TYPE *key, POINTER_TO_POINTER *valptr);
extern "C" PTR_TYPE GET_BEGIN_ITERATOR(PTR_TYPE *mapint);
extern "C" PTR_TYPE GET_END_ITERATOR(PTR_TYPE *mapint);
extern "C" void ITERATE_NEXT(PTR_TYPE *iterint);
extern "C" void GET_VALUE_FOR_ITERATOR(PTR_TYPE *iterint, POINTER_TO_POINTER *valptr);
extern "C" bool COMPARE_ITERATORS(PTR_TYPE *iter1int, PTR_TYPE *iter2int);
extern "C" void DELETE_ITERATOR(PTR_TYPE *iterint);

PTR_TYPE CREATE_MAP()
{
	return reinterpret_cast<PTR_TYPE>(new MAP());
}

void DELETE_MAP(PTR_TYPE *mapint)
{
	MAP *mapptr = reinterpret_cast<MAP*>(*mapint);
	delete mapptr;
}

void ADD_TO_MAP(PTR_TYPE *mapint, ID_TYPE *key, void *val)
{
	MAP *mapptr = reinterpret_cast<MAP*>(*mapint);
	mapptr->insert(pair<ID_TYPE,void *>(*key,val));
}

void REMOVE_FROM_MAP(PTR_TYPE *mapint, ID_TYPE *key)
{
	MAP *mapptr = reinterpret_cast<MAP*>(*mapint);
	mapptr->erase(*key);
}

void GET_VAL_FOR_KEY(PTR_TYPE *mapint, ID_TYPE *key, POINTER_TO_POINTER *valptr)
{ 
	MAP *mapptr = reinterpret_cast<MAP*>(*mapint);
	valptr->ptr = (*mapptr)[*key];
}

PTR_TYPE GET_BEGIN_ITERATOR(PTR_TYPE *mapint)
{
	MAP *mapptr = reinterpret_cast<MAP*>(*mapint);
	ITERATOR *begin = new ITERATOR(mapptr->begin());
	return reinterpret_cast<PTR_TYPE>(begin);
}

PTR_TYPE GET_END_ITERATOR(PTR_TYPE *mapint)
{
	MAP *mapptr = reinterpret_cast<MAP*>(*mapint);
	ITERATOR *end = new ITERATOR(mapptr->end());
	return reinterpret_cast<PTR_TYPE>(end);
}

void ITERATE_NEXT(PTR_TYPE *iterint)
{
	ITERATOR *iterptr = reinterpret_cast<ITERATOR*>(*iterint);
	++(*iterptr);
}

extern "C" void GET_VALUE_FOR_ITERATOR(PTR_TYPE *iterint, POINTER_TO_POINTER *valptr)
{
	ITERATOR *iterptr = reinterpret_cast<ITERATOR*>(*iterint);
	valptr->ptr = (*iterptr)->second;
}

bool COMPARE_ITERATORS(PTR_TYPE *iter1int, PTR_TYPE *iter2int)
{
	ITERATOR *iter1ptr = reinterpret_cast<ITERATOR*>(*iter1int);
	ITERATOR *iter2ptr = reinterpret_cast<ITERATOR*>(*iter2int);
	return (*iter1ptr == *iter2ptr);
}

void DELETE_ITERATOR(PTR_TYPE *iterint)
{
	ITERATOR *iterptr = reinterpret_cast<ITERATOR*>(*iterint);
	delete iterptr;
}