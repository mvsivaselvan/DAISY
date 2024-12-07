// Fortran Interface to STL map class

// This class is an associative container for Fortran pointers. These pointers are converted to
// their integer representation and stored as <ID,pointer> pairs in the map

#include <iostream>
#include <map>

#include "common.h"

#define MAP map<ID_TYPE,void *>
#define ITERATOR MAP::iterator

using namespace std;

extern "C" MAP* CREATE_MAP();
extern "C" void DELETE_MAP(MAP *mapptr);
extern "C" void ADD_TO_MAP(MAP *mapptr, ID_TYPE key, void *val);
extern "C" void REMOVE_FROM_MAP(MAP *mapptr, ID_TYPE key);
extern "C" void* GET_VAL_FOR_KEY(MAP *mapptr, ID_TYPE key);
extern "C" ITERATOR* GET_BEGIN_ITERATOR(MAP *mapptr);
extern "C" ITERATOR* GET_END_ITERATOR(MAP *mapptr);
extern "C" void ITERATE_NEXT(ITERATOR *iterptr);
extern "C" void* GET_VALUE_FOR_ITERATOR(ITERATOR *iterptr);
extern "C" bool COMPARE_ITERATORS(ITERATOR *iter1ptr, ITERATOR *iter2ptr);
extern "C" void DELETE_ITERATOR(ITERATOR *iterptr);

MAP* CREATE_MAP()
{
	return new MAP();
}

void DELETE_MAP(MAP *mapptr)
{
	if (mapptr != NULL) delete mapptr;
}

void ADD_TO_MAP(MAP *mapptr, ID_TYPE key, void *val)
{
	mapptr->insert(pair<ID_TYPE,void *>(key,val));
}

void REMOVE_FROM_MAP(MAP *mapptr, ID_TYPE key)
{
	mapptr->erase(key);
}

void* GET_VAL_FOR_KEY(MAP *mapptr, ID_TYPE key)
{ 
	return (void *)((*mapptr)[key]);
}

ITERATOR* GET_BEGIN_ITERATOR(MAP *mapptr)
{
	return new ITERATOR(mapptr->begin());
}

ITERATOR* GET_END_ITERATOR(MAP *mapptr)
{
	return new ITERATOR(mapptr->end());
}

void ITERATE_NEXT(ITERATOR *iterptr)
{
	++(*iterptr);
}

extern "C" void* GET_VALUE_FOR_ITERATOR(ITERATOR *iterptr)
{
	return (void *)((*iterptr)->second);
}

bool COMPARE_ITERATORS(ITERATOR *iter1ptr, ITERATOR *iter2ptr)
{
	return (*iter1ptr == *iter2ptr);
}

void DELETE_ITERATOR(ITERATOR *iterptr)
{
	if (iterptr != NULL) delete iterptr;
}