//
// hashMap.h
// CellCoal
//
// Created by Senbai Kang on 11/20/19.
//

#ifndef MYSIMULATOR_HASHMAP_H
#define MYSIMULATOR_HASHMAP_H

#include "extraDefinitions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct NodeIntToIntArr {
  int key;
  int *val;
  struct NodeIntToIntArr *next;
};

struct HashMapIntToIntArr {
  int size;
  struct NodeIntToIntArr **list;
};

struct HashMapIntToIntArr *CreateHashMapIntToIntArr(int size);

int HashCodeIntToIntArr(const struct HashMapIntToIntArr *hashMap, int key);

int *GetValIntToIntArr(const struct HashMapIntToIntArr *hashMap, int key);

void InsertKeyValIntToIntArr(struct HashMapIntToIntArr *hashMap, int key,
                             int *val);

void DeleteKeyIntToIntArr(struct HashMapIntToIntArr *hashMap, int key);

struct NodeIntArrToInt {
  int *key;
  int val;
  struct NodeIntArrToInt *next;
};

struct HashMapIntArrToInt {
  int size;
  struct NodeIntArrToInt **list;
};

struct HashMapIntArrToInt *CreateHashMapIntArrToInt(int size);

long ToBinary(int a);

long GetBinaryKey(int *key);

int GetDecimalKey(long key);

int HashCodeIntArrToInt(const struct HashMapIntArrToInt *hashMap, int *key);

int KeySanityCheck(int *tmpKey, int *key);

int GetValIntArrToInt(const struct HashMapIntArrToInt *hashMap, int *key,
                      int matchStyle);

void InsertKeyValIntArrToInt(struct HashMapIntArrToInt *hashMap, int *key,
                             int val);

void DeleteKeyIntArrToInt(struct HashMapIntArrToInt *hashMap, int *key);

int GreaterThanInt(const void *a, const void *b);

int LessThanInt(const void *a, const void *b);

#endif // MYSIMULATOR_HASHMAP_H
