//
// hashMap.c
// CellCoal
//
// Created by Senbai Kang on 11/20/19.
//

#include "hashMap.h"

struct HashMapIntToIntArr *CreateHashMapIntToIntArr(int size) {
  struct HashMapIntToIntArr *m;

  m = (struct HashMapIntToIntArr *)malloc(sizeof(struct HashMapIntToIntArr));
  m->size = size;
  m->list = (struct NodeIntToIntArr **)malloc((size + 1) *
                                              sizeof(struct NodeIntToIntArr *));
  for (int i = 0; i < size; i++) {
    m->list[i] = NULL;
  }

  return m;
}

int HashCodeIntToIntArr(const struct HashMapIntToIntArr *hashMap, int key) {
  if (key < 0)
    return -key % hashMap->size;
  return key % hashMap->size;
}

int *GetValIntToIntArr(const struct HashMapIntToIntArr *hashMap, int key) {
  int location;
  struct NodeIntToIntArr *itr;

  location = HashCodeIntToIntArr(hashMap, key);
  itr = hashMap->list[location];

  while (itr) {
    if (itr->key == key) {
      return itr->val;
    }
    itr = itr->next;
  }

  fprintf(stderr,
          "\nWARNING! KEY NOT FOUND: Key %d is not in the hash map. \n\n", key);
  return NULL;
}

void InsertKeyValIntToIntArr(struct HashMapIntToIntArr *hashMap, int key,
                             int *val) {
  int location;
  struct NodeIntToIntArr *targetList, *itr, *newNode;

  location = HashCodeIntToIntArr(hashMap, key);
  itr = targetList = hashMap->list[location];

  while (itr) {
    if (itr->key == key) {
      fprintf(
          stderr,
          "\nWARNING! DUPLICATE KEY: Key %d is already in the hash map. \n\n",
          key);
      return;
    }
    itr = itr->next;
  }

  newNode = (struct NodeIntToIntArr *)malloc(sizeof(struct NodeIntToIntArr));

  newNode->key = key;
  newNode->val = val;
  newNode->next = targetList;

  hashMap->list[location] = newNode;
}

void DeleteKeyIntToIntArr(struct HashMapIntToIntArr *hashMap, int key) {
  int location;
  struct NodeIntToIntArr *itr, *paNode = NULL;

  location = HashCodeIntToIntArr(hashMap, key);
  itr = hashMap->list[location];

  while (itr) {
    if (itr->key == key) {
      if (paNode == NULL) {
        hashMap->list[location] = itr->next;
      } else {
        paNode->next = itr->next;
      }
      free(itr);
      return;
    } else {
      paNode = itr;
      itr = itr->next;
    }
  }

  fprintf(stderr,
          "\nWARNING! KEY NOT FOUND: Key %d is not in the hash map. \n\n", key);
}

struct HashMapIntArrToInt *CreateHashMapIntArrToInt(int size) {
  struct HashMapIntArrToInt *m;

  m = (struct HashMapIntArrToInt *)malloc(sizeof(struct HashMapIntArrToInt));
  m->size = size;
  m->list = (struct NodeIntArrToInt **)malloc((size + 1) *
                                              sizeof(struct NodeIntArrToInt *));
  for (int i = 0; i < size; i++) {
    m->list[i] = NULL;
  }

  return m;
}

long ToBinary(int a) {
  long b = 0;
  int reminder, i = 1;

  while (a != 0) {
    reminder = a % 2;
    a /= 2;
    b += reminder * i;
    i *= 10;
  }

  return b;
}

long GetBinaryKey(int *key) {
  long binaryKey = 0;
  int *tmpPtr = key;

  for (int i = 0; i < 3; ++i) {
    binaryKey +=
        (ToBinary((*tmpPtr >= 0) ? *tmpPtr : ((*tmpPtr == -1) ? 4 : 5)) *
         (long)powl(1000, 2 - i));
    tmpPtr++;
  }

  return binaryKey;
}

int GetDecimalKey(long key) {
  int decimalKey = 0, i = 0;
  long reminder;

  while (key != 0) {
    reminder = key % 10;
    key /= 10;
    decimalKey += reminder * (int)pow(2, i);
    ++i;
  }

  return decimalKey;
}

int HashCodeIntArrToInt(const struct HashMapIntArrToInt *hashMap, int *key) {
  int decimalKey;
  long binaryKey;

  binaryKey = GetBinaryKey(key);
  decimalKey = GetDecimalKey(binaryKey);

  return decimalKey % hashMap->size;
}

int KeySanityCheck(int *tmpKey, int *key) {
  int *tmpKeyPtr = tmpKey;
  int *keyPtr = key;
  int noNegativeLociNr = 0;

  for (int i = 0; i < 3; ++i) {
    *tmpKeyPtr = *keyPtr;

    if (*keyPtr >= 0) {
      noNegativeLociNr++;
    }

    tmpKeyPtr++;
    keyPtr++;
  }

  if (noNegativeLociNr == 3) {
    qsort(tmpKey, 3, sizeof(int), GreaterThanInt);
  } else if (noNegativeLociNr == 2) {
    qsort(tmpKey, 3, sizeof(int), LessThanInt);
    qsort(tmpKey, 2, sizeof(int), GreaterThanInt);
  } else if (noNegativeLociNr == 1) {
    qsort(tmpKey, 3, sizeof(int), LessThanInt);
  }

  for (int i = 0; i < 3; ++i) {
    if (tmpKey[i] == ADOExtd) {
      if (i != 2) {
        tmpKey[i] = DELETED;
      } else {
        tmpKey[i] = INVALID;
      }
    } else if (tmpKey[i] == INVALID) {
      if (i != 2) {
        tmpKey[i] = DELETED;
      }
    } else if (tmpKey[i] == DELETED) {
      if (i == 2) {
        tmpKey[i] = INVALID;
      }
    }
  }

  return noNegativeLociNr;
}

/*
 * matchStyle: match the INVALID and deleted allele rigorously or not.
 * matchStyle == 1: rigorously match
 * matchStyle == 0: not rigorously match (do not check whether two numbers are
 * the same or not when they are both negative)
 */
int GetValIntArrToInt(const struct HashMapIntArrToInt *hashMap, int *key,
                      int matchStyle) {
  int location, noNegativeLociNr = 3;
  struct NodeIntArrToInt *itr;
  int tmpKey[3] = {0, 0, 0};

  if (matchStyle == 0)
    noNegativeLociNr = KeySanityCheck(tmpKey, key);

  location = HashCodeIntArrToInt(hashMap, tmpKey);
  itr = hashMap->list[location];

  while (itr) {
    int matchNr = 0;
    for (int i = 0; i < 3; ++i) {
      if ((itr->key[i] == tmpKey[i]) &&
          (((itr->key[i] >= 0) && (tmpKey[i] >= 0)) ||
           ((itr->key[i] < 0) && (tmpKey[i] < 0) && matchStyle))) {
        matchNr++;
      }
    }

    if (matchNr == noNegativeLociNr) {
      return itr->val;
    } else {
      itr = itr->next;
    }
  }

  fprintf(
      stderr,
      "\nWARNING! KEY NOT FOUND: Key (%d, %d, %d) is not in the hash map. \n\n",
      tmpKey[0], tmpKey[1], tmpKey[2]);
  return -1;
}

void InsertKeyValIntArrToInt(struct HashMapIntArrToInt *hashMap, int *key,
                             int val) {
  int location;
  struct NodeIntArrToInt *targetList, *itr, *newNode;

  location = HashCodeIntArrToInt(hashMap, key);
  itr = targetList = hashMap->list[location];

  while (itr) {
    if ((itr->key[0] == key[0]) && (itr->key[1] == key[1]) &&
        (itr->key[2] == key[2])) {
      fprintf(stderr,
              "\nWARNING! DUPLICATE KEY: Key (%d, %d, %d) is already in the "
              "hash map. \n\n",
              key[0], key[1], key[2]);
      return;
    }
    itr = itr->next;
  }

  newNode = (struct NodeIntArrToInt *)malloc(sizeof(struct NodeIntArrToInt));

  newNode->key = key;
  newNode->val = val;
  newNode->next = targetList;

  hashMap->list[location] = newNode;
}

void DeleteKeyIntArrToInt(struct HashMapIntArrToInt *hashMap, int *key) {
  int location;
  struct NodeIntArrToInt *itr, *paNode = NULL;

  location = HashCodeIntArrToInt(hashMap, key);
  itr = hashMap->list[location];

  while (itr != NULL) {
    if ((itr->key[0] == key[0]) && (itr->key[1] == key[1]) &&
        (itr->key[2] == key[2])) {
      if (paNode == NULL) {
        hashMap->list[location] = itr->next;
      } else {
        paNode->next = itr->next;
      }
      free(itr);
      return;
    } else {
      paNode = itr;
      itr = itr->next;
    }
  }

  fprintf(
      stderr,
      "\nWARNING! KEY NOT FOUND: Key (%d, %d, %d) is not in the hash map. \n\n",
      key[0], key[1], key[2]);
}
