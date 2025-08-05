//
// Created by Senbai Kang on 11/24/19.
//

#include <unity.h>

void testPrintGenotypes(struct HashMapIntToIntArr *map) {
  struct NodeIntToIntArr *itr;

  for (int i = 0; i < 13; ++i) {
    itr = map->list[i];
    printf("List %d:\n", i);
    while (itr) {
      printf("Key: %d, Value: (%d, %d, %d)\n", itr->key, itr->val[0],
             itr->val[1], itr->val[2]);
      itr = itr->next;
    }
  }
}
