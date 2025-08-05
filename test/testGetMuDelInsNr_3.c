//
// Created by Senbai Kang on 11/25/19.
//

#include "testFunctions.h"
#include <unity.c>
#include "testGetMuDelInsNr.h"

void testPa33Wrapper(void);
void testPa34Wrapper(void);
void testPa35Wrapper(void);
void testPa36Wrapper(void);
void testPa37Wrapper(void);
void testPa38Wrapper(void);
void testPa39Wrapper(void);
void testPa40Wrapper(void);
void testPa41Wrapper(void);
void testPa42Wrapper(void);
void testPa43Wrapper(void);
void testPa44Wrapper(void);
void testPa45Wrapper(void);
void testPa46Wrapper(void);
void testPa47Wrapper(void);
void testPa48Wrapper(void);

int main(void) {
  map = initHashMapIntToIntArrMDI();

  AllocateData();

  UNITY_BEGIN();

  //  RUN_TEST(testPa33Wrapper);
  //  RUN_TEST(testPa34Wrapper);
  //  RUN_TEST(testPa35Wrapper);
  //  RUN_TEST(testPa36Wrapper);
  //  RUN_TEST(testPa37Wrapper);
  //  RUN_TEST(testPa38Wrapper);
  //  RUN_TEST(testPa39Wrapper);
  //  RUN_TEST(testPa40Wrapper);
  RUN_TEST(testPa41Wrapper);
  //  RUN_TEST(testPa42Wrapper);
  //  RUN_TEST(testPa43Wrapper);
  //  RUN_TEST(testPa44Wrapper);
  //  RUN_TEST(testPa45Wrapper);
  //  RUN_TEST(testPa46Wrapper);
  //  RUN_TEST(testPa47Wrapper);
  //  RUN_TEST(testPa48Wrapper);

  FreeData();
  FreeGenoArr(65);
  FreeIntToIntArrHashMapMDI();

  return UNITY_END();
}

void testPa33Wrapper(void) {
  SetDataStrMDI(dataPa, 32);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa34Wrapper(void) {
  SetDataStrMDI(dataPa, 33);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa35Wrapper(void) {
  SetDataStrMDI(dataPa, 34);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa36Wrapper(void) {
  SetDataStrMDI(dataPa, 35);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa37Wrapper(void) {
  SetDataStrMDI(dataPa, 36);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa38Wrapper(void) {
  SetDataStrMDI(dataPa, 37);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa39Wrapper(void) {
  SetDataStrMDI(dataPa, 38);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa40Wrapper(void) {
  SetDataStrMDI(dataPa, 39);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa41Wrapper(void) {
  SetDataStrMDI(dataPa, 40);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa42Wrapper(void) {
  SetDataStrMDI(dataPa, 41);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa43Wrapper(void) {
  SetDataStrMDI(dataPa, 42);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa44Wrapper(void) {
  SetDataStrMDI(dataPa, 43);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa45Wrapper(void) {
  SetDataStrMDI(dataPa, 44);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa46Wrapper(void) {
  SetDataStrMDI(dataPa, 45);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa47Wrapper(void) {
  SetDataStrMDI(dataPa, 46);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa48Wrapper(void) {
  SetDataStrMDI(dataPa, 47);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 12);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}
