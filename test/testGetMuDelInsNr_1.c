//
// Created by Senbai Kang on 11/25/19.
//

#include "testFunctions.h"
#include <unity.c>
#include "testGetMuDelInsNr.h"

void testPa1Wrapper(void);
void testPa2Wrapper(void);
void testPa3Wrapper(void);
void testPa4Wrapper(void);
void testPa5Wrapper(void);
void testPa6Wrapper(void);
void testPa7Wrapper(void);
void testPa8Wrapper(void);
void testPa9Wrapper(void);
void testPa10Wrapper(void);
void testPa11Wrapper(void);
void testPa12Wrapper(void);
void testPa13Wrapper(void);
void testPa14Wrapper(void);
void testPa15Wrapper(void);
void testPa16Wrapper(void);

int main(void) {
  map = initHashMapIntToIntArrMDI();

  AllocateData();

  UNITY_BEGIN();

  RUN_TEST(testPa1Wrapper);
  RUN_TEST(testPa2Wrapper);
  //  RUN_TEST(testPa3Wrapper);
  //  RUN_TEST(testPa4Wrapper);
  //  RUN_TEST(testPa5Wrapper);
  //  RUN_TEST(testPa6Wrapper);
  //  RUN_TEST(testPa7Wrapper);
  //  RUN_TEST(testPa8Wrapper);
  RUN_TEST(testPa9Wrapper);
  //  RUN_TEST(testPa10Wrapper);
  //  RUN_TEST(testPa11Wrapper);
  //  RUN_TEST(testPa12Wrapper);
  //  RUN_TEST(testPa13Wrapper);
  //  RUN_TEST(testPa14Wrapper);
  //  RUN_TEST(testPa15Wrapper);
  //  RUN_TEST(testPa16Wrapper);

  FreeData();
  FreeGenoArr(65);
  FreeIntToIntArrHashMapMDI();

  return UNITY_END();
}

void testPa1Wrapper(void) {
  SetDataStrMDI(dataPa, 0);

  // tested in the order of number of:
  // mutations, back mutations, potential back mutations, deletions, and
  // insertions
  SetDataStrMDI(dataCh, 1);
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

  SetDataStrMDI(dataCh, 2);
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

  SetDataStrMDI(dataCh, 3);
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

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

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

  SetDataStrMDI(dataCh, 5);
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

  SetDataStrMDI(dataCh, 6);
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

  SetDataStrMDI(dataCh, 7);
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

  SetDataStrMDI(dataCh, 8);
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

  SetDataStrMDI(dataCh, 9);
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

  SetDataStrMDI(dataCh, 10);
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

  SetDataStrMDI(dataCh, 11);
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

  SetDataStrMDI(dataCh, 12);
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

  SetDataStrMDI(dataCh, 13);
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

  SetDataStrMDI(dataCh, 14);
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

  SetDataStrMDI(dataCh, 15);
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

  SetDataStrMDI(dataCh, 16);
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

  SetDataStrMDI(dataCh, 17);
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

  SetDataStrMDI(dataCh, 18);
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

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);

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

  SetDataStrMDI(dataCh, 20);
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

  SetDataStrMDI(dataCh, 21);
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

  SetDataStrMDI(dataCh, 22);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 23);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 29);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 30);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 35);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 36);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
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

  SetDataStrMDI(dataCh, 41);
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

  SetDataStrMDI(dataCh, 42);
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

  SetDataStrMDI(dataCh, 43);
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

  SetDataStrMDI(dataCh, 44);
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

  SetDataStrMDI(dataCh, 45);
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

  SetDataStrMDI(dataCh, 46);
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

  SetDataStrMDI(dataCh, 47);
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

  SetDataStrMDI(dataCh, 48);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 55);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 56);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
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

  SetDataStrMDI(dataCh, 61);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa2Wrapper(void) {
  SetDataStrMDI(dataPa, 1);

  // tested in the order of number of:
  // mutations, back mutations, potential back mutations, deletions, and
  // insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
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

  SetDataStrMDI(dataCh, 3);
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

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
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

  SetDataStrMDI(dataCh, 8);
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

  SetDataStrMDI(dataCh, 9);
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

  SetDataStrMDI(dataCh, 10);
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

  SetDataStrMDI(dataCh, 11);
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

  SetDataStrMDI(dataCh, 12);
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

  SetDataStrMDI(dataCh, 13);
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

  SetDataStrMDI(dataCh, 14);
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

  SetDataStrMDI(dataCh, 15);
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

  SetDataStrMDI(dataCh, 16);
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

  SetDataStrMDI(dataCh, 17);
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

  SetDataStrMDI(dataCh, 18);
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

  SetDataStrMDI(dataCh, 19);
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

  SetDataStrMDI(dataCh, 20);
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

  SetDataStrMDI(dataCh, 21);
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

  SetDataStrMDI(dataCh, 22);
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

  SetDataStrMDI(dataCh, 23);
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

  SetDataStrMDI(dataCh, 24);
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

  SetDataStrMDI(dataCh, 25);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 26);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
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

  SetDataStrMDI(dataCh, 29);
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

  SetDataStrMDI(dataCh, 30);
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

  SetDataStrMDI(dataCh, 31);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 32);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
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

  SetDataStrMDI(dataCh, 35);
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

  SetDataStrMDI(dataCh, 36);
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

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
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

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
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

  SetDataStrMDI(dataCh, 43);
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

  SetDataStrMDI(dataCh, 44);
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

  SetDataStrMDI(dataCh, 45);
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

  SetDataStrMDI(dataCh, 46);
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

  SetDataStrMDI(dataCh, 47);
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

  SetDataStrMDI(dataCh, 48);
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

  SetDataStrMDI(dataCh, 49);
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

  SetDataStrMDI(dataCh, 50);
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

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
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

  SetDataStrMDI(dataCh, 55);
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

  SetDataStrMDI(dataCh, 56);
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

  SetDataStrMDI(dataCh, 57);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 58);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
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

  SetDataStrMDI(dataCh, 61);
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

  SetDataStrMDI(dataCh, 62);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa3Wrapper(void) {
  SetDataStrMDI(dataPa, 2);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa4Wrapper(void) {
  SetDataStrMDI(dataPa, 3);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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

void testPa5Wrapper(void) {
  SetDataStrMDI(dataPa, 4);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa6Wrapper(void) {
  SetDataStrMDI(dataPa, 5);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa7Wrapper(void) {
  SetDataStrMDI(dataPa, 6);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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

void testPa8Wrapper(void) {
  SetDataStrMDI(dataPa, 7);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 9);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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

  SetDataStrMDI(dataCh, 34);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa9Wrapper(void) {
  SetDataStrMDI(dataPa, 8);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 6);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
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

  SetDataStrMDI(dataCh, 9);
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

  SetDataStrMDI(dataCh, 10);
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

  SetDataStrMDI(dataCh, 11);
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

  SetDataStrMDI(dataCh, 12);
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

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 17);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 18);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 19);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 20);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 21);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 22);
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

  SetDataStrMDI(dataCh, 23);
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

  SetDataStrMDI(dataCh, 24);
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

  SetDataStrMDI(dataCh, 25);
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

  SetDataStrMDI(dataCh, 26);
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

  SetDataStrMDI(dataCh, 27);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 28);
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

  SetDataStrMDI(dataCh, 29);
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

  SetDataStrMDI(dataCh, 30);
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

  SetDataStrMDI(dataCh, 31);
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

  SetDataStrMDI(dataCh, 32);
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

  SetDataStrMDI(dataCh, 33);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 34);
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

  SetDataStrMDI(dataCh, 35);
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

  SetDataStrMDI(dataCh, 36);
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

  SetDataStrMDI(dataCh, 37);
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

  SetDataStrMDI(dataCh, 38);
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

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
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

  SetDataStrMDI(dataCh, 43);
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

  SetDataStrMDI(dataCh, 44);
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

  SetDataStrMDI(dataCh, 45);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 46);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 47);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 48);
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

  SetDataStrMDI(dataCh, 49);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 50);
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

  SetDataStrMDI(dataCh, 51);
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

  SetDataStrMDI(dataCh, 52);
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

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 54);
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

  SetDataStrMDI(dataCh, 55);
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

  SetDataStrMDI(dataCh, 56);
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

  SetDataStrMDI(dataCh, 57);
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

  SetDataStrMDI(dataCh, 58);
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

  SetDataStrMDI(dataCh, 59);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
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

  SetDataStrMDI(dataCh, 61);
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

  SetDataStrMDI(dataCh, 62);
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

  SetDataStrMDI(dataCh, 63);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(1, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(2, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  GetMuDelInsNrCore(3, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[1]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[2]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa10Wrapper(void) {
  SetDataStrMDI(dataPa, 9);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 7);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 8);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 43);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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

void testPa11Wrapper(void) {
  SetDataStrMDI(dataPa, 10);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 37);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 38);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 51);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa12Wrapper(void) {
  SetDataStrMDI(dataPa, 11);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 10);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 44);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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

void testPa13Wrapper(void) {
  SetDataStrMDI(dataPa, 12);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 11);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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

void testPa14Wrapper(void) {
  SetDataStrMDI(dataPa, 13);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 16);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 24);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 39);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 52);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa15Wrapper(void) {
  SetDataStrMDI(dataPa, 14);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 1);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 2);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 3);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 15);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 53);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 64);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(3, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);
}

void testPa16Wrapper(void) {
  SetDataStrMDI(dataPa, 15);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(dataCh, 0);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 4);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 5);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 13);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 14);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(2, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 40);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 41);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 42);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[0]);
  TEST_ASSERT_EQUAL_INT(1, evolutionaryEvents[3]);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[4]);

  SetDataStrMDI(dataCh, 60);
  GetMuDelInsNrCore(0, 3);
  TEST_ASSERT_EQUAL_INT(0, evolutionaryEvents[0]);
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
