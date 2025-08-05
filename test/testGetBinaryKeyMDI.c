//
// Created by Ethan Kang on 12/3/19.
//

#include "testFunctions.h"
#include <unity.c>

static struct HashMapIntToIntArr *map;

void testCase0(void);
void testCase1(void);
void testCase2(void);
void testCase3(void);
void testCase4(void);
void testCase5(void);
void testCase6(void);
void testCase7(void);
void testCase8(void);
void testCase9(void);
void testCase10(void);
void testCase11(void);
void testCase12(void);
void testCase13(void);
void testCase14(void);
void testCase15(void);
void testCase16(void);
void testCase17(void);
void testCase18(void);
void testCase19(void);
void testCase20(void);
void testCase21(void);
void testCase22(void);
void testCase23(void);
void testCase24(void);
void testCase25(void);
void testCase26(void);
void testCase27(void);
void testCase28(void);
void testCase29(void);
void testCase30(void);
void testCase31(void);
void testCase32(void);
void testCase33(void);
void testCase34(void);
void testCase35(void);
void testCase36(void);
void testCase37(void);
void testCase38(void);
void testCase39(void);
void testCase40(void);
void testCase41(void);
void testCase42(void);
void testCase43(void);
void testCase44(void);
void testCase45(void);
void testCase46(void);
void testCase47(void);
void testCase48(void);
void testCase49(void);
void testCase50(void);
void testCase51(void);
void testCase52(void);
void testCase53(void);
void testCase54(void);
void testCase55(void);
void testCase56(void);
void testCase57(void);
void testCase58(void);
void testCase59(void);
void testCase60(void);
void testCase61(void);
void testCase62(void);
void testCase63(void);
void testCase64(void);

int main(void) {
  map = initHashMapIntToIntArrMDI();

  UNITY_BEGIN();

  RUN_TEST(testCase0);
  RUN_TEST(testCase1);
  RUN_TEST(testCase2);
  RUN_TEST(testCase3);
  RUN_TEST(testCase4);
  RUN_TEST(testCase5);
  RUN_TEST(testCase6);
  RUN_TEST(testCase7);
  RUN_TEST(testCase8);
  RUN_TEST(testCase9);
  RUN_TEST(testCase10);
  RUN_TEST(testCase11);
  RUN_TEST(testCase12);
  RUN_TEST(testCase13);
  RUN_TEST(testCase14);
  RUN_TEST(testCase15);
  RUN_TEST(testCase16);
  RUN_TEST(testCase17);
  RUN_TEST(testCase18);
  RUN_TEST(testCase19);
  RUN_TEST(testCase20);
  RUN_TEST(testCase21);
  RUN_TEST(testCase22);
  RUN_TEST(testCase23);
  RUN_TEST(testCase24);
  RUN_TEST(testCase25);
  RUN_TEST(testCase26);
  RUN_TEST(testCase27);
  RUN_TEST(testCase28);
  RUN_TEST(testCase29);
  RUN_TEST(testCase30);
  RUN_TEST(testCase31);
  RUN_TEST(testCase32);
  RUN_TEST(testCase33);
  RUN_TEST(testCase34);
  RUN_TEST(testCase35);
  RUN_TEST(testCase36);
  RUN_TEST(testCase37);
  RUN_TEST(testCase38);
  RUN_TEST(testCase39);
  RUN_TEST(testCase40);
  RUN_TEST(testCase41);
  RUN_TEST(testCase42);
  RUN_TEST(testCase43);
  RUN_TEST(testCase44);
  RUN_TEST(testCase45);
  RUN_TEST(testCase46);
  RUN_TEST(testCase47);
  RUN_TEST(testCase48);
  RUN_TEST(testCase49);
  RUN_TEST(testCase50);
  RUN_TEST(testCase51);
  RUN_TEST(testCase52);
  RUN_TEST(testCase53);
  RUN_TEST(testCase54);
  RUN_TEST(testCase55);
  RUN_TEST(testCase56);
  RUN_TEST(testCase57);
  RUN_TEST(testCase58);
  RUN_TEST(testCase59);
  RUN_TEST(testCase60);
  RUN_TEST(testCase61);
  RUN_TEST(testCase62);
  RUN_TEST(testCase63);
  RUN_TEST(testCase64);

  FreeGenoArr(65);
  FreeIntToIntArrHashMapMDI();

  return UNITY_END();
}

void testCase0(void) {
  TEST_ASSERT_EQUAL(0, GetBinaryKey(GetValIntToIntArr(map, 0)));
}

void testCase1(void) {
  TEST_ASSERT_EQUAL(1, GetBinaryKey(GetValIntToIntArr(map, 1)));
}

void testCase2(void) {
  TEST_ASSERT_EQUAL(10, GetBinaryKey(GetValIntToIntArr(map, 2)));
}

void testCase3(void) {
  TEST_ASSERT_EQUAL(11, GetBinaryKey(GetValIntToIntArr(map, 3)));
}

void testCase4(void) {
  TEST_ASSERT_EQUAL(1000000, GetBinaryKey(GetValIntToIntArr(map, 4)));
}

void testCase5(void) {
  TEST_ASSERT_EQUAL(10000000, GetBinaryKey(GetValIntToIntArr(map, 5)));
}

void testCase6(void) {
  TEST_ASSERT_EQUAL(11000000, GetBinaryKey(GetValIntToIntArr(map, 6)));
}

void testCase7(void) {
  TEST_ASSERT_EQUAL(1001, GetBinaryKey(GetValIntToIntArr(map, 7)));
}

void testCase8(void) {
  TEST_ASSERT_EQUAL(1010, GetBinaryKey(GetValIntToIntArr(map, 8)));
}

void testCase9(void) {
  TEST_ASSERT_EQUAL(1011, GetBinaryKey(GetValIntToIntArr(map, 9)));
}

void testCase10(void) {
  TEST_ASSERT_EQUAL(10010, GetBinaryKey(GetValIntToIntArr(map, 10)));
}

void testCase11(void) {
  TEST_ASSERT_EQUAL(10011, GetBinaryKey(GetValIntToIntArr(map, 11)));
}

void testCase12(void) {
  TEST_ASSERT_EQUAL(11011, GetBinaryKey(GetValIntToIntArr(map, 12)));
}

void testCase13(void) {
  TEST_ASSERT_EQUAL(1000001, GetBinaryKey(GetValIntToIntArr(map, 13)));
}

void testCase14(void) {
  TEST_ASSERT_EQUAL(1000010, GetBinaryKey(GetValIntToIntArr(map, 14)));
}

void testCase15(void) {
  TEST_ASSERT_EQUAL(1000011, GetBinaryKey(GetValIntToIntArr(map, 15)));
}

void testCase16(void) {
  TEST_ASSERT_EQUAL(10000001, GetBinaryKey(GetValIntToIntArr(map, 16)));
}

void testCase17(void) {
  TEST_ASSERT_EQUAL(10000010, GetBinaryKey(GetValIntToIntArr(map, 17)));
}

void testCase18(void) {
  TEST_ASSERT_EQUAL(10000011, GetBinaryKey(GetValIntToIntArr(map, 18)));
}

void testCase19(void) {
  TEST_ASSERT_EQUAL(11000001, GetBinaryKey(GetValIntToIntArr(map, 19)));
}

void testCase20(void) {
  TEST_ASSERT_EQUAL(11000010, GetBinaryKey(GetValIntToIntArr(map, 20)));
}

void testCase21(void) {
  TEST_ASSERT_EQUAL(11000011, GetBinaryKey(GetValIntToIntArr(map, 21)));
}

void testCase22(void) {
  TEST_ASSERT_EQUAL(1001001, GetBinaryKey(GetValIntToIntArr(map, 22)));
}

void testCase23(void) {
  TEST_ASSERT_EQUAL(1001010, GetBinaryKey(GetValIntToIntArr(map, 23)));
}

void testCase24(void) {
  TEST_ASSERT_EQUAL(1001011, GetBinaryKey(GetValIntToIntArr(map, 24)));
}

void testCase25(void) {
  TEST_ASSERT_EQUAL(1010010, GetBinaryKey(GetValIntToIntArr(map, 25)));
}

void testCase26(void) {
  TEST_ASSERT_EQUAL(1010011, GetBinaryKey(GetValIntToIntArr(map, 26)));
}

void testCase27(void) {
  TEST_ASSERT_EQUAL(1011011, GetBinaryKey(GetValIntToIntArr(map, 27)));
}

void testCase28(void) {
  TEST_ASSERT_EQUAL(10001001, GetBinaryKey(GetValIntToIntArr(map, 28)));
}

void testCase29(void) {
  TEST_ASSERT_EQUAL(10001010, GetBinaryKey(GetValIntToIntArr(map, 29)));
}

void testCase30(void) {
  TEST_ASSERT_EQUAL(10001011, GetBinaryKey(GetValIntToIntArr(map, 30)));
}

void testCase31(void) {
  TEST_ASSERT_EQUAL(10010010, GetBinaryKey(GetValIntToIntArr(map, 31)));
}

void testCase32(void) {
  TEST_ASSERT_EQUAL(10010011, GetBinaryKey(GetValIntToIntArr(map, 32)));
}

void testCase33(void) {
  TEST_ASSERT_EQUAL(10011011, GetBinaryKey(GetValIntToIntArr(map, 33)));
}

void testCase34(void) {
  TEST_ASSERT_EQUAL(11001001, GetBinaryKey(GetValIntToIntArr(map, 34)));
}

void testCase35(void) {
  TEST_ASSERT_EQUAL(11001010, GetBinaryKey(GetValIntToIntArr(map, 35)));
}

void testCase36(void) {
  TEST_ASSERT_EQUAL(11001011, GetBinaryKey(GetValIntToIntArr(map, 36)));
}

void testCase37(void) {
  TEST_ASSERT_EQUAL(11010010, GetBinaryKey(GetValIntToIntArr(map, 37)));
}

void testCase38(void) {
  TEST_ASSERT_EQUAL(11010011, GetBinaryKey(GetValIntToIntArr(map, 38)));
}

void testCase39(void) {
  TEST_ASSERT_EQUAL(11011011, GetBinaryKey(GetValIntToIntArr(map, 39)));
}

void testCase40(void) {
  TEST_ASSERT_EQUAL(100, GetBinaryKey(GetValIntToIntArr(map, 40)));
}

void testCase41(void) {
  TEST_ASSERT_EQUAL(101000000, GetBinaryKey(GetValIntToIntArr(map, 41)));
}

void testCase42(void) {
  TEST_ASSERT_EQUAL(1100, GetBinaryKey(GetValIntToIntArr(map, 42)));
}

void testCase43(void) {
  TEST_ASSERT_EQUAL(10100, GetBinaryKey(GetValIntToIntArr(map, 43)));
}

void testCase44(void) {
  TEST_ASSERT_EQUAL(11100, GetBinaryKey(GetValIntToIntArr(map, 44)));
}

void testCase45(void) {
  TEST_ASSERT_EQUAL(101000001, GetBinaryKey(GetValIntToIntArr(map, 45)));
}

void testCase46(void) {
  TEST_ASSERT_EQUAL(101000010, GetBinaryKey(GetValIntToIntArr(map, 46)));
}

void testCase47(void) {
  TEST_ASSERT_EQUAL(101000011, GetBinaryKey(GetValIntToIntArr(map, 47)));
}

void testCase48(void) {
  TEST_ASSERT_EQUAL(1001100, GetBinaryKey(GetValIntToIntArr(map, 48)));
}

void testCase49(void) {
  TEST_ASSERT_EQUAL(1010100, GetBinaryKey(GetValIntToIntArr(map, 49)));
}

void testCase50(void) {
  TEST_ASSERT_EQUAL(1011100, GetBinaryKey(GetValIntToIntArr(map, 50)));
}

void testCase51(void) {
  TEST_ASSERT_EQUAL(10010100, GetBinaryKey(GetValIntToIntArr(map, 51)));
}

void testCase52(void) {
  TEST_ASSERT_EQUAL(10011100, GetBinaryKey(GetValIntToIntArr(map, 52)));
}

void testCase53(void) {
  TEST_ASSERT_EQUAL(11011100, GetBinaryKey(GetValIntToIntArr(map, 53)));
}

void testCase54(void) {
  TEST_ASSERT_EQUAL(101001001, GetBinaryKey(GetValIntToIntArr(map, 54)));
}

void testCase55(void) {
  TEST_ASSERT_EQUAL(101001010, GetBinaryKey(GetValIntToIntArr(map, 55)));
}

void testCase56(void) {
  TEST_ASSERT_EQUAL(101001011, GetBinaryKey(GetValIntToIntArr(map, 56)));
}

void testCase57(void) {
  TEST_ASSERT_EQUAL(101010010, GetBinaryKey(GetValIntToIntArr(map, 57)));
}

void testCase58(void) {
  TEST_ASSERT_EQUAL(101010011, GetBinaryKey(GetValIntToIntArr(map, 58)));
}

void testCase59(void) {
  TEST_ASSERT_EQUAL(101011011, GetBinaryKey(GetValIntToIntArr(map, 59)));
}

void testCase60(void) {
  TEST_ASSERT_EQUAL(101100, GetBinaryKey(GetValIntToIntArr(map, 60)));
}

void testCase61(void) {
  TEST_ASSERT_EQUAL(1101100, GetBinaryKey(GetValIntToIntArr(map, 61)));
}

void testCase62(void) {
  TEST_ASSERT_EQUAL(10101100, GetBinaryKey(GetValIntToIntArr(map, 62)));
}

void testCase63(void) {
  TEST_ASSERT_EQUAL(11101100, GetBinaryKey(GetValIntToIntArr(map, 63)));
}

void testCase64(void) {
  TEST_ASSERT_EQUAL(101101100, GetBinaryKey(GetValIntToIntArr(map, 64)));
}
