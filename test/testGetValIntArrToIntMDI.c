//
// Created by Senbai Kang on 12/3/19.
//

#include "testFunctions.h"
#include <unity.c>

static struct HashMapIntArrToInt *map;

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

int main(void) {
  map = initHashMapIntArrToIntMDI();

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

  return UNITY_END();
}

void testCase0(void) {
  int key_1[3] = {-3, 0, 0};
  TEST_ASSERT_EQUAL(41, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 0};
  TEST_ASSERT_EQUAL(40, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 0, -3};
  TEST_ASSERT_EQUAL(40, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 0};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase1(void) {
  int key_1[3] = {-3, 0, 1};
  TEST_ASSERT_EQUAL(45, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 1};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 0, -3};
  TEST_ASSERT_EQUAL(40, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase2(void) {
  int key_1[3] = {-3, 0, 2};
  TEST_ASSERT_EQUAL(46, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 2};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 0, -3};
  TEST_ASSERT_EQUAL(40, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase3(void) {
  int key_1[3] = {-3, 0, 3};
  TEST_ASSERT_EQUAL(47, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 3};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 0, -3};
  TEST_ASSERT_EQUAL(40, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase4(void) {
  int key_1[3] = {-3, 0, 0};
  TEST_ASSERT_EQUAL(41, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 0};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 0, -3};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 0};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase5(void) {
  int key_1[3] = {-3, 0, 0};
  TEST_ASSERT_EQUAL(41, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 0};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 0, -3};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 0};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase6(void) {
  int key_1[3] = {-3, 0, 0};
  TEST_ASSERT_EQUAL(41, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 0};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 0, -3};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 0};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase7(void) {
  int key_1[3] = {-3, 1, 1};
  TEST_ASSERT_EQUAL(54, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 1};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 1, -3};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase8(void) {
  int key_1[3] = {-3, 1, 2};
  TEST_ASSERT_EQUAL(55, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 2};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 1, -3};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase9(void) {
  int key_1[3] = {-3, 1, 3};
  TEST_ASSERT_EQUAL(56, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 3};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 1, -3};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase10(void) {
  int key_1[3] = {-3, 2, 2};
  TEST_ASSERT_EQUAL(57, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 2};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 2, -3};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase11(void) {
  int key_1[3] = {-3, 2, 3};
  TEST_ASSERT_EQUAL(58, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 3};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 2, -3};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase12(void) {
  int key_1[3] = {-3, 3, 3};
  TEST_ASSERT_EQUAL(59, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, 3};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {0, 3, -3};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {0, -3, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase13(void) {
  int key_1[3] = {-3, 0, 1};
  TEST_ASSERT_EQUAL(45, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 1};
  TEST_ASSERT_EQUAL(48, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 0, -3};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase14(void) {
  int key_1[3] = {-3, 0, 2};
  TEST_ASSERT_EQUAL(46, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 2};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 0, -3};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase15(void) {
  int key_1[3] = {-3, 0, 3};
  TEST_ASSERT_EQUAL(47, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 3};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 0, -3};
  TEST_ASSERT_EQUAL(42, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase16(void) {
  int key_1[3] = {-3, 0, 1};
  TEST_ASSERT_EQUAL(45, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 1};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 0, -3};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase17(void) {
  int key_1[3] = {-3, 0, 2};
  TEST_ASSERT_EQUAL(46, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 2};
  TEST_ASSERT_EQUAL(51, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 0, -3};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase18(void) {
  int key_1[3] = {-3, 0, 3};
  TEST_ASSERT_EQUAL(47, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 3};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 0, -3};
  TEST_ASSERT_EQUAL(43, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase19(void) {
  int key_1[3] = {-3, 0, 1};
  TEST_ASSERT_EQUAL(45, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 1};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 0, -3};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase20(void) {
  int key_1[3] = {-3, 0, 2};
  TEST_ASSERT_EQUAL(46, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 2};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 0, -3};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase21(void) {
  int key_1[3] = {-3, 0, 3};
  TEST_ASSERT_EQUAL(47, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 3};
  TEST_ASSERT_EQUAL(53, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 0, -3};
  TEST_ASSERT_EQUAL(44, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase22(void) {
  int key_1[3] = {-3, 1, 1};
  TEST_ASSERT_EQUAL(54, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 1};
  TEST_ASSERT_EQUAL(48, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 1, -3};
  TEST_ASSERT_EQUAL(48, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase23(void) {
  int key_1[3] = {-3, 1, 2};
  TEST_ASSERT_EQUAL(55, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 2};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 1, -3};
  TEST_ASSERT_EQUAL(48, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase24(void) {
  int key_1[3] = {-3, 1, 3};
  TEST_ASSERT_EQUAL(56, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 3};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 1, -3};
  TEST_ASSERT_EQUAL(48, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase25(void) {
  int key_1[3] = {-3, 2, 2};
  TEST_ASSERT_EQUAL(57, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 2};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 2, -3};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase26(void) {
  int key_1[3] = {-3, 2, 3};
  TEST_ASSERT_EQUAL(58, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 3};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 2, -3};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase27(void) {
  int key_1[3] = {-3, 3, 3};
  TEST_ASSERT_EQUAL(59, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, 3};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {1, 3, -3};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {1, -3, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase28(void) {
  int key_1[3] = {-3, 1, 1};
  TEST_ASSERT_EQUAL(54, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 1};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 1, -3};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase29(void) {
  int key_1[3] = {-3, 1, 2};
  TEST_ASSERT_EQUAL(55, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 2};
  TEST_ASSERT_EQUAL(51, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 1, -3};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase30(void) {
  int key_1[3] = {-3, 1, 3};
  TEST_ASSERT_EQUAL(56, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 3};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 1, -3};
  TEST_ASSERT_EQUAL(49, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase31(void) {
  int key_1[3] = {-3, 2, 2};
  TEST_ASSERT_EQUAL(57, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 2};
  TEST_ASSERT_EQUAL(51, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 2, -3};
  TEST_ASSERT_EQUAL(51, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase32(void) {
  int key_1[3] = {-3, 2, 3};
  TEST_ASSERT_EQUAL(58, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 3};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 2, -3};
  TEST_ASSERT_EQUAL(51, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase33(void) {
  int key_1[3] = {-3, 3, 3};
  TEST_ASSERT_EQUAL(59, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, 3};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {2, 3, -3};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {2, -3, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase34(void) {
  int key_1[3] = {-3, 1, 1};
  TEST_ASSERT_EQUAL(54, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 1};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 1, -3};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase35(void) {
  int key_1[3] = {-3, 1, 2};
  TEST_ASSERT_EQUAL(55, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 2};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 1, -3};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase36(void) {
  int key_1[3] = {-3, 1, 3};
  TEST_ASSERT_EQUAL(56, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 3};
  TEST_ASSERT_EQUAL(53, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 1, -3};
  TEST_ASSERT_EQUAL(50, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase37(void) {
  int key_1[3] = {-3, 2, 2};
  TEST_ASSERT_EQUAL(57, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 2};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 2, -3};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase38(void) {
  int key_1[3] = {-3, 2, 3};
  TEST_ASSERT_EQUAL(58, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 3};
  TEST_ASSERT_EQUAL(53, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 2, -3};
  TEST_ASSERT_EQUAL(52, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase39(void) {
  int key_1[3] = {-3, 3, 3};
  TEST_ASSERT_EQUAL(59, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, 3};
  TEST_ASSERT_EQUAL(53, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {3, 3, -3};
  TEST_ASSERT_EQUAL(53, GetValIntArrToInt(map, key_3, 0));

  int key_4[3] = {-3, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_4, 0));

  int key_5[3] = {-3, 3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_5, 0));

  int key_6[3] = {3, -3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_6, 0));

  int key_7[3] = {-3, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_7, 0));
}

void testCase40(void) {
  int key_1[3] = {-3, 0, -1};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, -1};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase41(void) {
  int key_1[3] = {-1, -3, 0};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase42(void) {
  int key_1[3] = {-3, 1, -1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, -1};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase43(void) {
  int key_1[3] = {-3, 2, -1};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, -1};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase44(void) {
  int key_1[3] = {-3, 3, -1};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {0, -3, -1};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase45(void) {
  int key_1[3] = {-1, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase46(void) {
  int key_1[3] = {-1, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase47(void) {
  int key_1[3] = {-1, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 0, -3};
  TEST_ASSERT_EQUAL(60, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase48(void) {
  int key_1[3] = {-3, 1, -1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, -1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase49(void) {
  int key_1[3] = {-3, 2, -1};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, -1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase50(void) {
  int key_1[3] = {-3, 3, -1};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {1, -3, -1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase51(void) {
  int key_1[3] = {-3, 2, -1};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, -1};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase52(void) {
  int key_1[3] = {-3, 3, -1};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {2, -3, -1};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase53(void) {
  int key_1[3] = {-3, 3, -1};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {3, -3, -1};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-3, -3, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase54(void) {
  int key_1[3] = {-1, -3, 1};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase55(void) {
  int key_1[3] = {-1, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase56(void) {
  int key_1[3] = {-1, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 1, -3};
  TEST_ASSERT_EQUAL(61, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase57(void) {
  int key_1[3] = {-1, -3, 2};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase58(void) {
  int key_1[3] = {-1, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 2, -3};
  TEST_ASSERT_EQUAL(62, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase59(void) {
  int key_1[3] = {-1, -3, 3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_1, 0));

  int key_2[3] = {-1, 3, -3};
  TEST_ASSERT_EQUAL(63, GetValIntArrToInt(map, key_2, 0));

  int key_3[3] = {-1, -3, -3};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_3, 0));
}

void testCase60(void) {
  int key_1[3] = {-3, -2, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_1, 0));
}

void testCase61(void) {
  int key_1[3] = {-3, -2, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_1, 0));
}

void testCase62(void) {
  int key_1[3] = {-3, -2, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_1, 0));
}

void testCase63(void) {
  int key_1[3] = {-3, -2, -1};
  TEST_ASSERT_EQUAL(64, GetValIntArrToInt(map, key_1, 0));
}
