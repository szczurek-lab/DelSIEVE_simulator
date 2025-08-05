//
// Created by Senbai Kang on 11/24/19.
//

#include "testFunctions.h"
#include <unity.c>

void testPaEqual_1(struct HashMapIntToIntArr *map);
void testPaEqual_2(struct HashMapIntToIntArr *map);
void testPaEqual_3(struct HashMapIntToIntArr *map);
void testPaEqual_4(struct HashMapIntToIntArr *map);
void testPaEqual_5(struct HashMapIntToIntArr *map);
void testPaEqual_6(struct HashMapIntToIntArr *map);
void testPaEqual_7(struct HashMapIntToIntArr *map);
void testPaEqual_8(struct HashMapIntToIntArr *map);
void testPaEqual_9(struct HashMapIntToIntArr *map);
void testPaEqual_10(struct HashMapIntToIntArr *map);
void testPaEqual_11(struct HashMapIntToIntArr *map);
void testPaEqual_12(struct HashMapIntToIntArr *map);
void testPaEqual_13(struct HashMapIntToIntArr *map);
void testPaEqual_14(struct HashMapIntToIntArr *map);
void testPaEqual_15(struct HashMapIntToIntArr *map);
void testPaEqual_16(struct HashMapIntToIntArr *map);
void testPaEqual_17(struct HashMapIntToIntArr *map);
void testPaEqual_18(struct HashMapIntToIntArr *map);
void testPaEqual_19(struct HashMapIntToIntArr *map);
void testPaEqual_20(struct HashMapIntToIntArr *map);
void testPaEqual_21(struct HashMapIntToIntArr *map);
void testPaEqual_22(struct HashMapIntToIntArr *map);
void testPaEqual_23(struct HashMapIntToIntArr *map);
void testPaEqual_24(struct HashMapIntToIntArr *map);
void testPaEqual_25(struct HashMapIntToIntArr *map);
void testPaEqual_26(struct HashMapIntToIntArr *map);
void testPaEqual_27(struct HashMapIntToIntArr *map);
void testPaEqual_28(struct HashMapIntToIntArr *map);
void testPaEqual_29(struct HashMapIntToIntArr *map);
void testPaEqual_30(struct HashMapIntToIntArr *map);
void testPaEqual_31(struct HashMapIntToIntArr *map);
void testPaEqual_32(struct HashMapIntToIntArr *map);
void testPaEqual_33(struct HashMapIntToIntArr *map);
void testPaEqual_34(struct HashMapIntToIntArr *map);
void testPaEqual_35(struct HashMapIntToIntArr *map);
void testPaEqual_36(struct HashMapIntToIntArr *map);
void testPaEqual_37(struct HashMapIntToIntArr *map);
void testPaEqual_38(struct HashMapIntToIntArr *map);
void testPaEqual_39(struct HashMapIntToIntArr *map);
void testPaEqual_40(struct HashMapIntToIntArr *map);
void testPaEqual_41(struct HashMapIntToIntArr *map);
void testPaEqual_42(struct HashMapIntToIntArr *map);
void testPaEqual_43(struct HashMapIntToIntArr *map);
void testPaEqual_44(struct HashMapIntToIntArr *map);
void testPaEqual_45(struct HashMapIntToIntArr *map);
void testPaEqual_46(struct HashMapIntToIntArr *map);
void testPaEqual_47(struct HashMapIntToIntArr *map);
void testPaEqual_48(struct HashMapIntToIntArr *map);
void testPaEqual_49(struct HashMapIntToIntArr *map);
void testPaEqual_50(struct HashMapIntToIntArr *map);
void testPaEqual_51(struct HashMapIntToIntArr *map);
void testPaEqual_52(struct HashMapIntToIntArr *map);
void testPaEqual_53(struct HashMapIntToIntArr *map);
void testPaEqual_54(struct HashMapIntToIntArr *map);
void testPaEqual_55(struct HashMapIntToIntArr *map);
void testPaEqual_56(struct HashMapIntToIntArr *map);
void testPaEqual_57(struct HashMapIntToIntArr *map);
void testPaEqual_58(struct HashMapIntToIntArr *map);
void testPaEqual_59(struct HashMapIntToIntArr *map);
void testPaEqual_60(struct HashMapIntToIntArr *map);
void testPaEqual_61(struct HashMapIntToIntArr *map);
void testPaEqual_62(struct HashMapIntToIntArr *map);
void testPaEqual_63(struct HashMapIntToIntArr *map);
void testPaEqual_64(struct HashMapIntToIntArr *map);

void testNumOverlap(void) {
  struct HashMapIntToIntArr *map = initHashMapIntToIntArrMDI();

  testPaEqual_1(map);
  testPaEqual_2(map);
  testPaEqual_3(map);
  testPaEqual_4(map);
  testPaEqual_5(map);
  testPaEqual_6(map);
  testPaEqual_7(map);
  testPaEqual_8(map);
  testPaEqual_9(map);
  testPaEqual_10(map);
  testPaEqual_11(map);
  testPaEqual_12(map);
  testPaEqual_13(map);
  testPaEqual_14(map);
  testPaEqual_15(map);
  testPaEqual_16(map);
  testPaEqual_17(map);
  testPaEqual_18(map);
  testPaEqual_19(map);
  testPaEqual_20(map);
  testPaEqual_21(map);
  testPaEqual_22(map);
  testPaEqual_23(map);
  testPaEqual_24(map);
  testPaEqual_25(map);
  testPaEqual_26(map);
  testPaEqual_27(map);
  testPaEqual_28(map);
  testPaEqual_29(map);
  testPaEqual_30(map);
  testPaEqual_31(map);
  testPaEqual_32(map);
  testPaEqual_33(map);
  testPaEqual_34(map);
  testPaEqual_35(map);
  testPaEqual_36(map);
  testPaEqual_37(map);
  testPaEqual_38(map);
  testPaEqual_39(map);
  testPaEqual_40(map);
  testPaEqual_41(map);
  testPaEqual_42(map);
  testPaEqual_43(map);
  testPaEqual_44(map);
  testPaEqual_45(map);
  testPaEqual_46(map);
  testPaEqual_47(map);
  testPaEqual_48(map);
  testPaEqual_49(map);
  testPaEqual_50(map);
  testPaEqual_51(map);
  testPaEqual_52(map);
  testPaEqual_53(map);
  testPaEqual_54(map);
  testPaEqual_55(map);
  testPaEqual_56(map);
  testPaEqual_57(map);
  testPaEqual_58(map);
  testPaEqual_59(map);
  testPaEqual_60(map);
  testPaEqual_61(map);
  testPaEqual_62(map);
  testPaEqual_63(map);
  testPaEqual_64(map);
}

int main(void) {
  UNITY_BEGIN();

  RUN_TEST(testNumOverlap);

  return UNITY_END();
}

void testPaEqual_1(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 1), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 2), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 3), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 4), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 5), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 6), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 7), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 8), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 9), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 0),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_2(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 2), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 3), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 4), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 5), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 6), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 7), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 8), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 9), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 1),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_3(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 3), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 4), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 5), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 6), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 7), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 8), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 9), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 2),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_4(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 4), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 5), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 6), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 7), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 8), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 9), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 3),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_5(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 5), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 6), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 7), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 8), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 9), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 4),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_6(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 6), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 7), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 8), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 9), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 5),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_7(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 7), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 8), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 9), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 6),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_8(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 8), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 9), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 7),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_9(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 9), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 8),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_10(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 10), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 9),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_11(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 11), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 10),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_12(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 12), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 11),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_13(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 13), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 12),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_14(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 14), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 13),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_15(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 15), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 14),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_16(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 16), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 15),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_17(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 17), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 16),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_18(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 18), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 17),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_19(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 19), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 18),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_20(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 20), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 19),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_21(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 21), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 20),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_22(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 22), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 21),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_23(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 23), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 22),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_24(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 24), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 23),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_25(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 25), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 24),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_26(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 26), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 25),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_27(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 27), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 26),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_28(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 28), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 27),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_29(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 29), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 28),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_30(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 30), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 29),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_31(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 31), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 30),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_32(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 32), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 31),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_33(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 33), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 32),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_34(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 34), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(3, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 33),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_35(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 35), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 34),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_36(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 36), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 35),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_37(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 37), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 36),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_38(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 38), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 37),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_39(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 39), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 38),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_40(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 40), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 39),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_41(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 41), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 40),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_42(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 42), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 41),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_43(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 43), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 42),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_44(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 44), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 43),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_45(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 45), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 44),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_46(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 46), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 45),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_47(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 47), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 46),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_48(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 48), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 47),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_49(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 49), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 48),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_50(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 50), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 49),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_51(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 51), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 50),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_52(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 52), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 51),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_53(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 53), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 52),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_54(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 54), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(2, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(01, NumOverlap(GetValIntToIntArr(map, 53),
                                       GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 53),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_55(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 55), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 54),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_56(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 55),
                                      GetValIntToIntArr(map, 56), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 55),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 55),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 55),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 55),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 55),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 55),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 55),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 55),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_57(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 56),
                                      GetValIntToIntArr(map, 57), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 56),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 56),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 56),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 56),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 56),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 56),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 56),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_58(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 57),
                                      GetValIntToIntArr(map, 58), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 57),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 57),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 57),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 57),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 57),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 57),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_59(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 58),
                                      GetValIntToIntArr(map, 59), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 58),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 58),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 58),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 58),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 58),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_60(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 59),
                                      GetValIntToIntArr(map, 60), 3, 0, 0));

  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 59),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 59),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(1, NumOverlap(GetValIntToIntArr(map, 59),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 59),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_61(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 60),
                                      GetValIntToIntArr(map, 61), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 60),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 60),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 60),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_62(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 61),
                                      GetValIntToIntArr(map, 62), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 61),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 61),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_63(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 62),
                                      GetValIntToIntArr(map, 63), 3, 0, 0));
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 62),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}

void testPaEqual_64(struct HashMapIntToIntArr *map) {
  TEST_ASSERT_EQUAL_INT(0, NumOverlap(GetValIntToIntArr(map, 63),
                                      GetValIntToIntArr(map, 64), 3, 0, 0));
}
