//
// Created by Senbai Kang on 11/25/19.
//

#include "testFunctions.h"
#include <unity.c>
#include "testGetMuDelInsNr.h"

void testPa49Wrapper(void);
void testPa50Wrapper(void);
void testPa51Wrapper(void);
void testPa52Wrapper(void);
void testPa53Wrapper(void);
void testPa54Wrapper(void);
void testPa55Wrapper(void);
void testPa56Wrapper(void);
void testPa57Wrapper(void);
void testPa58Wrapper(void);
void testPa59Wrapper(void);
void testPa60Wrapper(void);
void testPa61Wrapper(void);
void testPa62Wrapper(void);
void testPa63Wrapper(void);
void testPa64Wrapper(void);
void testPa65Wrapper(void);

int main(void) {
  map = initHashMapIntToIntArrMDI();

  dataPa = (DataStr *)calloc(1, sizeof(DataStr));
  dataPa->alleles = (int *)calloc(3, sizeof(int));
  dataCh = (DataStr *)calloc(1, sizeof(DataStr));
  dataCh->alleles = (int *)calloc(3, sizeof(int));

  UNITY_BEGIN();

  RUN_TEST(testPa49Wrapper);
  RUN_TEST(testPa50Wrapper);
  RUN_TEST(testPa51Wrapper);
  RUN_TEST(testPa52Wrapper);
  RUN_TEST(testPa53Wrapper);
  RUN_TEST(testPa54Wrapper);
  RUN_TEST(testPa55Wrapper);
  RUN_TEST(testPa56Wrapper);
  RUN_TEST(testPa57Wrapper);
  RUN_TEST(testPa58Wrapper);
  RUN_TEST(testPa59Wrapper);
  RUN_TEST(testPa60Wrapper);
  RUN_TEST(testPa61Wrapper);
  RUN_TEST(testPa62Wrapper);
  RUN_TEST(testPa63Wrapper);
  RUN_TEST(testPa64Wrapper);
  RUN_TEST(testPa65Wrapper);

  return UNITY_END();
}

void testPa49Wrapper(void) {
  SetDataStrMDI(map, dataPa, 48);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(1, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(1, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(1, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(1, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(1, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(1, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(1, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(1, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(1, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(1, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(1, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(1, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(1, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(1, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(1, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(1, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(1, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(1, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(1, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(1, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(1, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(1, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(1, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(1, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(1, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(1, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(1, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(1, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(1, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(1, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(1, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(1, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(1, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(1, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(1, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(1, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(1, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(1, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(1, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(1, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c42[0]);
  TEST_ASSERT_EQUAL_INT(1, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c46[0]);
  TEST_ASSERT_EQUAL_INT(1, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c47[0]);
  TEST_ASSERT_EQUAL_INT(1, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c48[0]);
  TEST_ASSERT_EQUAL_INT(1, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(1, c54[1]);
  TEST_ASSERT_EQUAL_INT(1, c54[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(1, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(1, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c57[0]);
  TEST_ASSERT_EQUAL_INT(1, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c58[0]);
  TEST_ASSERT_EQUAL_INT(1, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c59[0]);
  TEST_ASSERT_EQUAL_INT(1, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa50Wrapper(void) {
  SetDataStrMDI(map, dataPa, 49);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(1, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(1, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(1, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(1, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(1, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(1, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(1, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(1, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(1, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(1, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(1, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(1, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(1, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(1, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(1, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(1, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(1, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(1, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(1, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(1, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(1, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(1, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(1, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(1, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(1, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(1, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(1, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(1, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(1, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(1, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(1, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(1, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(1, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(1, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(1, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(1, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(1, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(1, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(1, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(1, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c42[0]);
  TEST_ASSERT_EQUAL_INT(1, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c46[0]);
  TEST_ASSERT_EQUAL_INT(1, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c47[0]);
  TEST_ASSERT_EQUAL_INT(1, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c48[0]);
  TEST_ASSERT_EQUAL_INT(1, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(1, c54[1]);
  TEST_ASSERT_EQUAL_INT(1, c54[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(1, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(1, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c57[0]);
  TEST_ASSERT_EQUAL_INT(1, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(1, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c59[0]);
  TEST_ASSERT_EQUAL_INT(1, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa51Wrapper(void) {
  SetDataStrMDI(map, dataPa, 50);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(1, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(1, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(1, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(1, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(1, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(1, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(1, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(1, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(1, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(1, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(1, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(1, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(1, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(1, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(1, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(1, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(1, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(1, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(1, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(1, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(1, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(1, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(1, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(1, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(1, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(1, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(1, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(1, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(1, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(1, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(1, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(1, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(1, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(1, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(1, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(1, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(1, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(1, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(1, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(1, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c42[0]);
  TEST_ASSERT_EQUAL_INT(1, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c46[0]);
  TEST_ASSERT_EQUAL_INT(1, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c47[0]);
  TEST_ASSERT_EQUAL_INT(1, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c48[0]);
  TEST_ASSERT_EQUAL_INT(1, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(1, c54[1]);
  TEST_ASSERT_EQUAL_INT(1, c54[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(1, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(1, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c57[0]);
  TEST_ASSERT_EQUAL_INT(1, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(1, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c59[0]);
  TEST_ASSERT_EQUAL_INT(1, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa52Wrapper(void) {
  SetDataStrMDI(map, dataPa, 51);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(1, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(1, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(1, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(1, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(1, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(1, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(1, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(1, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(1, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(1, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(1, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(1, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(1, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(1, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(1, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(1, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(1, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(1, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(1, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(1, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(1, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(1, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(1, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(1, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(1, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(1, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(1, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(1, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(1, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(1, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(1, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(1, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(1, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(1, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(1, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(1, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(1, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(1, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(1, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(1, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c42[0]);
  TEST_ASSERT_EQUAL_INT(1, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c46[0]);
  TEST_ASSERT_EQUAL_INT(1, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c47[0]);
  TEST_ASSERT_EQUAL_INT(1, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c48[0]);
  TEST_ASSERT_EQUAL_INT(1, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c54[0]);
  TEST_ASSERT_EQUAL_INT(1, c54[1]);
  TEST_ASSERT_EQUAL_INT(1, c54[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(1, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c56[0]);
  TEST_ASSERT_EQUAL_INT(1, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c57[0]);
  TEST_ASSERT_EQUAL_INT(1, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(1, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c59[0]);
  TEST_ASSERT_EQUAL_INT(1, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa53Wrapper(void) {
  SetDataStrMDI(map, dataPa, 52);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(1, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(1, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(1, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(1, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(1, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(1, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(1, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(1, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(1, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(1, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(1, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(1, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(1, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(1, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(1, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(1, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(1, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(1, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(1, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(1, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(1, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(1, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(1, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(1, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(1, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(1, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(1, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(1, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(1, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(1, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(1, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(1, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(1, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(1, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(1, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(1, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(1, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(1, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(1, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(1, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c42[0]);
  TEST_ASSERT_EQUAL_INT(1, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c46[0]);
  TEST_ASSERT_EQUAL_INT(1, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c47[0]);
  TEST_ASSERT_EQUAL_INT(1, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c48[0]);
  TEST_ASSERT_EQUAL_INT(1, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c54[0]);
  TEST_ASSERT_EQUAL_INT(1, c54[1]);
  TEST_ASSERT_EQUAL_INT(1, c54[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(1, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(1, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c57[0]);
  TEST_ASSERT_EQUAL_INT(1, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(1, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c59[0]);
  TEST_ASSERT_EQUAL_INT(1, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa54Wrapper(void) {
  SetDataStrMDI(map, dataPa, 53);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(1, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(1, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(1, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(1, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(1, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(1, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(1, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(1, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(1, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(1, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(1, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(1, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(1, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(1, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(1, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(1, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(1, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(1, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(1, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(1, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(1, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(1, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(1, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(1, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(1, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(1, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(1, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(1, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(1, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(3, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(1, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(1, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(1, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(1, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(1, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(1, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(1, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(1, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(1, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(1, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(1, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c42[0]);
  TEST_ASSERT_EQUAL_INT(1, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c46[0]);
  TEST_ASSERT_EQUAL_INT(1, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c47[0]);
  TEST_ASSERT_EQUAL_INT(1, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c48[0]);
  TEST_ASSERT_EQUAL_INT(1, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c54[0]);
  TEST_ASSERT_EQUAL_INT(1, c54[1]);
  TEST_ASSERT_EQUAL_INT(1, c54[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c55[0]);
  TEST_ASSERT_EQUAL_INT(1, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(1, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c57[0]);
  TEST_ASSERT_EQUAL_INT(1, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(1, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c59[0]);
  TEST_ASSERT_EQUAL_INT(1, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa55Wrapper(void) {
  SetDataStrMDI(map, dataPa, 54);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(0, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(0, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(0, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(0, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(0, c55[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(0, c56[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(0, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(0, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(0, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa56Wrapper(void) {
  SetDataStrMDI(map, dataPa, 55);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(0, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(0, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(0, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(0, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(0, c55[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(0, c56[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(0, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(0, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(0, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa57Wrapper(void) {
  SetDataStrMDI(map, dataPa, 56);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(0, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(0, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(0, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(0, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(0, c55[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(0, c56[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(0, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(0, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(0, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa58Wrapper(void) {
  SetDataStrMDI(map, dataPa, 57);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(0, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(0, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(0, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(0, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(0, c55[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(0, c56[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(0, c57[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(0, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(0, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa59Wrapper(void) {
  SetDataStrMDI(map, dataPa, 58);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(0, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(0, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(0, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(0, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(0, c55[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(0, c56[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(0, c57[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(0, c58[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(0, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa60Wrapper(void) {
  SetDataStrMDI(map, dataPa, 59);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(0, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(0, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(0, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(0, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(0, c55[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(0, c56[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(0, c57[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(0, c58[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(0, c59[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(1, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(1, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(1, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c63[0]);
  TEST_ASSERT_EQUAL_INT(1, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(2, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa61Wrapper(void) {
  SetDataStrMDI(map, dataPa, 60);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(0, c60[1]);
  TEST_ASSERT_EQUAL_INT(1, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(0, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(0, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(0, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(1, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa62Wrapper(void) {
  SetDataStrMDI(map, dataPa, 61);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(0, c60[1]);
  TEST_ASSERT_EQUAL_INT(1, c60[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(0, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(0, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(0, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(1, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa63Wrapper(void) {
  SetDataStrMDI(map, dataPa, 62);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c60[0]);
  TEST_ASSERT_EQUAL_INT(0, c60[1]);
  TEST_ASSERT_EQUAL_INT(1, c60[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(0, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(0, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(0, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(1, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa64Wrapper(void) {
  SetDataStrMDI(map, dataPa, 63);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(1, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(1, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(1, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(1, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(1, c55[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(2, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(1, c56[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(1, c57[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(1, c58[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(1, c59[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c60[0]);
  TEST_ASSERT_EQUAL_INT(0, c60[1]);
  TEST_ASSERT_EQUAL_INT(1, c60[2]);

  SetDataStrMDI(map, dataCh, 60);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c61[0]);
  TEST_ASSERT_EQUAL_INT(0, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c62[0]);
  TEST_ASSERT_EQUAL_INT(0, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(1, c63[0]);
  TEST_ASSERT_EQUAL_INT(0, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(1, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}

void testPa65Wrapper(void) {
  SetDataStrMDI(map, dataPa, 64);

  // tested in the order of number of mutations, deletions, and insertions
  SetDataStrMDI(map, dataCh, 0);
  int *c1 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c1[0]);
  TEST_ASSERT_EQUAL_INT(0, c1[1]);
  TEST_ASSERT_EQUAL_INT(0, c1[2]);

  SetDataStrMDI(map, dataCh, 1);
  int *c2 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c2[0]);
  TEST_ASSERT_EQUAL_INT(0, c2[1]);
  TEST_ASSERT_EQUAL_INT(0, c2[2]);

  SetDataStrMDI(map, dataCh, 2);
  int *c3 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c3[0]);
  TEST_ASSERT_EQUAL_INT(0, c3[1]);
  TEST_ASSERT_EQUAL_INT(0, c3[2]);

  SetDataStrMDI(map, dataCh, 3);
  int *c4 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c4[0]);
  TEST_ASSERT_EQUAL_INT(0, c4[1]);
  TEST_ASSERT_EQUAL_INT(0, c4[2]);

  SetDataStrMDI(map, dataCh, 4);
  int *c5 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c5[0]);
  TEST_ASSERT_EQUAL_INT(0, c5[1]);
  TEST_ASSERT_EQUAL_INT(0, c5[2]);

  SetDataStrMDI(map, dataCh, 5);
  int *c6 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c6[0]);
  TEST_ASSERT_EQUAL_INT(0, c6[1]);
  TEST_ASSERT_EQUAL_INT(0, c6[2]);

  SetDataStrMDI(map, dataCh, 6);
  int *c7 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c7[0]);
  TEST_ASSERT_EQUAL_INT(0, c7[1]);
  TEST_ASSERT_EQUAL_INT(0, c7[2]);

  SetDataStrMDI(map, dataCh, 7);
  int *c8 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c8[0]);
  TEST_ASSERT_EQUAL_INT(0, c8[1]);
  TEST_ASSERT_EQUAL_INT(0, c8[2]);

  SetDataStrMDI(map, dataCh, 8);
  int *c9 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c9[0]);
  TEST_ASSERT_EQUAL_INT(0, c9[1]);
  TEST_ASSERT_EQUAL_INT(0, c9[2]);

  SetDataStrMDI(map, dataCh, 9);
  int *c10 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c10[0]);
  TEST_ASSERT_EQUAL_INT(0, c10[1]);
  TEST_ASSERT_EQUAL_INT(0, c10[2]);

  SetDataStrMDI(map, dataCh, 10);
  int *c11 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c11[0]);
  TEST_ASSERT_EQUAL_INT(0, c11[1]);
  TEST_ASSERT_EQUAL_INT(0, c11[2]);

  SetDataStrMDI(map, dataCh, 11);
  int *c12 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c12[0]);
  TEST_ASSERT_EQUAL_INT(0, c12[1]);
  TEST_ASSERT_EQUAL_INT(0, c12[2]);

  SetDataStrMDI(map, dataCh, 12);
  int *c13 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c13[0]);
  TEST_ASSERT_EQUAL_INT(0, c13[1]);
  TEST_ASSERT_EQUAL_INT(0, c13[2]);

  SetDataStrMDI(map, dataCh, 13);
  int *c14 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c14[0]);
  TEST_ASSERT_EQUAL_INT(0, c14[1]);
  TEST_ASSERT_EQUAL_INT(0, c14[2]);

  SetDataStrMDI(map, dataCh, 14);
  int *c15 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c15[0]);
  TEST_ASSERT_EQUAL_INT(0, c15[1]);
  TEST_ASSERT_EQUAL_INT(0, c15[2]);

  SetDataStrMDI(map, dataCh, 15);
  int *c16 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c16[0]);
  TEST_ASSERT_EQUAL_INT(0, c16[1]);
  TEST_ASSERT_EQUAL_INT(0, c16[2]);

  SetDataStrMDI(map, dataCh, 16);
  int *c17 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c17[0]);
  TEST_ASSERT_EQUAL_INT(0, c17[1]);
  TEST_ASSERT_EQUAL_INT(0, c17[2]);

  SetDataStrMDI(map, dataCh, 17);
  int *c18 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c18[0]);
  TEST_ASSERT_EQUAL_INT(0, c18[1]);
  TEST_ASSERT_EQUAL_INT(0, c18[2]);

  SetDataStrMDI(map, dataCh, 18);
  int *c19 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c19[0]);
  TEST_ASSERT_EQUAL_INT(0, c19[1]);
  TEST_ASSERT_EQUAL_INT(0, c19[2]);

  SetDataStrMDI(map, dataCh, 19);
  int *c20 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c20[0]);
  TEST_ASSERT_EQUAL_INT(0, c20[1]);
  TEST_ASSERT_EQUAL_INT(0, c20[2]);

  SetDataStrMDI(map, dataCh, 20);
  int *c21 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c21[0]);
  TEST_ASSERT_EQUAL_INT(0, c21[1]);
  TEST_ASSERT_EQUAL_INT(0, c21[2]);

  SetDataStrMDI(map, dataCh, 21);
  int *c22 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c22[0]);
  TEST_ASSERT_EQUAL_INT(0, c22[1]);
  TEST_ASSERT_EQUAL_INT(0, c22[2]);

  SetDataStrMDI(map, dataCh, 22);
  int *c23 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c23[0]);
  TEST_ASSERT_EQUAL_INT(0, c23[1]);
  TEST_ASSERT_EQUAL_INT(0, c23[2]);

  SetDataStrMDI(map, dataCh, 23);
  int *c24 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c24[0]);
  TEST_ASSERT_EQUAL_INT(0, c24[1]);
  TEST_ASSERT_EQUAL_INT(0, c24[2]);

  SetDataStrMDI(map, dataCh, 24);
  int *c25 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c25[0]);
  TEST_ASSERT_EQUAL_INT(0, c25[1]);
  TEST_ASSERT_EQUAL_INT(0, c25[2]);

  SetDataStrMDI(map, dataCh, 25);
  int *c26 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c26[0]);
  TEST_ASSERT_EQUAL_INT(0, c26[1]);
  TEST_ASSERT_EQUAL_INT(0, c26[2]);

  SetDataStrMDI(map, dataCh, 26);
  int *c27 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c27[0]);
  TEST_ASSERT_EQUAL_INT(0, c27[1]);
  TEST_ASSERT_EQUAL_INT(0, c27[2]);

  SetDataStrMDI(map, dataCh, 27);
  int *c28 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c28[0]);
  TEST_ASSERT_EQUAL_INT(0, c28[1]);
  TEST_ASSERT_EQUAL_INT(0, c28[2]);

  SetDataStrMDI(map, dataCh, 28);
  int *c29 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c29[0]);
  TEST_ASSERT_EQUAL_INT(0, c29[1]);
  TEST_ASSERT_EQUAL_INT(0, c29[2]);

  SetDataStrMDI(map, dataCh, 29);
  int *c30 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c30[0]);
  TEST_ASSERT_EQUAL_INT(0, c30[1]);
  TEST_ASSERT_EQUAL_INT(0, c30[2]);

  SetDataStrMDI(map, dataCh, 30);
  int *c31 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c31[0]);
  TEST_ASSERT_EQUAL_INT(0, c31[1]);
  TEST_ASSERT_EQUAL_INT(0, c31[2]);

  SetDataStrMDI(map, dataCh, 31);
  int *c32 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c32[0]);
  TEST_ASSERT_EQUAL_INT(0, c32[1]);
  TEST_ASSERT_EQUAL_INT(0, c32[2]);

  SetDataStrMDI(map, dataCh, 32);
  int *c33 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c33[0]);
  TEST_ASSERT_EQUAL_INT(0, c33[1]);
  TEST_ASSERT_EQUAL_INT(0, c33[2]);

  SetDataStrMDI(map, dataCh, 33);
  int *c34 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c34[0]);
  TEST_ASSERT_EQUAL_INT(0, c34[1]);
  TEST_ASSERT_EQUAL_INT(0, c34[2]);

  SetDataStrMDI(map, dataCh, 34);
  int *c35 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c35[0]);
  TEST_ASSERT_EQUAL_INT(0, c35[1]);
  TEST_ASSERT_EQUAL_INT(0, c35[2]);

  SetDataStrMDI(map, dataCh, 35);
  int *c36 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c36[0]);
  TEST_ASSERT_EQUAL_INT(0, c36[1]);
  TEST_ASSERT_EQUAL_INT(0, c36[2]);

  SetDataStrMDI(map, dataCh, 36);
  int *c37 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c37[0]);
  TEST_ASSERT_EQUAL_INT(0, c37[1]);
  TEST_ASSERT_EQUAL_INT(0, c37[2]);

  SetDataStrMDI(map, dataCh, 37);
  int *c38 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c38[0]);
  TEST_ASSERT_EQUAL_INT(0, c38[1]);
  TEST_ASSERT_EQUAL_INT(0, c38[2]);

  SetDataStrMDI(map, dataCh, 38);
  int *c39 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c39[0]);
  TEST_ASSERT_EQUAL_INT(0, c39[1]);
  TEST_ASSERT_EQUAL_INT(0, c39[2]);

  SetDataStrMDI(map, dataCh, 39);
  int *c40 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c40[0]);
  TEST_ASSERT_EQUAL_INT(0, c40[1]);
  TEST_ASSERT_EQUAL_INT(0, c40[2]);

  SetDataStrMDI(map, dataCh, 40);
  int *c41 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c41[0]);
  TEST_ASSERT_EQUAL_INT(0, c41[1]);
  TEST_ASSERT_EQUAL_INT(0, c41[2]);

  SetDataStrMDI(map, dataCh, 41);
  int *c42 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c42[0]);
  TEST_ASSERT_EQUAL_INT(0, c42[1]);
  TEST_ASSERT_EQUAL_INT(0, c42[2]);

  SetDataStrMDI(map, dataCh, 42);
  int *c43 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c43[0]);
  TEST_ASSERT_EQUAL_INT(0, c43[1]);
  TEST_ASSERT_EQUAL_INT(0, c43[2]);

  SetDataStrMDI(map, dataCh, 43);
  int *c44 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c44[0]);
  TEST_ASSERT_EQUAL_INT(0, c44[1]);
  TEST_ASSERT_EQUAL_INT(0, c44[2]);

  SetDataStrMDI(map, dataCh, 44);
  int *c45 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c45[0]);
  TEST_ASSERT_EQUAL_INT(0, c45[1]);
  TEST_ASSERT_EQUAL_INT(0, c45[2]);

  SetDataStrMDI(map, dataCh, 45);
  int *c46 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c46[0]);
  TEST_ASSERT_EQUAL_INT(0, c46[1]);
  TEST_ASSERT_EQUAL_INT(0, c46[2]);

  SetDataStrMDI(map, dataCh, 46);
  int *c47 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c47[0]);
  TEST_ASSERT_EQUAL_INT(0, c47[1]);
  TEST_ASSERT_EQUAL_INT(0, c47[2]);

  SetDataStrMDI(map, dataCh, 47);
  int *c48 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c48[0]);
  TEST_ASSERT_EQUAL_INT(0, c48[1]);
  TEST_ASSERT_EQUAL_INT(0, c48[2]);

  SetDataStrMDI(map, dataCh, 48);
  int *c49 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c49[0]);
  TEST_ASSERT_EQUAL_INT(0, c49[1]);
  TEST_ASSERT_EQUAL_INT(0, c49[2]);

  SetDataStrMDI(map, dataCh, 49);
  int *c50 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c50[0]);
  TEST_ASSERT_EQUAL_INT(0, c50[1]);
  TEST_ASSERT_EQUAL_INT(0, c50[2]);

  SetDataStrMDI(map, dataCh, 50);
  int *c51 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c51[0]);
  TEST_ASSERT_EQUAL_INT(0, c51[1]);
  TEST_ASSERT_EQUAL_INT(0, c51[2]);

  SetDataStrMDI(map, dataCh, 51);
  int *c52 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c52[0]);
  TEST_ASSERT_EQUAL_INT(0, c52[1]);
  TEST_ASSERT_EQUAL_INT(0, c52[2]);

  SetDataStrMDI(map, dataCh, 52);
  int *c53 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c53[0]);
  TEST_ASSERT_EQUAL_INT(0, c53[1]);
  TEST_ASSERT_EQUAL_INT(0, c53[2]);

  SetDataStrMDI(map, dataCh, 53);
  int *c54 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c54[0]);
  TEST_ASSERT_EQUAL_INT(0, c54[1]);
  TEST_ASSERT_EQUAL_INT(0, c54[2]);

  SetDataStrMDI(map, dataCh, 54);
  int *c55 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c55[0]);
  TEST_ASSERT_EQUAL_INT(0, c55[1]);
  TEST_ASSERT_EQUAL_INT(0, c55[2]);

  SetDataStrMDI(map, dataCh, 55);
  int *c56 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c56[0]);
  TEST_ASSERT_EQUAL_INT(0, c56[1]);
  TEST_ASSERT_EQUAL_INT(0, c56[2]);

  SetDataStrMDI(map, dataCh, 56);
  int *c57 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c57[0]);
  TEST_ASSERT_EQUAL_INT(0, c57[1]);
  TEST_ASSERT_EQUAL_INT(0, c57[2]);

  SetDataStrMDI(map, dataCh, 57);
  int *c58 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c58[0]);
  TEST_ASSERT_EQUAL_INT(0, c58[1]);
  TEST_ASSERT_EQUAL_INT(0, c58[2]);

  SetDataStrMDI(map, dataCh, 58);
  int *c59 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c59[0]);
  TEST_ASSERT_EQUAL_INT(0, c59[1]);
  TEST_ASSERT_EQUAL_INT(0, c59[2]);

  SetDataStrMDI(map, dataCh, 59);
  int *c60 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c60[0]);
  TEST_ASSERT_EQUAL_INT(0, c60[1]);
  TEST_ASSERT_EQUAL_INT(0, c60[2]);

  SetDataStrMDI(map, dataCh, 61);
  int *c61 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c61[0]);
  TEST_ASSERT_EQUAL_INT(0, c61[1]);
  TEST_ASSERT_EQUAL_INT(0, c61[2]);

  SetDataStrMDI(map, dataCh, 62);
  int *c62 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c62[0]);
  TEST_ASSERT_EQUAL_INT(0, c62[1]);
  TEST_ASSERT_EQUAL_INT(0, c62[2]);

  SetDataStrMDI(map, dataCh, 63);
  int *c63 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c63[0]);
  TEST_ASSERT_EQUAL_INT(0, c63[1]);
  TEST_ASSERT_EQUAL_INT(0, c63[2]);

  SetDataStrMDI(map, dataCh, 64);
  int *c64 = GetMuDelInsNrCore(dataPa, dataCh, 3);
  TEST_ASSERT_EQUAL_INT(0, c64[0]);
  TEST_ASSERT_EQUAL_INT(0, c64[1]);
  TEST_ASSERT_EQUAL_INT(0, c64[2]);
}