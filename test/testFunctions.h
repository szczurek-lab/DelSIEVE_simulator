//
// Created by Senbai Kang on 11/25/19.
//

#include "../src/finiteMuInDelSimulator.h"

#ifndef MYSIMULATOR_TESTFUNCTIONS_H
#define MYSIMULATOR_TESTFUNCTIONS_H

int **genoArr = NULL;
struct HashMapIntToIntArr *hashMapIntToIntArr = NULL;
struct HashMapIntArrToInt *hashMapIntArrToInt = NULL;

void AllocateGenoArr(int n) {
  genoArr = (int **)calloc(n + 1, sizeof(int *));
  if (!genoArr) {
    fprintf(stderr, "Could not allocate the genoArr structure\n");
    exit(-1);
  }
  for (int i = 0; i < n; ++i) {
    genoArr[i] = (int *)calloc(4, sizeof(int));
    if (!genoArr[i]) {
      fprintf(stderr, "Could not allocate the genoArr[] structure\n");
      exit(-1);
    }
  }
}

void FreeGenoArr(int n) {
  for (int i = 0; i < n; ++i) {
    free(genoArr[i]);
    genoArr[i] = NULL;
  }
  free(genoArr);
  genoArr = NULL;
}

struct HashMapIntToIntArr *initHashMapIntToIntArrMD() {
  int n = 15;
  double muRate = 1.0;
  double delRate = 1.0;
  double Mij[n][n];

  if (genoArr != NULL) {
    FreeGenoArr(n);
  }
  AllocateGenoArr(n);

  // normalize the deletionRate and insRate regarding to the mutationRate
  delRate /= (muRate + 1e-15);

  // initialize all entries to 0.0
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      Mij[i][j] = 0.0;
    }

  // 0    1    2    3    4    5    6    7    8    9    10   11   12   13   14
  // A/A, A/C, A/G, A/T, C/C, C/G, C/T, G/G, G/T, T/T, A/-, C/-, G/-, T/-, -

  // 0, A/A
  Mij[0][1] = Mij[0][2] = Mij[0][3] = 1.0 / 3;
  Mij[0][10] = delRate;

  genoArr[0][0] = 0;
  genoArr[0][1] = 0;
  genoArr[0][2] = INVALID;

  // 1, A/C
  Mij[1][0] = Mij[1][2] = Mij[1][3] = Mij[1][4] = Mij[1][5] = Mij[1][6] =
      1.0 / 6;
  Mij[1][10] = Mij[1][11] = delRate / 2;

  genoArr[1][0] = 0;
  genoArr[1][1] = 1;
  genoArr[1][2] = INVALID;

  // 2, A/G
  Mij[2][0] = Mij[2][1] = Mij[2][3] = Mij[2][5] = Mij[2][7] = Mij[2][8] =
      1.0 / 6;
  Mij[2][10] = Mij[2][12] = delRate / 2;

  genoArr[2][0] = 0;
  genoArr[2][1] = 2;
  genoArr[2][2] = INVALID;

  // 3, A/T
  Mij[3][0] = Mij[3][1] = Mij[3][2] = Mij[3][6] = Mij[3][8] = Mij[3][9] =
      1.0 / 6;
  Mij[3][10] = Mij[3][13] = delRate / 2;

  genoArr[3][0] = 0;
  genoArr[3][1] = 3;
  genoArr[3][2] = INVALID;

  // 4, C/C
  Mij[4][1] = Mij[4][5] = Mij[4][6] = 1.0 / 3;
  Mij[4][11] = delRate;

  genoArr[4][0] = 1;
  genoArr[4][1] = 1;
  genoArr[4][2] = INVALID;

  // 5, C/G
  Mij[5][1] = Mij[5][2] = Mij[5][4] = Mij[5][6] = Mij[5][7] = Mij[5][8] =
      1.0 / 6;
  Mij[5][11] = Mij[5][12] = delRate / 2;

  genoArr[5][0] = 1;
  genoArr[5][1] = 2;
  genoArr[5][2] = INVALID;

  // 6, C/T
  Mij[6][1] = Mij[6][3] = Mij[6][4] = Mij[6][5] = Mij[6][8] = Mij[6][9] =
      1.0 / 6;
  Mij[6][11] = Mij[6][13] = delRate / 2;

  genoArr[6][0] = 1;
  genoArr[6][1] = 3;
  genoArr[6][2] = INVALID;

  // 7, G/G
  Mij[7][2] = Mij[7][5] = Mij[7][8] = 1.0 / 3;
  Mij[7][12] = delRate;

  genoArr[7][0] = 2;
  genoArr[7][1] = 2;
  genoArr[7][2] = INVALID;

  // 8, G/T
  Mij[8][2] = Mij[8][3] = Mij[8][5] = Mij[8][6] = Mij[8][7] = Mij[8][9] =
      1.0 / 6;
  Mij[8][12] = Mij[8][13] = delRate / 2;

  genoArr[8][0] = 2;
  genoArr[8][1] = 3;
  genoArr[8][2] = INVALID;

  // 9, T/T
  Mij[9][3] = Mij[9][6] = Mij[9][8] = 1.0 / 3;
  Mij[9][13] = delRate;

  genoArr[9][0] = 3;
  genoArr[9][1] = 3;
  genoArr[9][2] = INVALID;

  // 10, A/-
  Mij[10][11] = Mij[10][12] = Mij[10][13] = 1.0 / 6;
  Mij[10][14] = delRate / 2;

  genoArr[10][0] = 0;
  genoArr[10][1] = DELETED;
  genoArr[10][2] = INVALID;

  // 11, C/-
  Mij[11][10] = Mij[11][12] = Mij[11][13] = 1.0 / 6;
  Mij[11][14] = delRate / 2;

  genoArr[11][0] = 1;
  genoArr[11][1] = DELETED;
  genoArr[11][2] = INVALID;

  // 12, G/-
  Mij[12][10] = Mij[12][11] = Mij[12][13] = 1.0 / 6;
  Mij[12][14] = delRate / 2;

  genoArr[12][0] = 2;
  genoArr[12][1] = DELETED;
  genoArr[12][2] = INVALID;

  // 13, T/-
  Mij[13][10] = Mij[13][11] = Mij[13][12] = 1.0 / 6;
  Mij[13][14] = delRate / 2;

  genoArr[13][0] = 3;
  genoArr[13][1] = DELETED;
  genoArr[13][2] = INVALID;

  // 14, -
  genoArr[14][0] = genoArr[14][1] = DELETED;
  genoArr[14][2] = INVALID;

  hashMapIntToIntArr = CreateHashMapIntToIntArr(3);
  for (int i = 0; i < n; ++i)
    InsertKeyValIntToIntArr(hashMapIntToIntArr, i, genoArr[i]);

  return hashMapIntToIntArr;
}

struct HashMapIntArrToInt *initHashMapIntArrToIntMD() {
  int n = 15;
  double muRate = 1.0;
  double delRate = 1.0;
  double Mij[n][n];

  if (genoArr != NULL) {
    FreeGenoArr(n);
  }
  AllocateGenoArr(n);

  // normalize the deletionRate and insRate regarding to the mutationRate
  delRate /= (muRate + 1e-15);

  // initialize all entries to 0.0
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      Mij[i][j] = 0.0;
    }

  // 0    1    2    3    4    5    6    7    8    9    10   11   12   13   14
  // A/A, A/C, A/G, A/T, C/C, C/G, C/T, G/G, G/T, T/T, A/-, C/-, G/-, T/-, -

  // 0, A/A
  Mij[0][1] = Mij[0][2] = Mij[0][3] = 1.0 / 3;
  Mij[0][10] = delRate;

  genoArr[0][0] = 0;
  genoArr[0][1] = 0;
  genoArr[0][2] = INVALID;

  // 1, A/C
  Mij[1][0] = Mij[1][2] = Mij[1][3] = Mij[1][4] = Mij[1][5] = Mij[1][6] =
      1.0 / 6;
  Mij[1][10] = Mij[1][11] = delRate / 2;

  genoArr[1][0] = 0;
  genoArr[1][1] = 1;
  genoArr[1][2] = INVALID;

  // 2, A/G
  Mij[2][0] = Mij[2][1] = Mij[2][3] = Mij[2][5] = Mij[2][7] = Mij[2][8] =
      1.0 / 6;
  Mij[2][10] = Mij[2][12] = delRate / 2;

  genoArr[2][0] = 0;
  genoArr[2][1] = 2;
  genoArr[2][2] = INVALID;

  // 3, A/T
  Mij[3][0] = Mij[3][1] = Mij[3][2] = Mij[3][6] = Mij[3][8] = Mij[3][9] =
      1.0 / 6;
  Mij[3][10] = Mij[3][13] = delRate / 2;

  genoArr[3][0] = 0;
  genoArr[3][1] = 3;
  genoArr[3][2] = INVALID;

  // 4, C/C
  Mij[4][1] = Mij[4][5] = Mij[4][6] = 1.0 / 3;
  Mij[4][11] = delRate;

  genoArr[4][0] = 1;
  genoArr[4][1] = 1;
  genoArr[4][2] = INVALID;

  // 5, C/G
  Mij[5][1] = Mij[5][2] = Mij[5][4] = Mij[5][6] = Mij[5][7] = Mij[5][8] =
      1.0 / 6;
  Mij[5][11] = Mij[5][12] = delRate / 2;

  genoArr[5][0] = 1;
  genoArr[5][1] = 2;
  genoArr[5][2] = INVALID;

  // 6, C/T
  Mij[6][1] = Mij[6][3] = Mij[6][4] = Mij[6][5] = Mij[6][8] = Mij[6][9] =
      1.0 / 6;
  Mij[6][11] = Mij[6][13] = delRate / 2;

  genoArr[6][0] = 1;
  genoArr[6][1] = 3;
  genoArr[6][2] = INVALID;

  // 7, G/G
  Mij[7][2] = Mij[7][5] = Mij[7][8] = 1.0 / 3;
  Mij[7][12] = delRate;

  genoArr[7][0] = 2;
  genoArr[7][1] = 2;
  genoArr[7][2] = INVALID;

  // 8, G/T
  Mij[8][2] = Mij[8][3] = Mij[8][5] = Mij[8][6] = Mij[8][7] = Mij[8][9] =
      1.0 / 6;
  Mij[8][12] = Mij[8][13] = delRate / 2;

  genoArr[8][0] = 2;
  genoArr[8][1] = 3;
  genoArr[8][2] = INVALID;

  // 9, T/T
  Mij[9][3] = Mij[9][6] = Mij[9][8] = 1.0 / 3;
  Mij[9][13] = delRate;

  genoArr[9][0] = 3;
  genoArr[9][1] = 3;
  genoArr[9][2] = INVALID;

  // 10, A/-
  Mij[10][11] = Mij[10][12] = Mij[10][13] = 1.0 / 6;
  Mij[10][14] = delRate / 2;

  genoArr[10][0] = 0;
  genoArr[10][1] = DELETED;
  genoArr[10][2] = INVALID;

  // 11, C/-
  Mij[11][10] = Mij[11][12] = Mij[11][13] = 1.0 / 6;
  Mij[11][14] = delRate / 2;

  genoArr[11][0] = 1;
  genoArr[11][1] = DELETED;
  genoArr[11][2] = INVALID;

  // 12, G/-
  Mij[12][10] = Mij[12][11] = Mij[12][13] = 1.0 / 6;
  Mij[12][14] = delRate / 2;

  genoArr[12][0] = 2;
  genoArr[12][1] = DELETED;
  genoArr[12][2] = INVALID;

  // 13, T/-
  Mij[13][10] = Mij[13][11] = Mij[13][12] = 1.0 / 6;
  Mij[13][14] = delRate / 2;

  genoArr[13][0] = 3;
  genoArr[13][1] = DELETED;
  genoArr[13][2] = INVALID;

  // 14, -
  genoArr[14][0] = genoArr[14][1] = DELETED;
  genoArr[14][2] = INVALID;

  hashMapIntArrToInt = CreateHashMapIntArrToInt(3);
  for (int i = 0; i < n; ++i)
    InsertKeyValIntArrToInt(hashMapIntArrToInt, genoArr[i], i);

  return hashMapIntArrToInt;
}

struct HashMapIntToIntArr *initHashMapIntToIntArrMDI() {
  int n = 65;
  double muRate = 1.0;
  double delRate = 1.0;
  double insRate = 1e-1;
  double Mij[n][n];

  if (genoArr != NULL) {
    FreeGenoArr(n);
  }
  AllocateGenoArr(n);

  // normalize the deletionRate and insRate regarding to the mutationRate
  delRate /= (muRate + 1e-15);
  insRate /= (muRate + 1e-15);

  // initialize all entries to 0.0
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      Mij[i][j] = 0.0;
    }

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 0, A/AA
  Mij[0][1] = Mij[0][2] = Mij[0][3] = 2 / 9.0;
  Mij[0][4] = Mij[0][5] = Mij[0][6] = 1 / 9.0;
  Mij[0][40] = 2 * delRate / 3.0;
  Mij[0][41] = delRate / 3.0;

  genoArr[0][0] = 0;
  genoArr[0][1] = 0;
  genoArr[0][2] = 0;

  // 1, A/AC
  Mij[1][0] = Mij[1][2] = Mij[1][3] = Mij[1][7] = Mij[1][8] = Mij[1][9] =
      Mij[1][13] = Mij[1][16] = Mij[1][19] = 1 / 9.0;
  Mij[1][40] = Mij[1][42] = Mij[1][45] = delRate / 3.0;

  genoArr[1][0] = 0;
  genoArr[1][1] = 0;
  genoArr[1][2] = 1;

  // 2, A/AG
  Mij[2][0] = Mij[2][1] = Mij[2][3] = Mij[2][8] = Mij[2][10] = Mij[2][11] =
      Mij[2][14] = Mij[2][17] = Mij[2][20] = 1 / 9.0;
  Mij[2][40] = Mij[2][43] = Mij[2][46] = delRate / 3.0;

  genoArr[2][0] = 0;
  genoArr[2][1] = 0;
  genoArr[2][2] = 2;

  // 3, A/AT
  Mij[3][0] = Mij[3][1] = Mij[3][2] = Mij[3][9] = Mij[3][11] = Mij[3][12] =
      Mij[3][15] = Mij[3][18] = Mij[3][21] = 1 / 9.0;
  Mij[3][40] = Mij[3][44] = Mij[3][47] = delRate / 3.0;

  genoArr[3][0] = 0;
  genoArr[3][1] = 0;
  genoArr[3][2] = 3;

  // 4, C/AA
  Mij[4][0] = Mij[4][5] = Mij[4][6] = 1 / 9.0;
  Mij[4][13] = Mij[4][14] = Mij[4][15] = 2 / 9.0;
  Mij[4][41] = delRate / 3.0;
  Mij[4][42] = 2 * delRate / 3.0;

  genoArr[4][0] = 1;
  genoArr[4][1] = 0;
  genoArr[4][2] = 0;

  // 5, G/AA
  Mij[5][0] = Mij[5][4] = Mij[5][6] = 1 / 9.0;
  Mij[5][16] = Mij[5][17] = Mij[5][18] = 2 / 9.0;
  Mij[5][41] = delRate / 3.0;
  Mij[5][43] = 2 * delRate / 3.0;

  genoArr[5][0] = 2;
  genoArr[5][1] = 0;
  genoArr[5][2] = 0;

  // 6, T/AA
  Mij[6][0] = Mij[6][4] = Mij[6][5] = 1 / 9.0;
  Mij[6][19] = Mij[6][20] = Mij[6][21] = 2 / 9.0;
  Mij[6][41] = delRate / 3.0;
  Mij[6][44] = 2 * delRate / 3.0;

  genoArr[6][0] = 3;
  genoArr[6][1] = 0;
  genoArr[6][2] = 0;

  // 7, A/CC
  Mij[7][1] = Mij[7][8] = Mij[7][9] = 2 / 9.0;
  Mij[7][22] = Mij[7][28] = Mij[7][34] = 1 / 9.0;
  Mij[7][42] = 2 * delRate / 3.0;
  Mij[7][54] = delRate / 3.0;

  genoArr[7][0] = 0;
  genoArr[7][1] = 1;
  genoArr[7][2] = 1;

  // 8, A/CG
  Mij[8][1] = Mij[8][2] = Mij[8][7] = Mij[8][9] = Mij[8][10] = Mij[8][11] =
      Mij[8][23] = Mij[8][29] = Mij[8][35] = 1 / 9.0;
  Mij[8][42] = Mij[8][43] = Mij[8][55] = delRate / 3.0;

  genoArr[8][0] = 0;
  genoArr[8][1] = 1;
  genoArr[8][2] = 2;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 9, A/CT
  Mij[9][1] = Mij[9][3] = Mij[9][7] = Mij[9][8] = Mij[9][11] = Mij[9][12] =
      Mij[9][24] = Mij[9][30] = Mij[9][36] = 1 / 9.0;
  Mij[9][42] = Mij[9][44] = Mij[9][56] = delRate / 3.0;

  genoArr[9][0] = 0;
  genoArr[9][1] = 1;
  genoArr[9][2] = 3;

  // 10, A/GG
  Mij[10][2] = Mij[10][8] = Mij[10][11] = 2 / 9.0;
  Mij[10][25] = Mij[10][31] = Mij[10][37] = 1 / 9.0;
  Mij[10][43] = 2 * delRate / 3.0;
  Mij[10][57] = delRate / 3.0;

  genoArr[10][0] = 0;
  genoArr[10][1] = 2;
  genoArr[10][2] = 2;

  // 11, A/GT
  Mij[11][2] = Mij[11][3] = Mij[11][8] = Mij[11][9] = Mij[11][10] =
      Mij[11][12] = Mij[11][26] = Mij[11][32] = Mij[11][38] = 1 / 9.0;
  Mij[11][43] = Mij[11][44] = Mij[11][58] = delRate / 3.0;

  genoArr[11][0] = 0;
  genoArr[11][1] = 2;
  genoArr[11][2] = 3;

  // 12, A/TT
  Mij[12][3] = Mij[12][9] = Mij[12][11] = 2 / 9.0;
  Mij[12][27] = Mij[12][33] = Mij[12][39] = 1 / 9.0;
  Mij[12][44] = 2 * delRate / 3.0;
  Mij[12][59] = delRate / 3.0;

  genoArr[12][0] = 0;
  genoArr[12][1] = 3;
  genoArr[12][2] = 3;

  // 13, C/AC
  Mij[13][1] = Mij[13][4] = Mij[13][14] = Mij[13][15] = Mij[13][16] =
      Mij[13][19] = Mij[13][22] = Mij[13][23] = Mij[13][24] = 1 / 9.0;
  Mij[13][42] = Mij[13][45] = Mij[13][48] = delRate / 3.0;

  genoArr[13][0] = 1;
  genoArr[13][1] = 0;
  genoArr[13][2] = 1;

  // 14, C/AG
  Mij[14][2] = Mij[14][4] = Mij[14][13] = Mij[14][15] = Mij[14][17] =
      Mij[14][20] = Mij[14][23] = Mij[14][25] = Mij[14][26] = 1 / 9.0;
  Mij[14][42] = Mij[14][46] = Mij[14][49] = delRate / 3.0;

  genoArr[14][0] = 1;
  genoArr[14][1] = 0;
  genoArr[14][2] = 2;

  // 15, C/AT
  Mij[15][3] = Mij[15][4] = Mij[15][13] = Mij[15][14] = Mij[15][18] =
      Mij[15][21] = Mij[15][24] = Mij[15][26] = Mij[15][27] = 1 / 9.0;
  Mij[15][42] = Mij[15][47] = Mij[15][50] = delRate / 3.0;

  genoArr[15][0] = 1;
  genoArr[15][1] = 0;
  genoArr[15][2] = 3;

  // 16, G/AC
  Mij[16][1] = Mij[16][5] = Mij[16][13] = Mij[16][17] = Mij[16][18] =
      Mij[16][19] = Mij[16][28] = Mij[16][29] = Mij[16][30] = 1 / 9.0;
  Mij[16][43] = Mij[16][45] = Mij[16][49] = delRate / 3.0;

  genoArr[16][0] = 2;
  genoArr[16][1] = 0;
  genoArr[16][2] = 1;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 17, G/AG
  Mij[17][2] = Mij[17][5] = Mij[17][14] = Mij[17][16] = Mij[17][18] =
      Mij[17][20] = Mij[17][29] = Mij[17][31] = Mij[17][32] = 1 / 9.0;
  Mij[17][43] = Mij[17][46] = Mij[17][51] = delRate / 3.0;

  genoArr[17][0] = 2;
  genoArr[17][1] = 0;
  genoArr[17][2] = 2;

  // 18, G/AT
  Mij[18][3] = Mij[18][5] = Mij[18][15] = Mij[18][16] = Mij[18][17] =
      Mij[18][21] = Mij[18][30] = Mij[18][32] = Mij[18][33] = 1 / 9.0;
  Mij[18][43] = Mij[18][47] = Mij[18][52] = delRate / 3.0;

  genoArr[18][0] = 2;
  genoArr[18][1] = 0;
  genoArr[18][2] = 3;

  // 19, T/AC
  Mij[19][1] = Mij[19][6] = Mij[19][13] = Mij[19][16] = Mij[19][20] =
      Mij[19][21] = Mij[19][34] = Mij[19][35] = Mij[19][36] = 1 / 9.0;
  Mij[19][44] = Mij[19][45] = Mij[19][50] = delRate / 3.0;

  genoArr[19][0] = 3;
  genoArr[19][1] = 0;
  genoArr[19][2] = 1;

  // 20, T/AG
  Mij[20][2] = Mij[20][6] = Mij[20][14] = Mij[20][17] = Mij[20][19] =
      Mij[20][21] = Mij[20][35] = Mij[20][37] = Mij[20][38] = 1 / 9.0;
  Mij[20][44] = Mij[20][46] = Mij[20][52] = delRate / 3.0;

  genoArr[20][0] = 3;
  genoArr[20][1] = 0;
  genoArr[20][2] = 2;

  // 21, T/AT
  Mij[21][3] = Mij[21][6] = Mij[21][15] = Mij[21][18] = Mij[21][19] =
      Mij[21][20] = Mij[21][36] = Mij[21][38] = Mij[21][39] = 1 / 9.0;
  Mij[21][44] = Mij[21][47] = Mij[21][53] = delRate / 3.0;

  genoArr[21][0] = 3;
  genoArr[21][1] = 0;
  genoArr[21][2] = 3;

  // 22, C/CC
  Mij[22][7] = Mij[22][28] = Mij[22][34] = 1 / 9.0;
  Mij[22][13] = Mij[22][23] = Mij[22][24] = 2 / 9.0;
  Mij[22][48] = 2 * delRate / 3.0;
  Mij[22][54] = delRate / 3.0;

  genoArr[22][0] = 1;
  genoArr[22][1] = 1;
  genoArr[22][2] = 1;

  // 23, C/CG
  Mij[23][8] = Mij[23][13] = Mij[23][14] = Mij[23][22] = Mij[23][24] =
      Mij[23][25] = Mij[23][26] = Mij[23][29] = Mij[23][35] = 1 / 9.0;
  Mij[23][48] = Mij[23][49] = Mij[23][55] = delRate / 3.0;

  genoArr[23][0] = 1;
  genoArr[23][1] = 1;
  genoArr[23][2] = 2;

  // 24, C/CT
  Mij[24][9] = Mij[24][13] = Mij[24][15] = Mij[24][22] = Mij[24][23] =
      Mij[24][26] = Mij[24][27] = Mij[24][30] = Mij[24][36] = 1 / 9.0;
  Mij[24][48] = Mij[24][50] = Mij[24][56] = delRate / 3.0;

  genoArr[24][0] = 1;
  genoArr[24][1] = 1;
  genoArr[24][2] = 3;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 25, C/GG
  Mij[25][10] = Mij[25][31] = Mij[25][37] = 1 / 9.0;
  Mij[25][14] = Mij[25][23] = Mij[25][26] = 2 / 9.0;
  Mij[25][49] = 2 * delRate / 3.0;
  Mij[25][57] = delRate / 3.0;

  genoArr[25][0] = 1;
  genoArr[25][1] = 2;
  genoArr[25][2] = 2;

  // 26, C/GT
  Mij[26][11] = Mij[26][14] = Mij[26][15] = Mij[26][23] = Mij[26][24] =
      Mij[26][25] = Mij[26][27] = Mij[26][32] = Mij[26][38] = 1 / 9.0;
  Mij[26][49] = Mij[26][50] = Mij[26][58] = delRate / 3.0;

  genoArr[26][0] = 1;
  genoArr[26][1] = 2;
  genoArr[26][2] = 3;

  // 27, C/TT
  Mij[27][12] = Mij[27][33] = Mij[27][39] = 1 / 9.0;
  Mij[27][15] = Mij[27][24] = Mij[27][26] = 2 / 9.0;
  Mij[27][50] = 2 * delRate / 3.0;
  Mij[27][59] = delRate / 3.0;

  genoArr[27][0] = 1;
  genoArr[27][1] = 3;
  genoArr[27][2] = 3;

  // 28, G/CC
  Mij[28][7] = Mij[28][22] = Mij[28][34] = 1 / 9.0;
  Mij[28][16] = Mij[28][29] = Mij[28][30] = 2 / 9.0;
  Mij[28][49] = 2 * delRate / 3.0;
  Mij[28][54] = delRate / 3.0;

  genoArr[28][0] = 2;
  genoArr[28][1] = 1;
  genoArr[28][2] = 1;

  // 29, G/CG
  Mij[29][8] = Mij[29][16] = Mij[29][17] = Mij[29][23] = Mij[29][28] =
      Mij[29][30] = Mij[29][31] = Mij[29][32] = Mij[29][35] = 1 / 9.0;
  Mij[29][49] = Mij[29][51] = Mij[29][55] = delRate / 3.0;

  genoArr[29][0] = 2;
  genoArr[29][1] = 1;
  genoArr[29][2] = 2;

  // 30, G/CT
  Mij[30][9] = Mij[30][16] = Mij[30][18] = Mij[30][24] = Mij[30][28] =
      Mij[30][29] = Mij[30][32] = Mij[30][33] = Mij[30][36] = 1 / 9.0;
  Mij[30][49] = Mij[30][52] = Mij[30][56] = delRate / 3.0;

  genoArr[30][0] = 2;
  genoArr[30][1] = 1;
  genoArr[30][2] = 3;

  // 31, G/GG
  Mij[31][10] = Mij[31][25] = Mij[31][37] = 1 / 9.0;
  Mij[31][17] = Mij[31][29] = Mij[31][32] = 2 / 9.0;
  Mij[31][51] = 2 * delRate / 3.0;
  Mij[31][57] = delRate / 3.0;

  genoArr[31][0] = 2;
  genoArr[31][1] = 2;
  genoArr[31][2] = 2;

  // 32, G/GT
  Mij[32][11] = Mij[32][17] = Mij[32][18] = Mij[32][26] = Mij[32][29] =
      Mij[32][30] = Mij[32][31] = Mij[32][33] = Mij[32][38] = 1 / 9.0;
  Mij[32][51] = Mij[32][52] = Mij[32][58] = delRate / 3.0;

  genoArr[32][0] = 2;
  genoArr[32][1] = 2;
  genoArr[32][2] = 3;

  // 33, G/TT
  Mij[33][12] = Mij[33][27] = Mij[33][39] = 1 / 9.0;
  Mij[33][18] = Mij[33][30] = Mij[33][32] = 2 / 9.0;
  Mij[33][52] = 2 * delRate / 3.0;
  Mij[33][59] = delRate / 3.0;

  genoArr[33][0] = 2;
  genoArr[33][1] = 3;
  genoArr[33][2] = 3;

  // 34, T/CC
  Mij[34][7] = Mij[34][22] = Mij[34][28] = 1 / 9.0;
  Mij[34][19] = Mij[34][35] = Mij[34][36] = 2 / 9.0;
  Mij[34][50] = 2 * delRate / 3.0;
  Mij[34][54] = delRate / 3.0;

  genoArr[34][0] = 3;
  genoArr[34][1] = 1;
  genoArr[34][2] = 1;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 35, T/CG
  Mij[35][8] = Mij[35][19] = Mij[35][20] = Mij[35][23] = Mij[35][29] =
      Mij[35][34] = Mij[35][36] = Mij[35][37] = Mij[35][38] = 1 / 9.0;
  Mij[35][50] = Mij[35][52] = Mij[35][55] = delRate / 3.0;

  genoArr[35][0] = 3;
  genoArr[35][1] = 1;
  genoArr[35][2] = 2;

  // 36, T/CT
  Mij[36][9] = Mij[36][19] = Mij[36][21] = Mij[36][24] = Mij[36][30] =
      Mij[36][34] = Mij[36][35] = Mij[36][38] = Mij[36][39] = 1 / 9.0;
  Mij[36][50] = Mij[36][53] = Mij[36][56] = delRate / 3.0;

  genoArr[36][0] = 3;
  genoArr[36][1] = 1;
  genoArr[36][2] = 3;

  // 37, T/GG
  Mij[37][10] = Mij[37][25] = Mij[37][31] = 1 / 9.0;
  Mij[37][20] = Mij[37][35] = Mij[37][38] = 2 / 9.0;
  Mij[37][52] = 2 * delRate / 3.0;
  Mij[37][57] = delRate / 3.0;

  genoArr[37][0] = 3;
  genoArr[37][1] = 2;
  genoArr[37][2] = 2;

  // 38, T/GT
  Mij[38][11] = Mij[38][20] = Mij[38][21] = Mij[38][26] = Mij[38][32] =
      Mij[38][35] = Mij[38][36] = Mij[38][37] = Mij[38][39] = 1 / 9.0;
  Mij[38][52] = Mij[38][53] = Mij[38][58] = delRate / 3.0;

  genoArr[38][0] = 3;
  genoArr[38][1] = 2;
  genoArr[38][2] = 3;

  // 39, T/TT
  Mij[39][12] = Mij[39][27] = Mij[39][33] = 1 / 9.0;
  Mij[39][21] = Mij[39][36] = Mij[39][38] = 2 / 9.0;
  Mij[39][53] = 2 * delRate / 3.0;
  Mij[39][59] = delRate / 3.0;

  genoArr[39][0] = 3;
  genoArr[39][1] = 3;
  genoArr[39][2] = 3;

  // 40, A/A
  Mij[40][0] = 2 * insRate / 3.0;
  Mij[40][42] = Mij[40][43] = Mij[40][44] = 2 / 9.0;
  Mij[40][60] = 2 * delRate / 3.0;

  genoArr[40][0] = 0;
  genoArr[40][1] = 0;
  genoArr[40][2] = INVALID;

  // 41, -/AA
  Mij[41][45] = Mij[41][46] = Mij[41][47] = 2 / 9.0;
  Mij[41][60] = 2 * delRate / 3.0;

  genoArr[41][0] = DELETED;
  genoArr[41][1] = 0;
  genoArr[41][2] = 0;

  // 42, A/C
  Mij[42][4] = Mij[42][7] = insRate / 3.0;
  Mij[42][40] = Mij[42][43] = Mij[42][44] = Mij[42][48] = Mij[42][49] =
      Mij[42][50] = 1 / 9.0;
  Mij[42][60] = Mij[42][61] = delRate / 3.0;

  genoArr[42][0] = 0;
  genoArr[42][1] = 1;
  genoArr[42][2] = INVALID;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 43, A/G
  Mij[43][5] = Mij[43][10] = insRate / 3.0;
  Mij[43][40] = Mij[43][42] = Mij[43][44] = Mij[43][49] = Mij[43][51] =
      Mij[43][52] = 1 / 9.0;
  Mij[43][60] = Mij[43][62] = delRate / 3.0;

  genoArr[43][0] = 0;
  genoArr[43][1] = 2;
  genoArr[43][2] = INVALID;

  // 44, A/T
  Mij[44][6] = Mij[44][12] = insRate / 3.0;
  Mij[44][40] = Mij[44][42] = Mij[44][43] = Mij[44][50] = Mij[44][52] =
      Mij[44][53] = 1 / 9.0;
  Mij[44][60] = Mij[44][63] = delRate / 3.0;

  genoArr[44][0] = 0;
  genoArr[44][1] = 3;
  genoArr[44][2] = INVALID;

  // 45, -/AC
  Mij[45][41] = Mij[45][46] = Mij[45][47] = Mij[45][54] = Mij[45][55] =
      Mij[45][56] = 1 / 9.0;
  Mij[45][60] = Mij[45][61] = delRate / 3.0;

  genoArr[45][0] = DELETED;
  genoArr[45][1] = 0;
  genoArr[45][2] = 1;

  // 46, -/AG
  Mij[46][41] = Mij[46][45] = Mij[46][47] = Mij[46][55] = Mij[46][57] =
      Mij[46][58] = 1 / 9.0;
  Mij[46][60] = Mij[46][62] = delRate / 3.0;

  genoArr[46][0] = DELETED;
  genoArr[46][1] = 0;
  genoArr[46][2] = 2;

  // 47, -/AT
  Mij[47][41] = Mij[47][45] = Mij[47][46] = Mij[47][56] = Mij[47][58] =
      Mij[47][59] = 1 / 9.0;
  Mij[47][60] = Mij[47][63] = delRate / 3.0;

  genoArr[47][0] = DELETED;
  genoArr[47][1] = 0;
  genoArr[47][2] = 3;

  // 48, C/C
  Mij[48][22] = 2 * insRate / 3.0;
  Mij[48][42] = Mij[48][49] = Mij[48][50] = 2 / 9.0;
  Mij[48][61] = 2 * delRate / 3.0;

  genoArr[48][0] = 1;
  genoArr[48][1] = 1;
  genoArr[48][2] = INVALID;

  // 49, C/G
  Mij[49][25] = Mij[49][28] = insRate / 3.0;
  Mij[49][42] = Mij[49][43] = Mij[49][48] = Mij[49][50] = Mij[49][51] =
      Mij[49][52] = 1 / 9.0;
  Mij[49][61] = Mij[49][62] = delRate / 3.0;

  genoArr[49][0] = 1;
  genoArr[49][1] = 2;
  genoArr[49][2] = INVALID;

  // 50, C/T
  Mij[50][27] = Mij[50][34] = insRate / 3.0;
  Mij[50][42] = Mij[50][44] = Mij[50][48] = Mij[50][49] = Mij[50][52] =
      Mij[50][53] = 1 / 9.0;
  Mij[50][61] = Mij[50][63] = delRate / 3.0;

  genoArr[50][0] = 1;
  genoArr[50][1] = 3;
  genoArr[50][2] = INVALID;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 51, G/G
  Mij[51][31] = 2 * insRate / 3.0;
  Mij[51][43] = Mij[51][49] = Mij[51][52] = 2 / 9.0;
  Mij[51][62] = delRate / 3.0;

  genoArr[51][0] = 2;
  genoArr[51][1] = 2;
  genoArr[51][2] = INVALID;

  // 52, G/T
  Mij[52][33] = Mij[52][37] = insRate / 3.0;
  Mij[52][43] = Mij[52][44] = Mij[52][49] = Mij[52][50] = Mij[52][51] =
      Mij[52][53] = 1 / 9.0;
  Mij[52][62] = Mij[52][63] = delRate / 3.0;

  genoArr[52][0] = 2;
  genoArr[52][1] = 3;
  genoArr[52][2] = INVALID;

  // 53, T/T
  Mij[53][39] = 2 * insRate / 3.0;
  Mij[53][44] = Mij[53][50] = Mij[53][52] = 2 / 9.0;
  Mij[53][63] = 2 * delRate / 3.0;

  genoArr[53][0] = 3;
  genoArr[53][1] = 3;
  genoArr[53][2] = INVALID;

  // 54, -/CC
  Mij[54][45] = Mij[54][55] = Mij[54][56] = 2 / 9.0;
  Mij[54][61] = 2 * delRate / 3.0;

  genoArr[54][0] = DELETED;
  genoArr[54][1] = 1;
  genoArr[54][2] = 1;

  // 55, -/CG
  Mij[55][45] = Mij[55][46] = Mij[55][54] = Mij[55][56] = Mij[55][57] =
      Mij[55][58] = 1 / 9.0;
  Mij[55][61] = Mij[55][62] = delRate / 3.0;

  genoArr[55][0] = DELETED;
  genoArr[55][1] = 1;
  genoArr[55][2] = 2;

  // 56, -/CT
  Mij[56][45] = Mij[56][47] = Mij[56][54] = Mij[56][55] = Mij[56][58] =
      Mij[56][59] = 1 / 9.0;
  Mij[56][61] = Mij[56][63] = delRate / 3.0;

  genoArr[56][0] = DELETED;
  genoArr[56][1] = 1;
  genoArr[56][2] = 3;

  // 57, -/GG
  Mij[57][46] = Mij[57][55] = Mij[57][58] = 2 / 9.0;
  Mij[57][62] = 2 * delRate / 3.0;

  genoArr[57][0] = DELETED;
  genoArr[57][1] = 2;
  genoArr[57][2] = 2;

  // 58, -/GT
  Mij[58][46] = Mij[58][47] = Mij[58][55] = Mij[58][56] = Mij[58][57] =
      Mij[58][59] = 1 / 9.0;
  Mij[58][62] = Mij[58][63] = delRate / 3.0;

  genoArr[58][0] = DELETED;
  genoArr[58][1] = 2;
  genoArr[58][2] = 3;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 59, -/TT
  Mij[59][47] = Mij[59][56] = Mij[59][58] = 2 / 9.0;
  Mij[59][63] = 2 * delRate / 3.0;

  genoArr[59][0] = DELETED;
  genoArr[59][1] = 3;
  genoArr[59][2] = 3;

  // 60, A/-
  Mij[60][41] = insRate / 3.0;
  Mij[60][61] = Mij[60][62] = Mij[60][63] = 1 / 9.0;
  Mij[60][64] = delRate / 3.0;

  genoArr[60][0] = 0;
  genoArr[60][1] = DELETED;
  genoArr[60][2] = INVALID;

  // 61, C/-
  Mij[61][54] = insRate / 3.0;
  Mij[61][60] = Mij[61][62] = Mij[61][63] = 1 / 9.0;
  Mij[61][64] = delRate / 3.0;

  genoArr[61][0] = 1;
  genoArr[61][1] = DELETED;
  genoArr[61][2] = INVALID;

  // 62, G/-
  Mij[62][57] = insRate / 3.0;
  Mij[62][60] = Mij[62][61] = Mij[62][63] = 1 / 9.0;
  Mij[62][64] = delRate / 3.0;

  genoArr[62][0] = 2;
  genoArr[62][1] = DELETED;
  genoArr[62][2] = INVALID;

  // 63, T/-
  Mij[63][59] = insRate / 3.0;
  Mij[63][60] = Mij[63][61] = Mij[63][62] = 1 / 9.0;
  Mij[63][64] = delRate / 3.0;

  genoArr[63][0] = 3;
  genoArr[63][1] = DELETED;
  genoArr[63][2] = INVALID;

  // 64, -
  genoArr[64][0] = genoArr[64][1] = DELETED;
  genoArr[64][2] = INVALID;

  hashMapIntToIntArr = CreateHashMapIntToIntArr(13);
  for (int i = 0; i < n; ++i)
    InsertKeyValIntToIntArr(hashMapIntToIntArr, i, genoArr[i]);

  return hashMapIntToIntArr;
}

struct HashMapIntArrToInt *initHashMapIntArrToIntMDI() {
  int n = 65;
  double muRate = 1.0;
  double delRate = 1.0;
  double insRate = 1e-1;
  double Mij[n][n];

  if (genoArr != NULL) {
    FreeGenoArr(n);
  }
  AllocateGenoArr(n);

  // normalize the deletionRate and insRate regarding to the mutationRate
  delRate /= (muRate + 1e-15);
  insRate /= (muRate + 1e-15);

  // initialize all entries to 0.0
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      Mij[i][j] = 0.0;
    }

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 0, A/AA
  Mij[0][1] = Mij[0][2] = Mij[0][3] = 2 / 9.0;
  Mij[0][4] = Mij[0][5] = Mij[0][6] = 1 / 9.0;
  Mij[0][40] = 2 * delRate / 3.0;
  Mij[0][41] = delRate / 3.0;

  genoArr[0][0] = 0;
  genoArr[0][1] = 0;
  genoArr[0][2] = 0;

  // 1, A/AC
  Mij[1][0] = Mij[1][2] = Mij[1][3] = Mij[1][7] = Mij[1][8] = Mij[1][9] =
      Mij[1][13] = Mij[1][16] = Mij[1][19] = 1 / 9.0;
  Mij[1][40] = Mij[1][42] = Mij[1][45] = delRate / 3.0;

  genoArr[1][0] = 0;
  genoArr[1][1] = 0;
  genoArr[1][2] = 1;

  // 2, A/AG
  Mij[2][0] = Mij[2][1] = Mij[2][3] = Mij[2][8] = Mij[2][10] = Mij[2][11] =
      Mij[2][14] = Mij[2][17] = Mij[2][20] = 1 / 9.0;
  Mij[2][40] = Mij[2][43] = Mij[2][46] = delRate / 3.0;

  genoArr[2][0] = 0;
  genoArr[2][1] = 0;
  genoArr[2][2] = 2;

  // 3, A/AT
  Mij[3][0] = Mij[3][1] = Mij[3][2] = Mij[3][9] = Mij[3][11] = Mij[3][12] =
      Mij[3][15] = Mij[3][18] = Mij[3][21] = 1 / 9.0;
  Mij[3][40] = Mij[3][44] = Mij[3][47] = delRate / 3.0;

  genoArr[3][0] = 0;
  genoArr[3][1] = 0;
  genoArr[3][2] = 3;

  // 4, C/AA
  Mij[4][0] = Mij[4][5] = Mij[4][6] = 1 / 9.0;
  Mij[4][13] = Mij[4][14] = Mij[4][15] = 2 / 9.0;
  Mij[4][41] = delRate / 3.0;
  Mij[4][42] = 2 * delRate / 3.0;

  genoArr[4][0] = 1;
  genoArr[4][1] = 0;
  genoArr[4][2] = 0;

  // 5, G/AA
  Mij[5][0] = Mij[5][4] = Mij[5][6] = 1 / 9.0;
  Mij[5][16] = Mij[5][17] = Mij[5][18] = 2 / 9.0;
  Mij[5][41] = delRate / 3.0;
  Mij[5][43] = 2 * delRate / 3.0;

  genoArr[5][0] = 2;
  genoArr[5][1] = 0;
  genoArr[5][2] = 0;

  // 6, T/AA
  Mij[6][0] = Mij[6][4] = Mij[6][5] = 1 / 9.0;
  Mij[6][19] = Mij[6][20] = Mij[6][21] = 2 / 9.0;
  Mij[6][41] = delRate / 3.0;
  Mij[6][44] = 2 * delRate / 3.0;

  genoArr[6][0] = 3;
  genoArr[6][1] = 0;
  genoArr[6][2] = 0;

  // 7, A/CC
  Mij[7][1] = Mij[7][8] = Mij[7][9] = 2 / 9.0;
  Mij[7][22] = Mij[7][28] = Mij[7][34] = 1 / 9.0;
  Mij[7][42] = 2 * delRate / 3.0;
  Mij[7][54] = delRate / 3.0;

  genoArr[7][0] = 0;
  genoArr[7][1] = 1;
  genoArr[7][2] = 1;

  // 8, A/CG
  Mij[8][1] = Mij[8][2] = Mij[8][7] = Mij[8][9] = Mij[8][10] = Mij[8][11] =
      Mij[8][23] = Mij[8][29] = Mij[8][35] = 1 / 9.0;
  Mij[8][42] = Mij[8][43] = Mij[8][55] = delRate / 3.0;

  genoArr[8][0] = 0;
  genoArr[8][1] = 1;
  genoArr[8][2] = 2;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 9, A/CT
  Mij[9][1] = Mij[9][3] = Mij[9][7] = Mij[9][8] = Mij[9][11] = Mij[9][12] =
      Mij[9][24] = Mij[9][30] = Mij[9][36] = 1 / 9.0;
  Mij[9][42] = Mij[9][44] = Mij[9][56] = delRate / 3.0;

  genoArr[9][0] = 0;
  genoArr[9][1] = 1;
  genoArr[9][2] = 3;

  // 10, A/GG
  Mij[10][2] = Mij[10][8] = Mij[10][11] = 2 / 9.0;
  Mij[10][25] = Mij[10][31] = Mij[10][37] = 1 / 9.0;
  Mij[10][43] = 2 * delRate / 3.0;
  Mij[10][57] = delRate / 3.0;

  genoArr[10][0] = 0;
  genoArr[10][1] = 2;
  genoArr[10][2] = 2;

  // 11, A/GT
  Mij[11][2] = Mij[11][3] = Mij[11][8] = Mij[11][9] = Mij[11][10] =
      Mij[11][12] = Mij[11][26] = Mij[11][32] = Mij[11][38] = 1 / 9.0;
  Mij[11][43] = Mij[11][44] = Mij[11][58] = delRate / 3.0;

  genoArr[11][0] = 0;
  genoArr[11][1] = 2;
  genoArr[11][2] = 3;

  // 12, A/TT
  Mij[12][3] = Mij[12][9] = Mij[12][11] = 2 / 9.0;
  Mij[12][27] = Mij[12][33] = Mij[12][39] = 1 / 9.0;
  Mij[12][44] = 2 * delRate / 3.0;
  Mij[12][59] = delRate / 3.0;

  genoArr[12][0] = 0;
  genoArr[12][1] = 3;
  genoArr[12][2] = 3;

  // 13, C/AC
  Mij[13][1] = Mij[13][4] = Mij[13][14] = Mij[13][15] = Mij[13][16] =
      Mij[13][19] = Mij[13][22] = Mij[13][23] = Mij[13][24] = 1 / 9.0;
  Mij[13][42] = Mij[13][45] = Mij[13][48] = delRate / 3.0;

  genoArr[13][0] = 1;
  genoArr[13][1] = 0;
  genoArr[13][2] = 1;

  // 14, C/AG
  Mij[14][2] = Mij[14][4] = Mij[14][13] = Mij[14][15] = Mij[14][17] =
      Mij[14][20] = Mij[14][23] = Mij[14][25] = Mij[14][26] = 1 / 9.0;
  Mij[14][42] = Mij[14][46] = Mij[14][49] = delRate / 3.0;

  genoArr[14][0] = 1;
  genoArr[14][1] = 0;
  genoArr[14][2] = 2;

  // 15, C/AT
  Mij[15][3] = Mij[15][4] = Mij[15][13] = Mij[15][14] = Mij[15][18] =
      Mij[15][21] = Mij[15][24] = Mij[15][26] = Mij[15][27] = 1 / 9.0;
  Mij[15][42] = Mij[15][47] = Mij[15][50] = delRate / 3.0;

  genoArr[15][0] = 1;
  genoArr[15][1] = 0;
  genoArr[15][2] = 3;

  // 16, G/AC
  Mij[16][1] = Mij[16][5] = Mij[16][13] = Mij[16][17] = Mij[16][18] =
      Mij[16][19] = Mij[16][28] = Mij[16][29] = Mij[16][30] = 1 / 9.0;
  Mij[16][43] = Mij[16][45] = Mij[16][49] = delRate / 3.0;

  genoArr[16][0] = 2;
  genoArr[16][1] = 0;
  genoArr[16][2] = 1;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 17, G/AG
  Mij[17][2] = Mij[17][5] = Mij[17][14] = Mij[17][16] = Mij[17][18] =
      Mij[17][20] = Mij[17][29] = Mij[17][31] = Mij[17][32] = 1 / 9.0;
  Mij[17][43] = Mij[17][46] = Mij[17][51] = delRate / 3.0;

  genoArr[17][0] = 2;
  genoArr[17][1] = 0;
  genoArr[17][2] = 2;

  // 18, G/AT
  Mij[18][3] = Mij[18][5] = Mij[18][15] = Mij[18][16] = Mij[18][17] =
      Mij[18][21] = Mij[18][30] = Mij[18][32] = Mij[18][33] = 1 / 9.0;
  Mij[18][43] = Mij[18][47] = Mij[18][52] = delRate / 3.0;

  genoArr[18][0] = 2;
  genoArr[18][1] = 0;
  genoArr[18][2] = 3;

  // 19, T/AC
  Mij[19][1] = Mij[19][6] = Mij[19][13] = Mij[19][16] = Mij[19][20] =
      Mij[19][21] = Mij[19][34] = Mij[19][35] = Mij[19][36] = 1 / 9.0;
  Mij[19][44] = Mij[19][45] = Mij[19][50] = delRate / 3.0;

  genoArr[19][0] = 3;
  genoArr[19][1] = 0;
  genoArr[19][2] = 1;

  // 20, T/AG
  Mij[20][2] = Mij[20][6] = Mij[20][14] = Mij[20][17] = Mij[20][19] =
      Mij[20][21] = Mij[20][35] = Mij[20][37] = Mij[20][38] = 1 / 9.0;
  Mij[20][44] = Mij[20][46] = Mij[20][52] = delRate / 3.0;

  genoArr[20][0] = 3;
  genoArr[20][1] = 0;
  genoArr[20][2] = 2;

  // 21, T/AT
  Mij[21][3] = Mij[21][6] = Mij[21][15] = Mij[21][18] = Mij[21][19] =
      Mij[21][20] = Mij[21][36] = Mij[21][38] = Mij[21][39] = 1 / 9.0;
  Mij[21][44] = Mij[21][47] = Mij[21][53] = delRate / 3.0;

  genoArr[21][0] = 3;
  genoArr[21][1] = 0;
  genoArr[21][2] = 3;

  // 22, C/CC
  Mij[22][7] = Mij[22][28] = Mij[22][34] = 1 / 9.0;
  Mij[22][13] = Mij[22][23] = Mij[22][24] = 2 / 9.0;
  Mij[22][48] = 2 * delRate / 3.0;
  Mij[22][54] = delRate / 3.0;

  genoArr[22][0] = 1;
  genoArr[22][1] = 1;
  genoArr[22][2] = 1;

  // 23, C/CG
  Mij[23][8] = Mij[23][13] = Mij[23][14] = Mij[23][22] = Mij[23][24] =
      Mij[23][25] = Mij[23][26] = Mij[23][29] = Mij[23][35] = 1 / 9.0;
  Mij[23][48] = Mij[23][49] = Mij[23][55] = delRate / 3.0;

  genoArr[23][0] = 1;
  genoArr[23][1] = 1;
  genoArr[23][2] = 2;

  // 24, C/CT
  Mij[24][9] = Mij[24][13] = Mij[24][15] = Mij[24][22] = Mij[24][23] =
      Mij[24][26] = Mij[24][27] = Mij[24][30] = Mij[24][36] = 1 / 9.0;
  Mij[24][48] = Mij[24][50] = Mij[24][56] = delRate / 3.0;

  genoArr[24][0] = 1;
  genoArr[24][1] = 1;
  genoArr[24][2] = 3;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 25, C/GG
  Mij[25][10] = Mij[25][31] = Mij[25][37] = 1 / 9.0;
  Mij[25][14] = Mij[25][23] = Mij[25][26] = 2 / 9.0;
  Mij[25][49] = 2 * delRate / 3.0;
  Mij[25][57] = delRate / 3.0;

  genoArr[25][0] = 1;
  genoArr[25][1] = 2;
  genoArr[25][2] = 2;

  // 26, C/GT
  Mij[26][11] = Mij[26][14] = Mij[26][15] = Mij[26][23] = Mij[26][24] =
      Mij[26][25] = Mij[26][27] = Mij[26][32] = Mij[26][38] = 1 / 9.0;
  Mij[26][49] = Mij[26][50] = Mij[26][58] = delRate / 3.0;

  genoArr[26][0] = 1;
  genoArr[26][1] = 2;
  genoArr[26][2] = 3;

  // 27, C/TT
  Mij[27][12] = Mij[27][33] = Mij[27][39] = 1 / 9.0;
  Mij[27][15] = Mij[27][24] = Mij[27][26] = 2 / 9.0;
  Mij[27][50] = 2 * delRate / 3.0;
  Mij[27][59] = delRate / 3.0;

  genoArr[27][0] = 1;
  genoArr[27][1] = 3;
  genoArr[27][2] = 3;

  // 28, G/CC
  Mij[28][7] = Mij[28][22] = Mij[28][34] = 1 / 9.0;
  Mij[28][16] = Mij[28][29] = Mij[28][30] = 2 / 9.0;
  Mij[28][49] = 2 * delRate / 3.0;
  Mij[28][54] = delRate / 3.0;

  genoArr[28][0] = 2;
  genoArr[28][1] = 1;
  genoArr[28][2] = 1;

  // 29, G/CG
  Mij[29][8] = Mij[29][16] = Mij[29][17] = Mij[29][23] = Mij[29][28] =
      Mij[29][30] = Mij[29][31] = Mij[29][32] = Mij[29][35] = 1 / 9.0;
  Mij[29][49] = Mij[29][51] = Mij[29][55] = delRate / 3.0;

  genoArr[29][0] = 2;
  genoArr[29][1] = 1;
  genoArr[29][2] = 2;

  // 30, G/CT
  Mij[30][9] = Mij[30][16] = Mij[30][18] = Mij[30][24] = Mij[30][28] =
      Mij[30][29] = Mij[30][32] = Mij[30][33] = Mij[30][36] = 1 / 9.0;
  Mij[30][49] = Mij[30][52] = Mij[30][56] = delRate / 3.0;

  genoArr[30][0] = 2;
  genoArr[30][1] = 1;
  genoArr[30][2] = 3;

  // 31, G/GG
  Mij[31][10] = Mij[31][25] = Mij[31][37] = 1 / 9.0;
  Mij[31][17] = Mij[31][29] = Mij[31][32] = 2 / 9.0;
  Mij[31][51] = 2 * delRate / 3.0;
  Mij[31][57] = delRate / 3.0;

  genoArr[31][0] = 2;
  genoArr[31][1] = 2;
  genoArr[31][2] = 2;

  // 32, G/GT
  Mij[32][11] = Mij[32][17] = Mij[32][18] = Mij[32][26] = Mij[32][29] =
      Mij[32][30] = Mij[32][31] = Mij[32][33] = Mij[32][38] = 1 / 9.0;
  Mij[32][51] = Mij[32][52] = Mij[32][58] = delRate / 3.0;

  genoArr[32][0] = 2;
  genoArr[32][1] = 2;
  genoArr[32][2] = 3;

  // 33, G/TT
  Mij[33][12] = Mij[33][27] = Mij[33][39] = 1 / 9.0;
  Mij[33][18] = Mij[33][30] = Mij[33][32] = 2 / 9.0;
  Mij[33][52] = 2 * delRate / 3.0;
  Mij[33][59] = delRate / 3.0;

  genoArr[33][0] = 2;
  genoArr[33][1] = 3;
  genoArr[33][2] = 3;

  // 34, T/CC
  Mij[34][7] = Mij[34][22] = Mij[34][28] = 1 / 9.0;
  Mij[34][19] = Mij[34][35] = Mij[34][36] = 2 / 9.0;
  Mij[34][50] = 2 * delRate / 3.0;
  Mij[34][54] = delRate / 3.0;

  genoArr[34][0] = 3;
  genoArr[34][1] = 1;
  genoArr[34][2] = 1;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 35, T/CG
  Mij[35][8] = Mij[35][19] = Mij[35][20] = Mij[35][23] = Mij[35][29] =
      Mij[35][34] = Mij[35][36] = Mij[35][37] = Mij[35][38] = 1 / 9.0;
  Mij[35][50] = Mij[35][52] = Mij[35][55] = delRate / 3.0;

  genoArr[35][0] = 3;
  genoArr[35][1] = 1;
  genoArr[35][2] = 2;

  // 36, T/CT
  Mij[36][9] = Mij[36][19] = Mij[36][21] = Mij[36][24] = Mij[36][30] =
      Mij[36][34] = Mij[36][35] = Mij[36][38] = Mij[36][39] = 1 / 9.0;
  Mij[36][50] = Mij[36][53] = Mij[36][56] = delRate / 3.0;

  genoArr[36][0] = 3;
  genoArr[36][1] = 1;
  genoArr[36][2] = 3;

  // 37, T/GG
  Mij[37][10] = Mij[37][25] = Mij[37][31] = 1 / 9.0;
  Mij[37][20] = Mij[37][35] = Mij[37][38] = 2 / 9.0;
  Mij[37][52] = 2 * delRate / 3.0;
  Mij[37][57] = delRate / 3.0;

  genoArr[37][0] = 3;
  genoArr[37][1] = 2;
  genoArr[37][2] = 2;

  // 38, T/GT
  Mij[38][11] = Mij[38][20] = Mij[38][21] = Mij[38][26] = Mij[38][32] =
      Mij[38][35] = Mij[38][36] = Mij[38][37] = Mij[38][39] = 1 / 9.0;
  Mij[38][52] = Mij[38][53] = Mij[38][58] = delRate / 3.0;

  genoArr[38][0] = 3;
  genoArr[38][1] = 2;
  genoArr[38][2] = 3;

  // 39, T/TT
  Mij[39][12] = Mij[39][27] = Mij[39][33] = 1 / 9.0;
  Mij[39][21] = Mij[39][36] = Mij[39][38] = 2 / 9.0;
  Mij[39][53] = 2 * delRate / 3.0;
  Mij[39][59] = delRate / 3.0;

  genoArr[39][0] = 3;
  genoArr[39][1] = 3;
  genoArr[39][2] = 3;

  // 40, A/A
  Mij[40][0] = 2 * insRate / 3.0;
  Mij[40][42] = Mij[40][43] = Mij[40][44] = 2 / 9.0;
  Mij[40][60] = 2 * delRate / 3.0;

  genoArr[40][0] = 0;
  genoArr[40][1] = 0;
  genoArr[40][2] = INVALID;

  // 41, -/AA
  Mij[41][45] = Mij[41][46] = Mij[41][47] = 2 / 9.0;
  Mij[41][60] = 2 * delRate / 3.0;

  genoArr[41][0] = DELETED;
  genoArr[41][1] = 0;
  genoArr[41][2] = 0;

  // 42, A/C
  Mij[42][4] = Mij[42][7] = insRate / 3.0;
  Mij[42][40] = Mij[42][43] = Mij[42][44] = Mij[42][48] = Mij[42][49] =
      Mij[42][50] = 1 / 9.0;
  Mij[42][60] = Mij[42][61] = delRate / 3.0;

  genoArr[42][0] = 0;
  genoArr[42][1] = 1;
  genoArr[42][2] = INVALID;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 43, A/G
  Mij[43][5] = Mij[43][10] = insRate / 3.0;
  Mij[43][40] = Mij[43][42] = Mij[43][44] = Mij[43][49] = Mij[43][51] =
      Mij[43][52] = 1 / 9.0;
  Mij[43][60] = Mij[43][62] = delRate / 3.0;

  genoArr[43][0] = 0;
  genoArr[43][1] = 2;
  genoArr[43][2] = INVALID;

  // 44, A/T
  Mij[44][6] = Mij[44][12] = insRate / 3.0;
  Mij[44][40] = Mij[44][42] = Mij[44][43] = Mij[44][50] = Mij[44][52] =
      Mij[44][53] = 1 / 9.0;
  Mij[44][60] = Mij[44][63] = delRate / 3.0;

  genoArr[44][0] = 0;
  genoArr[44][1] = 3;
  genoArr[44][2] = INVALID;

  // 45, -/AC
  Mij[45][41] = Mij[45][46] = Mij[45][47] = Mij[45][54] = Mij[45][55] =
      Mij[45][56] = 1 / 9.0;
  Mij[45][60] = Mij[45][61] = delRate / 3.0;

  genoArr[45][0] = DELETED;
  genoArr[45][1] = 0;
  genoArr[45][2] = 1;

  // 46, -/AG
  Mij[46][41] = Mij[46][45] = Mij[46][47] = Mij[46][55] = Mij[46][57] =
      Mij[46][58] = 1 / 9.0;
  Mij[46][60] = Mij[46][62] = delRate / 3.0;

  genoArr[46][0] = DELETED;
  genoArr[46][1] = 0;
  genoArr[46][2] = 2;

  // 47, -/AT
  Mij[47][41] = Mij[47][45] = Mij[47][46] = Mij[47][56] = Mij[47][58] =
      Mij[47][59] = 1 / 9.0;
  Mij[47][60] = Mij[47][63] = delRate / 3.0;

  genoArr[47][0] = DELETED;
  genoArr[47][1] = 0;
  genoArr[47][2] = 3;

  // 48, C/C
  Mij[48][22] = 2 * insRate / 3.0;
  Mij[48][42] = Mij[48][49] = Mij[48][50] = 2 / 9.0;
  Mij[48][61] = 2 * delRate / 3.0;

  genoArr[48][0] = 1;
  genoArr[48][1] = 1;
  genoArr[48][2] = INVALID;

  // 49, C/G
  Mij[49][25] = Mij[49][28] = insRate / 3.0;
  Mij[49][42] = Mij[49][43] = Mij[49][48] = Mij[49][50] = Mij[49][51] =
      Mij[49][52] = 1 / 9.0;
  Mij[49][61] = Mij[49][62] = delRate / 3.0;

  genoArr[49][0] = 1;
  genoArr[49][1] = 2;
  genoArr[49][2] = INVALID;

  // 50, C/T
  Mij[50][27] = Mij[50][34] = insRate / 3.0;
  Mij[50][42] = Mij[50][44] = Mij[50][48] = Mij[50][49] = Mij[50][52] =
      Mij[50][53] = 1 / 9.0;
  Mij[50][61] = Mij[50][63] = delRate / 3.0;

  genoArr[50][0] = 1;
  genoArr[50][1] = 3;
  genoArr[50][2] = INVALID;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 51, G/G
  Mij[51][31] = 2 * insRate / 3.0;
  Mij[51][43] = Mij[51][49] = Mij[51][52] = 2 / 9.0;
  Mij[51][62] = delRate / 3.0;

  genoArr[51][0] = 2;
  genoArr[51][1] = 2;
  genoArr[51][2] = INVALID;

  // 52, G/T
  Mij[52][33] = Mij[52][37] = insRate / 3.0;
  Mij[52][43] = Mij[52][44] = Mij[52][49] = Mij[52][50] = Mij[52][51] =
      Mij[52][53] = 1 / 9.0;
  Mij[52][62] = Mij[52][63] = delRate / 3.0;

  genoArr[52][0] = 2;
  genoArr[52][1] = 3;
  genoArr[52][2] = INVALID;

  // 53, T/T
  Mij[53][39] = 2 * insRate / 3.0;
  Mij[53][44] = Mij[53][50] = Mij[53][52] = 2 / 9.0;
  Mij[53][63] = 2 * delRate / 3.0;

  genoArr[53][0] = 3;
  genoArr[53][1] = 3;
  genoArr[53][2] = INVALID;

  // 54, -/CC
  Mij[54][45] = Mij[54][55] = Mij[54][56] = 2 / 9.0;
  Mij[54][61] = 2 * delRate / 3.0;

  genoArr[54][0] = DELETED;
  genoArr[54][1] = 1;
  genoArr[54][2] = 1;

  // 55, -/CG
  Mij[55][45] = Mij[55][46] = Mij[55][54] = Mij[55][56] = Mij[55][57] =
      Mij[55][58] = 1 / 9.0;
  Mij[55][61] = Mij[55][62] = delRate / 3.0;

  genoArr[55][0] = DELETED;
  genoArr[55][1] = 1;
  genoArr[55][2] = 2;

  // 56, -/CT
  Mij[56][45] = Mij[56][47] = Mij[56][54] = Mij[56][55] = Mij[56][58] =
      Mij[56][59] = 1 / 9.0;
  Mij[56][61] = Mij[56][63] = delRate / 3.0;

  genoArr[56][0] = DELETED;
  genoArr[56][1] = 1;
  genoArr[56][2] = 3;

  // 57, -/GG
  Mij[57][46] = Mij[57][55] = Mij[57][58] = 2 / 9.0;
  Mij[57][62] = 2 * delRate / 3.0;

  genoArr[57][0] = DELETED;
  genoArr[57][1] = 2;
  genoArr[57][2] = 2;

  // 58, -/GT
  Mij[58][46] = Mij[58][47] = Mij[58][55] = Mij[58][56] = Mij[58][57] =
      Mij[58][59] = 1 / 9.0;
  Mij[58][62] = Mij[58][63] = delRate / 3.0;

  genoArr[58][0] = DELETED;
  genoArr[58][1] = 2;
  genoArr[58][2] = 3;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 59, -/TT
  Mij[59][47] = Mij[59][56] = Mij[59][58] = 2 / 9.0;
  Mij[59][63] = 2 * delRate / 3.0;

  genoArr[59][0] = DELETED;
  genoArr[59][1] = 3;
  genoArr[59][2] = 3;

  // 60, A/-
  Mij[60][41] = insRate / 3.0;
  Mij[60][61] = Mij[60][62] = Mij[60][63] = 1 / 9.0;
  Mij[60][64] = delRate / 3.0;

  genoArr[60][0] = 0;
  genoArr[60][1] = DELETED;
  genoArr[60][2] = INVALID;

  // 61, C/-
  Mij[61][54] = insRate / 3.0;
  Mij[61][60] = Mij[61][62] = Mij[61][63] = 1 / 9.0;
  Mij[61][64] = delRate / 3.0;

  genoArr[61][0] = 1;
  genoArr[61][1] = DELETED;
  genoArr[61][2] = INVALID;

  // 62, G/-
  Mij[62][57] = insRate / 3.0;
  Mij[62][60] = Mij[62][61] = Mij[62][63] = 1 / 9.0;
  Mij[62][64] = delRate / 3.0;

  genoArr[62][0] = 2;
  genoArr[62][1] = DELETED;
  genoArr[62][2] = INVALID;

  // 63, T/-
  Mij[63][59] = insRate / 3.0;
  Mij[63][60] = Mij[63][61] = Mij[63][62] = 1 / 9.0;
  Mij[63][64] = delRate / 3.0;

  genoArr[63][0] = 3;
  genoArr[63][1] = DELETED;
  genoArr[63][2] = INVALID;

  // 64, -
  genoArr[64][0] = genoArr[64][1] = DELETED;
  genoArr[64][2] = INVALID;

  hashMapIntArrToInt = CreateHashMapIntArrToInt(13);
  for (int i = 0; i < n; ++i)
    InsertKeyValIntArrToInt(hashMapIntArrToInt, genoArr[i], i);

  return hashMapIntArrToInt;
}

void FreeIntToIntArrHashMapMD() {
  for (int i = 0; i < 3; ++i) {
    struct NodeIntToIntArr *tmp = NULL;

    while (hashMapIntToIntArr->list[i]) {
      tmp = hashMapIntToIntArr->list[i];
      hashMapIntToIntArr->list[i] = hashMapIntToIntArr->list[i]->next;
      free(tmp);
      tmp = NULL;
    }
  }
  free(hashMapIntToIntArr->list);
  hashMapIntToIntArr->list = NULL;
  free(hashMapIntToIntArr);
  hashMapIntToIntArr = NULL;
}

void FreeIntArrToIntHashMapMD() {
  for (int i = 0; i < 3; ++i) {
    struct NodeIntArrToInt *tmp = NULL;

    while (hashMapIntArrToInt->list[i]) {
      tmp = hashMapIntArrToInt->list[i];
      hashMapIntArrToInt->list[i] = hashMapIntArrToInt->list[i]->next;
      free(tmp);
      tmp = NULL;
    }
  }
  free(hashMapIntArrToInt->list);
  hashMapIntArrToInt->list = NULL;
  free(hashMapIntArrToInt);
  hashMapIntArrToInt = NULL;
}

void FreeIntToIntArrHashMapMDI() {
  for (int i = 0; i < 13; ++i) {
    struct NodeIntToIntArr *tmp = NULL;

    while (hashMapIntToIntArr->list[i]) {
      tmp = hashMapIntToIntArr->list[i];
      hashMapIntToIntArr->list[i] = hashMapIntToIntArr->list[i]->next;
      free(tmp);
      tmp = NULL;
    }
  }
  free(hashMapIntToIntArr->list);
  hashMapIntToIntArr->list = NULL;
  free(hashMapIntToIntArr);
  hashMapIntToIntArr = NULL;
}

void FreeIntArrToIntHashMapMDI() {
  for (int i = 0; i < 13; ++i) {
    struct NodeIntArrToInt *tmp = NULL;

    while (hashMapIntArrToInt->list[i]) {
      tmp = hashMapIntArrToInt->list[i];
      hashMapIntArrToInt->list[i] = hashMapIntArrToInt->list[i]->next;
      free(tmp);
      tmp = NULL;
    }
  }
  free(hashMapIntArrToInt->list);
  hashMapIntArrToInt->list = NULL;
  free(hashMapIntArrToInt);
  hashMapIntArrToInt = NULL;
}

#endif // MYSIMULATOR_TESTFUNCTIONS_H
