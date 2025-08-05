//
// finiteMuInDelSimulator.h
// CellCoal
//
// Created by Senbai Kang on 11/20/19.
//

#ifndef MYSIMULATOR_FINITEMUINDELSIMULATOR_H
#define MYSIMULATOR_FINITEMUINDELSIMULATOR_H

#include "definitions.h"
#include "hashMap.h"
#include "rand.h"
#include "matrix_exponential.h"
#include <limits.h>
#include <omp.h>
#include <pthread.h>
#include <string.h>

typedef struct {
  long id;
  TreeNode *root;
  int numGenotypes;
} TreeRootRandomizer;

typedef struct {
  long id;
  CellStrExtd *cellStr;
} ThreadedCellStr;

typedef struct {
  FILE *fp;
  char fileName[MAX_NAME];
} RawDataFile;

typedef struct {
  long id;
  unsigned int numCells;
  RandomDevice *rd;

  unsigned int *cellIndices;
  RawDataFile *files;
} ThreadedRawDataPrinter;

int GetMaxNumExistingAllelesPerStrand();
int GetMaxInt(const int *arr, size_t size);
int GetSumInt(const int *arr, size_t size);
int GetShapeController(const int *alleles, int numExistingAlleles);
void PrepareGlobalFilesExtd(int argc, char **argv);
void PrepareSeparateFoldersAndFiles(int replicate);
void PrepareSeparateFilesExtd(int replicate);
void PrepareRawReadSeqDataFiles(int replicate, int numTumCells);
void InitializeGenomesExtd(TreeNode *p, RandomDevice *rd);
//extern int EigenREVMD(double RootExtd[], double CijkExtd[]);
//extern int EigenREVMDI(double RootExtd[], double CijkExtd[]);
void EvolveSitesOnTreeExtd(TreeNode *treeRoot, RandomDevice *rd);
void ThreadedSimulateFiniteMuInDel(TreeNode *treeRoot, int nrOfGenotypes);
void *ThreadedSimulateFiniteMuInDelForSite(void *arg);
void SimulateFiniteMuInDelForSite(TreeNode *treeRoot, ThreadStr *threadStr,
                                  int site, int nrOfGenotypes);
int ChooseUniformStateExtd(const double *prob, int start_pos, RandomDevice *rd);
void GenerateReadCountsExtd(int replicate, RandomDevice *rd);
void SiteReadCountsExtd(CellSiteStrExtd *c, double sizeFactor, int j,
                        double *probs, double **ngsEij, double **ampEij,
                        RandomDevice *rd);
void MakeDoubletsExtd(double *probs, double **ngsEij, double **ampEij, RandomDevice *rd);
int GenerateSequencingCoverage(int numExistingAlleles, double alleleCoverage,
                               double sizeFactor, double alleleRawVariance,
                               RandomDevice *rd);
// int RandomBetaBinomial(double mean, double overdispersion, int numTrials,
//                        RandomDevice *rd);
// void RandomDirichletMultinomial(double *means, double overdispersion,
//                                 int numTrials, int pivot, size_t size,
//                                 RandomDevice *rd, int *readCounts);
// void GenotypeLikelihoodsExtd(CellSiteStrExtd *c, int i, int j, double
// **ngsEij,
//                             double **ampEij);
int CountTrueEvents();
void ProcessCoverageProb(int replicate);
void PrintTumCellNames(FILE *fp, char col_sep, int includeNormalCell);
void PrintSNVSitesNames(FILE *fp);
void PrintSNVAllelicInformation(FILE *fp);
void PrintAdoStates(FILE *fp);
void PrintParallelMut(FILE *fp);
void PrintSizeFactors(FILE *fp);
void PrintTrueSNV(FILE *fpAl, FILE *fpNr, FILE *fpMt);
void PrintEvolutionaryEventsSummary(FILE *fp);
void PrintRawData(RawDataFile *files);
void PrintRawDataSam();
void *ThreadedPrintRawDataSam(void *arg);
char GetANucleotide(int *readCounts, RandomDevice *rd);
void PrintTrueSNVCov(FILE *fp);
void PrintMPileUp(FILE *fpNM, FILE *fpNU, FILE *fpM, FILE *fpU);
void PrintSNVGenotypesExtd(FILE *fp, int toStderrFlag);
// int GenerateAnAlternativeNucleotide(int referenceAllele, RandomDevice *rd);
void PrintFullGenotypesExtd(FILE *fp, int toStderrFlag);
void PrintSiteInfoExtd(FILE *fp, int i);
void PrintRunInformationExtd(FILE *fp);
void PrintSummaryToLog(FILE *fp, int replicate, int toStderrFlag);

/**
 * 0: read length of both ends
 * 1: read length of internal fragments
 */
static int readLengths[2] = {0, 150};
int numReadBlocks;
int *startSiteBlocks;

pthread_t *threads;
int numThreads;
int *sitesPointsForThreads;
ThreadStr *threadStr;
RandomDevice *randDevs; // If using multithreading, one more random device
                        // should be left for the main thread.
RandomDevice *mainThreadRandDev;
int fixedTree;

ThreadedRawDataPrinter *rawDataPrinters;

// always print out
char cellNamesFile[MAX_NAME], cellNamesDir[MAX_NAME], SNVSitesNamesFile[MAX_NAME],
    SNVSitesNamesDir[MAX_NAME], SNVAllelicSeqInfoFile[MAX_NAME],
    SNVAllelicSeqInfoDir[MAX_NAME];
FILE *fpSNVSitesNames, *fpSNVAllelicSeqInfo;

char trueSNVGenotypesAlFile[MAX_NAME], trueSNVGenotypesNrFile[MAX_NAME],
    trueSNVTernaryMatrixFile[MAX_NAME], trueSNVDir[MAX_NAME],
    evolutionaryEventsFile[MAX_NAME], rawSeqDataRefFile[MAX_NAME],
    rawSeqDataDir[MAX_NAME], dummyBamDataDir[MAX_NAME],
    mpileupDataNMFile[MAX_NAME], mpileupDataMFile[MAX_NAME],
    mpileupDataNUFile[MAX_NAME], mpileupDataUFile[MAX_NAME],
    mpileupDataDir[MAX_NAME], parallelMutDir[MAX_NAME],
    parallelMutFile[MAX_NAME], adoStatesDir[MAX_NAME], adoStatesFile[MAX_NAME],
    trueSNVCovFile[MAX_NAME];
FILE *fpTrueSNVGenotypesAl, *fpTrueSNVGenotypesNr, *fpTrueSNVTernaryMatrix,
    *fpMpileupDataNM, *fpMpileupDataM, *fpMpileupDataNU, *fpMpileupDataU,
    *fpParallelMut, *fpAdoStates,
    *fpTrueSNVCov, *fpEvolutionaryEvents;

/* size: (numCells + 2) without doublet or (2 * numCells + 2) with doublet
 * 0: reference sequence;
 * 1 to numCells: cancerous cells
 * (numCells + 1): healthy cell
 * (numCells + 2) to (2 * numCells + 1): duplicate cancerous cells
 */
RawDataFile *rawDataFiles;

/**
 * Tree structure represented as a matrix with non-negative elements.
 * Each element represents the relative parent-child relationship between its
 * row and column. 0: no parent-child relationship. 1: the row node itself. 2:
 * the column node is a direct child of the row node. 3: the column node is a
 * grandchild of the row node. 4: ... and so on.
 */
unsigned **treeStructure;

char **cellNamesExtd;
int ADOtype;
double MijMD[15][15] = {0.0};  // mutation, deletion
double MijMDI[35][35] = {0.0}; // mutation, deletion, insertion
int constrained = NO; // 0: not constrained (allowing more than 1 deletion); 1: constrained (allowing at most 1 deletion)
int seqCovCellCoalStyle = NO; // whether to simulate sequencing coverage following CellCoal's style?
double allelicCovMean, allelicCovVar, allelicRawVarMean, allelicRawVarVar;
double *allelicCov, *allelicRawVar; // for all sites
double *sizeFactors;
double meanSizeFactor, varSizeFactor;
double shapeController[2] = {0.0, 0.0};
static int cumNucleotides[4] = {1, 2, 3, 4};
static int fixedAmplificationError = 1;
static double numericalStabilizer = 1e-15;
static int doPrintHealthyTip = NO;
static int doCoalescenceTree = YES;
static int doPrintFullNucleotideReads =
    NO; // save reads in a form of: 0, (A, C, G, T), or: 1, (c, m)

FILE *fpTreeWithoutTrunk, *fpTreeWithHealthyTip, *fpCellNames, *fpSizeFactors;
char treeWithTrunkDir[MAX_NAME], treeWithoutTrunkDir[MAX_NAME],
    treeWithoutTrunkFile[MAX_NAME], treeWithoutTrunkFileCopy[MAX_NAME],
    treeWithHealthyTipDir[MAX_NAME], treeWithHealthyTipFile[MAX_NAME],
    sizeFactorsDir[MAX_NAME], sizeFactorsFile[MAX_NAME], unmergedDoubletsDir[MAX_NAME];

static const int candidateMLGenotypesMD[NUM_CAN_ML_GENOTYPES] = {0, 1, 2, 3, 4,
                                                                 5, 6, 7, 8, 9};
static const char genotypesInCellPhy[NUM_CAN_ML_GENOTYPES] = {
    'A', 'M', 'R', 'W', 'C', 'S', 'Y', 'G', 'K', 'T'};
static double insertionRate = 0.0;
struct HashMapIntToIntArr *genotypeHashMapIntToIntArr;
struct HashMapIntArrToInt *genotypeHashMapIntArrToInt;
int genotypeArrMD[15][3];
int genotypeArrMDI[35][3];
int numINS;
double cumNumINS, meanNumINS, varNumINS, cumNumINSSq;
int *mutatedSites, *insertedSites, *deletedSites;
CellStrExtd *cellExtd;
long completeDeletions, completeADOs, noCoverage, noReadsPositiveCoverage,
    zeroCoverage, zeroProportionVAF, completeDeletionsAll, completeADOsAll,
    noCoverageAll, noReadsPositiveCoverageAll, zeroCoverageAll,
    zeroProportionVAFAll;
/* for SNV sites
 * 1. the index of each row is the number of existing alleles
 * 2. the fist column is the total number of zero coverage occurrences of a
 * specific number of existing alleles for the current data set (for different
 * unfixed trees) or for all data sets (for a fixed tree)
 * 3. the second column has the same meaning as the first column does, only
 * those sites counted with no reads (representing non-coverage) */
int allelicZeroCoverageProportion[4][2], allelicZeroCoverageProportionAll[4][2];
/* the probability of a successful trial for SNV sites
 * 1. the first dimension is the number of existing alleles (0, 1, 2, 3).
 * 2. the second dimension contains the minimum (0) and the maximum (1) number
 * of probability or number of successful trials, and the corresponding mean (2)
 * as well as the variance (3). */
double probNBInfo[4][4], probNBInfoAll[4][4];
/* details of the probability of a successful trial for SNV sites
 * 1. the first dimension is the number of existing alleles (0, 1, 2, 3).
 * 2. the second dimension is the concrete value of probability for each site in
 * each cell. The space will be dynamically allocated at run time. */
double *probNBDetails[4], *probNBDetailsAll[4];

/*********************
 *      Methods      *
 *********************/

void AllocateTreeStructure(int numNodes) {
  treeStructure = (unsigned **)calloc(numNodes + 1, sizeof(unsigned *));
  for (int i = 0; i < numNodes; i++) {
    treeStructure[i] = (unsigned *)calloc(numNodes + 1, sizeof(unsigned));
  }
} // AllocateTreeStructure

void InitRandomDevices(const unsigned long seed) {
  const int size = numThreads > 1 ? (numThreads + 1) : 1;
  randDevs = AllocateMem(size + 2);
  if (!randDevs) {
    fprintf(stderr, "Could not allocate the randDevs structure\n");
    exit(-1);
  }

  //  fprintf(stderr, "%p\n", GetRandDevElement(randDevs, 0));
  //  fprintf(stderr, "%p\n", GetRandDevElement(randDevs, 1));
  //  fprintf(stderr, "%p\n", GetRandDevElement(randDevs, 2));
  //  fprintf(stderr, "%p\n", GetRandDevElement(randDevs, 3));
  //  fprintf(stderr, "%p\n", GetRandDevElement(randDevs, 4));

  unsigned long seeds[size];
  seeds[0] = seed;
  if (size > 1) {
    RandomDevice *tmp = GetRandDevElement(randDevs, size);
    Call_C_SetRandomDevice(seed, tmp);
    RandIntegers(size, seeds, tmp);
  }

  for (int i = 0; i < size; i++) {
    Call_C_SetRandomDevice(seeds[i], GetRandDevElement(randDevs, i));
  }

  mainThreadRandDev = GetRandDevElement(randDevs, size - 1);
} // InitRandomDevices

void InitTreeStructure(int numNodes) {
  for (int i = 0; i < numNodes; i++) {
    for (int j = 0; j < numNodes; j++) {
      if (i == j) {
        treeStructure[i][i] = 1;
      } else {
        treeStructure[i][j] = 0;
      }
    }
  }
} // InitTreeStructure

void SetTreeStructure(const TreeNode *node, int numNodes) {
  if (node == NULL)
    return;

  SetTreeStructure(node->left, numNodes);
  SetTreeStructure(node->right, numNodes);

  if (node->left != NULL) {
    for (int i = 0; i < numNodes; i++) {
      if (treeStructure[node->left->label][i] > 0) {
        treeStructure[node->label][i] +=
            treeStructure[node->left->label][i] + 1;
      }
    }
  }

  if (node->right != NULL) {
    for (int i = 0; i < numNodes; i++) {
      if (treeStructure[node->right->label][i] > 0) {
        treeStructure[node->label][i] +=
            treeStructure[node->right->label][i] + 1;
      }
    }
  }
} // SetTreeStructure

void PrintTreeStructure(int numNodes) {
  fprintf(stdout, "\n |");
  for (int i = 0; i < numNodes; i++) {
    fprintf(stdout, "%d ", i);
  }
  fprintf(stdout, "\n");

  char sep[2 * numNodes + 2];
  memset(sep, '_', sizeof(sep) + 1);
  sep[2 * numNodes + 1] = '\0';
  fprintf(stdout, "%s\n", sep);

  for (int i = 0; i < numNodes; i++) {
    fprintf(stdout, "%d|", i);
    for (int j = 0; j < numNodes; j++) {
      fprintf(stdout, "%d ", treeStructure[i][j]);
    }
    fprintf(stdout, "\n");
  }
} // PrintTreeStructure

void FreeTreeStructure(int numNodes) {
  for (int i = 0; i < numNodes; i++) {
    free(treeStructure[i]);
    treeStructure[i] = NULL;
  }
  free(treeStructure);
  treeStructure = NULL;
} // AllocateTreeStructure

void SetReadLengths(int siteNr) {
  if (siteNr % readLengths[1] == 0) {
    readLengths[0] = readLengths[1];
    numReadBlocks = siteNr / readLengths[0];
  } else {
    int readsNr = (int)(((double)siteNr) / readLengths[1] - 1);
    readLengths[0] = (siteNr - readsNr * readLengths[1]) / 2;
    numReadBlocks = readsNr + 2;
  }

  startSiteBlocks = (int *)malloc((numReadBlocks + 1) * sizeof(int));
  startSiteBlocks[0] = 0;
  for (int i = 0; i < numReadBlocks; i++) {
    startSiteBlocks[i + 1] =
        startSiteBlocks[i] +
        (i == 0 || i == numReadBlocks - 1 ? readLengths[0] : readLengths[1]);
  }
}

/* Get the total number of genotype states. */
int GetGenotypeTotalNum(int model, int constrained) {
  if (model == FiniteMu) {
    return 15;
  } else if (model == FiniteMuDel) {
    if (constrained == YES) {
      return 14;
    } else {
      return 15;
    }
  } else if (model == FiniteMuInDel) {
    if (constrained == YES) {
      return 34;
    } else {
      return 35;
    }
  } else {
    fprintf(stderr, "Unknown model type: %d\n", model);
    exit(-1);
  }
}

/* Get the rate matrix according to the evolutionary model. */
const double* GetRateMtx(int model) {
  if (model == FiniteMu || model == FiniteMuDel) {
    return &(MijMD[0][0]);
  } else if (model == FiniteMuInDel) {
    return &(MijMDI[0][0]);
  } else {
    fprintf(stderr, "Unknown model type: %d\n", model);
    exit(-1);
  }
}

void ProbMtxSanityCheck(const double *Pij, int nrOfGenotypes) {
  for (int i = 0; i < nrOfGenotypes; i++) {
    double rowSum = 0.0;

    for (int j = 0; j < nrOfGenotypes; j++) {
      rowSum += Pij[i * nrOfGenotypes + j];
    }

    if (fabs(rowSum - 1.0) > 1e-12) {
      fprintf(stderr,
              "\n!!! PROBABILITY MATRIX ERROR: "
              "The sum of row %d is %.12f, not 1. \n\n",
              i, rowSum);
      exit(1);
    }
  }
}

double *GetTransitionMtx(const int model, int stateNum, const void *rateMtx, double branchLength) {
  double *EPij, mtxProd[stateNum * stateNum];
  const double (*mtx)[GetGenotypeTotalNum(model, NO)] = rateMtx;

  for (int i = 0; i < stateNum; i++)
    for (int j = 0; j < stateNum; j++)
      mtxProd[i * stateNum + j] = branchLength * mtx[i][j];

  EPij = r8mat_expm1(stateNum, mtxProd);
  ProbMtxSanityCheck(EPij, stateNum);
  return EPij;
}

/*
 * Initialize the hash map store the information of mapping codes to genotypes.
 * For example, 0 -> A/A
 * */
void InitGenotypeHashMapMD() {
  genotypeHashMapIntToIntArr = CreateHashMapIntToIntArr(3);
  genotypeHashMapIntArrToInt = CreateHashMapIntArrToInt(3);

  for (int i = 0; i < 15; ++i) {
    InsertKeyValIntToIntArr(genotypeHashMapIntToIntArr, i, genotypeArrMD[i]);
    InsertKeyValIntArrToInt(genotypeHashMapIntArrToInt, genotypeArrMD[i], i);
  }
}

/*
 * mutation, deletion
 * define the substitution rate matrix internally.
 * users can configure the mutation rate, deletion rate and as
 * they wish. mutationRate is defined as the rate of a reference nucleotide
 * mutating to any one of the three alternative nucleotides.
 * If deletions is constrained, the genotype state - is not allowed.
 */
void InitRateMatrixMD(double mutationRate, double deletionRate) {
  // normalize the deletionRate and insertionRate regarding the mutationRatE
  deletionRate /= (mutationRate + 1e-15);

  // initialize all entries to 0.0
  for (int i = 0; i < 15; ++i)
    for (int j = 0; j < 15; ++j)
      MijMD[i][j] = 0.0;

  // 0    1    2    3    4    5    6    7    8    9    10   11   12   13   14
  // A/A, A/C, A/G, A/T, C/C, C/G, C/T, G/G, G/T, T/T, A/-, C/-, G/-, T/-, -

  // 0, A/A
  MijMD[0][1] = MijMD[0][2] = MijMD[0][3] = 1.0 / 3;
  MijMD[0][10] = deletionRate;

  genotypeArrMD[0][0] = 0;
  genotypeArrMD[0][1] = 0;
  genotypeArrMD[0][2] = INVALID;

  // 1, A/C
  MijMD[1][0] = MijMD[1][2] = MijMD[1][3] = MijMD[1][4] = MijMD[1][5] =
      MijMD[1][6] = 1.0 / 6;
  MijMD[1][10] = MijMD[1][11] = deletionRate / 2;

  genotypeArrMD[1][0] = 0;
  genotypeArrMD[1][1] = 1;
  genotypeArrMD[1][2] = INVALID;

  // 2, A/G
  MijMD[2][0] = MijMD[2][1] = MijMD[2][3] = MijMD[2][5] = MijMD[2][7] =
      MijMD[2][8] = 1.0 / 6;
  MijMD[2][10] = MijMD[2][12] = deletionRate / 2;

  genotypeArrMD[2][0] = 0;
  genotypeArrMD[2][1] = 2;
  genotypeArrMD[2][2] = INVALID;

  // 3, A/T
  MijMD[3][0] = MijMD[3][1] = MijMD[3][2] = MijMD[3][6] = MijMD[3][8] =
      MijMD[3][9] = 1.0 / 6;
  MijMD[3][10] = MijMD[3][13] = deletionRate / 2;

  genotypeArrMD[3][0] = 0;
  genotypeArrMD[3][1] = 3;
  genotypeArrMD[3][2] = INVALID;

  // 4, C/C
  MijMD[4][1] = MijMD[4][5] = MijMD[4][6] = 1.0 / 3;
  MijMD[4][11] = deletionRate;

  genotypeArrMD[4][0] = 1;
  genotypeArrMD[4][1] = 1;
  genotypeArrMD[4][2] = INVALID;

  // 5, C/G
  MijMD[5][1] = MijMD[5][2] = MijMD[5][4] = MijMD[5][6] = MijMD[5][7] =
      MijMD[5][8] = 1.0 / 6;
  MijMD[5][11] = MijMD[5][12] = deletionRate / 2;

  genotypeArrMD[5][0] = 1;
  genotypeArrMD[5][1] = 2;
  genotypeArrMD[5][2] = INVALID;

  // 6, C/T
  MijMD[6][1] = MijMD[6][3] = MijMD[6][4] = MijMD[6][5] = MijMD[6][8] =
      MijMD[6][9] = 1.0 / 6;
  MijMD[6][11] = MijMD[6][13] = deletionRate / 2;

  genotypeArrMD[6][0] = 1;
  genotypeArrMD[6][1] = 3;
  genotypeArrMD[6][2] = INVALID;

  // 7, G/G
  MijMD[7][2] = MijMD[7][5] = MijMD[7][8] = 1.0 / 3;
  MijMD[7][12] = deletionRate;

  genotypeArrMD[7][0] = 2;
  genotypeArrMD[7][1] = 2;
  genotypeArrMD[7][2] = INVALID;

  // 8, G/T
  MijMD[8][2] = MijMD[8][3] = MijMD[8][5] = MijMD[8][6] = MijMD[8][7] =
      MijMD[8][9] = 1.0 / 6;
  MijMD[8][12] = MijMD[8][13] = deletionRate / 2;

  genotypeArrMD[8][0] = 2;
  genotypeArrMD[8][1] = 3;
  genotypeArrMD[8][2] = INVALID;

  // 9, T/T
  MijMD[9][3] = MijMD[9][6] = MijMD[9][8] = 1.0 / 3;
  MijMD[9][13] = deletionRate;

  genotypeArrMD[9][0] = 3;
  genotypeArrMD[9][1] = 3;
  genotypeArrMD[9][2] = INVALID;

  // 10, A/-
  MijMD[10][11] = MijMD[10][12] = MijMD[10][13] = 1.0 / 6;
  if (constrained == NO) {
    MijMD[10][14] = deletionRate / 2;
  }

  genotypeArrMD[10][0] = 0;
  genotypeArrMD[10][1] = DELETED;
  genotypeArrMD[10][2] = INVALID;

  // 11, C/-
  MijMD[11][10] = MijMD[11][12] = MijMD[11][13] = 1.0 / 6;
  if (constrained == NO) {
    MijMD[11][14] = deletionRate / 2;
  }

  genotypeArrMD[11][0] = 1;
  genotypeArrMD[11][1] = DELETED;
  genotypeArrMD[11][2] = INVALID;

  // 12, G/-
  MijMD[12][10] = MijMD[12][11] = MijMD[12][13] = 1.0 / 6;
  if (constrained == NO) {
    MijMD[12][14] = deletionRate / 2;
  }

  genotypeArrMD[12][0] = 2;
  genotypeArrMD[12][1] = DELETED;
  genotypeArrMD[12][2] = INVALID;

  // 13, T/-
  MijMD[13][10] = MijMD[13][11] = MijMD[13][12] = 1.0 / 6;
  if (constrained == NO) {
    MijMD[13][14] = deletionRate / 2;
  }

  genotypeArrMD[13][0] = 3;
  genotypeArrMD[13][1] = DELETED;
  genotypeArrMD[13][2] = INVALID;

  // 14, -
  if (constrained == NO) {
    genotypeArrMD[14][0] = genotypeArrMD[14][1] = DELETED;
  } else {
    genotypeArrMD[14][0] = genotypeArrMD[14][1] = INVALID;
  }
  genotypeArrMD[14][2] = INVALID;

  for (int i = 0; i < 15; i++) {
    for (int j = 0; j < 15; j++) {
      if (i != j) {
        MijMD[i][i] -= MijMD[i][j];
      }
    }
  }

  InitGenotypeHashMapMD();
}

/*
 * Initialize the hash map store the information of mapping codes to genotypes.
 * For example, 0 -> AAA
 * */
void InitGenotypeHashMapMDI() {
  genotypeHashMapIntToIntArr = CreateHashMapIntToIntArr(7);
  genotypeHashMapIntArrToInt = CreateHashMapIntArrToInt(7);

  for (int i = 0; i < 35; ++i) {
    InsertKeyValIntToIntArr(genotypeHashMapIntToIntArr, i, genotypeArrMDI[i]);
    InsertKeyValIntArrToInt(genotypeHashMapIntArrToInt, genotypeArrMDI[i], i);
  }
}

/*
 * mutation, deletion, insertion
 * define the substitution rate matrix internally.
 * users can configure the mutation rate, deletion rate and insertion rate as
 * they wish. mutationRate is defined as the rate of a reference nucleotide
 * mutating to any one of the three alternative nucleotides.
 * If deletions is constrained, the genotype state - is not allowed.
 */
void InitRateMatrixMDI(double mutationRate, double deletionRate,
                       double insertionRate) {
  // normalize the deletionRate and insertionRate regarding the mutationRate
  deletionRate /= (mutationRate + 1e-15);
  insertionRate /= (mutationRate + 1e-15);

  // initialize all entries to 0.0
  for (int i = 0; i < 35; ++i)
    for (int j = 0; j < 35; ++j)
      MijMDI[i][j] = 0.0;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // AAA   CCC   GGG   TTT   AAC   AAG   AAT   ACC   CCG   CCT   AGG   CGG

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // GGT   ATT   CTT   GTT   ACG   ACT   AGT   CGT   AA    AC    AG    AT

  // 24    25    26    27    28    29    30    31    32    33    34
  // CC    CG    CT    GG    GT    TT    A-    C-    G-    T-    -

  // 0, AAA
  MijMDI[0][4] = MijMDI[0][5] = MijMDI[0][6] = 1.0 / 2;
  MijMDI[0][20] = 3.0 * deletionRate / 2;

  genotypeArrMDI[0][0] = 0;
  genotypeArrMDI[0][1] = 0;
  genotypeArrMDI[0][2] = 0;

  // 1, CCC
  MijMDI[1][7] = MijMDI[1][8] = MijMDI[1][9] = 1.0 / 2;
  MijMDI[1][24] = 3.0 * deletionRate / 2;

  genotypeArrMDI[1][0] = 1;
  genotypeArrMDI[1][1] = 1;
  genotypeArrMDI[1][2] = 1;

  // 2, GGG
  MijMDI[2][10] = MijMDI[2][11] = MijMDI[2][12] = 1.0 / 2;
  MijMDI[2][27] = 3.0 * deletionRate / 2;

  genotypeArrMDI[2][0] = 2;
  genotypeArrMDI[2][1] = 2;
  genotypeArrMDI[2][2] = 2;

  // 3, TTT
  MijMDI[3][13] = MijMDI[3][14] = MijMDI[3][15] = 1.0 / 2;
  MijMDI[3][29] = 3.0 * deletionRate / 2;

  genotypeArrMDI[3][0] = 3;
  genotypeArrMDI[3][1] = 3;
  genotypeArrMDI[3][2] = 3;

  // 4, AAC
  MijMDI[4][0] = MijMDI[4][5] = MijMDI[4][6] = 1.0 / 6;
  MijMDI[4][7] = MijMDI[4][16] = MijMDI[4][17] = 1.0 / 3;
  MijMDI[4][20] = deletionRate / 2;
  MijMDI[4][21] = deletionRate;

  genotypeArrMDI[4][0] = 0;
  genotypeArrMDI[4][1] = 0;
  genotypeArrMDI[4][2] = 1;

  // 5, AAG
  MijMDI[5][0] = MijMDI[5][4] = MijMDI[5][6] = 1.0 / 6;
  MijMDI[5][10] = MijMDI[5][16] = MijMDI[5][18] = 1.0 / 3;
  MijMDI[5][20] = deletionRate / 2;
  MijMDI[5][22] = deletionRate;

  genotypeArrMDI[5][0] = 0;
  genotypeArrMDI[5][1] = 0;
  genotypeArrMDI[5][2] = 2;

  // 6, AAT
  MijMDI[6][0] = MijMDI[6][4] = MijMDI[6][5] = 1.0 / 6;
  MijMDI[6][13] = MijMDI[6][17] = MijMDI[6][18] = 1.0 / 3;
  MijMDI[6][20] = deletionRate / 2;
  MijMDI[6][23] = deletionRate;

  genotypeArrMDI[6][0] = 0;
  genotypeArrMDI[6][1] = 0;
  genotypeArrMDI[6][2] = 3;

  // 7, ACC
  MijMDI[7][1] = MijMDI[7][8] = MijMDI[7][9] = 1.0 / 6;
  MijMDI[7][4] = MijMDI[7][16] = MijMDI[7][17] = 1.0 / 3;
  MijMDI[7][21] = deletionRate;
  MijMDI[7][24] = deletionRate / 2;

  genotypeArrMDI[7][0] = 0;
  genotypeArrMDI[7][1] = 1;
  genotypeArrMDI[7][2] = 1;

  // 8, CCG
  MijMDI[8][1] = MijMDI[8][7] = MijMDI[8][9] = 1.0 / 6;
  MijMDI[8][11] = MijMDI[8][16] = MijMDI[8][19] = 1.0 / 3;
  MijMDI[8][24] = deletionRate / 2;
  MijMDI[8][25] = deletionRate;

  genotypeArrMDI[8][0] = 1;
  genotypeArrMDI[8][1] = 1;
  genotypeArrMDI[8][2] = 2;

  // 9, CCT
  MijMDI[9][1] = MijMDI[9][7] = MijMDI[9][8] = 1.0 / 6;
  MijMDI[9][14] = MijMDI[9][17] = MijMDI[9][19] = 1.0 / 3;
  MijMDI[9][24] = deletionRate / 2;
  MijMDI[9][26] = deletionRate;

  genotypeArrMDI[9][0] = 1;
  genotypeArrMDI[9][1] = 1;
  genotypeArrMDI[9][2] = 3;

  // 10, AGG
  MijMDI[10][2] = MijMDI[10][11] = MijMDI[10][12] = 1.0 / 6;
  MijMDI[10][5] = MijMDI[10][16] = MijMDI[10][18] = 1.0 / 3;
  MijMDI[10][22] = deletionRate;
  MijMDI[10][27] = deletionRate / 2;

  genotypeArrMDI[10][0] = 0;
  genotypeArrMDI[10][1] = 2;
  genotypeArrMDI[10][2] = 2;

  // 11, CGG
  MijMDI[11][2] = MijMDI[11][10] = MijMDI[11][12] = 1.0 / 6;
  MijMDI[11][8] = MijMDI[11][16] = MijMDI[11][19] = 1.0 / 3;
  MijMDI[11][25] = deletionRate;
  MijMDI[11][27] = deletionRate / 2;

  genotypeArrMDI[11][0] = 1;
  genotypeArrMDI[11][1] = 2;
  genotypeArrMDI[11][2] = 2;

  // 12, GGT
  MijMDI[12][2] = MijMDI[12][10] = MijMDI[12][11] = 1.0 / 6;
  MijMDI[12][15] = MijMDI[12][18] = MijMDI[12][19] = 1.0 / 3;
  MijMDI[12][27] = deletionRate / 2;
  MijMDI[12][28] = deletionRate;

  genotypeArrMDI[12][0] = 2;
  genotypeArrMDI[12][1] = 2;
  genotypeArrMDI[12][2] = 3;

  // 13, ATT
  MijMDI[13][3] = MijMDI[13][14] = MijMDI[13][15] = 1.0 / 6;
  MijMDI[13][6] = MijMDI[13][17] = MijMDI[13][18] = 1.0 / 3;
  MijMDI[13][23] = deletionRate;
  MijMDI[13][29] = deletionRate / 2;

  genotypeArrMDI[13][0] = 0;
  genotypeArrMDI[13][1] = 3;
  genotypeArrMDI[13][2] = 3;

  // 14, CTT
  MijMDI[14][3] = MijMDI[14][13] = MijMDI[14][15] = 1.0 / 6;
  MijMDI[14][9] = MijMDI[14][17] = MijMDI[14][19] = 1.0 / 3;
  MijMDI[14][26] = deletionRate;
  MijMDI[14][29] = deletionRate / 2;

  genotypeArrMDI[14][0] = 1;
  genotypeArrMDI[14][1] = 3;
  genotypeArrMDI[14][2] = 3;

  // 15, GTT
  MijMDI[15][3] = MijMDI[15][13] = MijMDI[15][14] = 1.0 / 6;
  MijMDI[15][12] = MijMDI[15][18] = MijMDI[15][19] = 1.0 / 3;
  MijMDI[15][28] = deletionRate;
  MijMDI[15][29] = deletionRate / 2;

  genotypeArrMDI[15][0] = 2;
  genotypeArrMDI[15][1] = 3;
  genotypeArrMDI[15][2] = 3;

  // 16, ACG
  MijMDI[16][4] = MijMDI[16][5] = MijMDI[16][7] = MijMDI[16][8] =
      MijMDI[16][10] = MijMDI[16][11] = MijMDI[16][17] = MijMDI[16][18] =
          MijMDI[16][19] = 1.0 / 6;
  MijMDI[16][21] = MijMDI[16][22] = MijMDI[16][25] = deletionRate / 2;

  genotypeArrMDI[16][0] = 0;
  genotypeArrMDI[16][1] = 1;
  genotypeArrMDI[16][2] = 2;

  // 17, ACT
  MijMDI[17][4] = MijMDI[17][6] = MijMDI[17][7] = MijMDI[17][9] =
      MijMDI[17][13] = MijMDI[17][14] = MijMDI[17][16] = MijMDI[17][18] =
          MijMDI[17][19] = 1.0 / 6;
  MijMDI[17][21] = MijMDI[17][23] = MijMDI[17][26] = deletionRate / 2;

  genotypeArrMDI[17][0] = 0;
  genotypeArrMDI[17][1] = 1;
  genotypeArrMDI[17][2] = 3;

  // 18, AGT
  MijMDI[18][5] = MijMDI[18][6] = MijMDI[18][10] = MijMDI[18][12] =
      MijMDI[18][13] = MijMDI[18][15] = MijMDI[18][16] = MijMDI[18][17] =
          MijMDI[18][19] = 1.0 / 6;
  MijMDI[18][22] = MijMDI[18][23] = MijMDI[18][28] = deletionRate / 2;

  genotypeArrMDI[18][0] = 0;
  genotypeArrMDI[18][1] = 2;
  genotypeArrMDI[18][2] = 3;

  // 19, CGT
  MijMDI[19][8] = MijMDI[19][9] = MijMDI[19][11] = MijMDI[19][12] =
      MijMDI[19][14] = MijMDI[19][15] = MijMDI[19][16] = MijMDI[19][17] =
          MijMDI[19][18] = 1.0 / 6;
  MijMDI[19][25] = MijMDI[19][26] = MijMDI[19][28] = deletionRate / 2;

  genotypeArrMDI[19][0] = 1;
  genotypeArrMDI[19][1] = 2;
  genotypeArrMDI[19][2] = 3;

  // 20, AA
  MijMDI[20][0] = insertionRate;
  MijMDI[20][21] = MijMDI[20][22] = MijMDI[0][23] = 1.0 / 3;
  MijMDI[20][30] = deletionRate;

  genotypeArrMDI[20][0] = 0;
  genotypeArrMDI[20][1] = 0;
  genotypeArrMDI[20][2] = INVALID;

  // 21, AC
  MijMDI[21][4] = MijMDI[21][7] = insertionRate / 2;
  MijMDI[21][20] = MijMDI[21][22] = MijMDI[21][23] = MijMDI[21][24] =
      MijMDI[21][25] = MijMDI[21][26] = 1.0 / 6;
  MijMDI[21][30] = MijMDI[21][31] = deletionRate / 2;

  genotypeArrMDI[21][0] = 0;
  genotypeArrMDI[21][1] = 1;
  genotypeArrMDI[21][2] = INVALID;

  // 22, AG
  MijMDI[22][5] = MijMDI[22][10] = insertionRate / 2;
  MijMDI[22][20] = MijMDI[22][21] = MijMDI[22][23] = MijMDI[22][25] =
      MijMDI[22][27] = MijMDI[22][28] = 1.0 / 6;
  MijMDI[22][30] = MijMDI[22][32] = deletionRate / 2;

  genotypeArrMDI[22][0] = 0;
  genotypeArrMDI[22][1] = 2;
  genotypeArrMDI[22][2] = INVALID;

  // 23, AT
  MijMDI[23][6] = MijMDI[23][13] = insertionRate / 2;
  MijMDI[23][20] = MijMDI[23][21] = MijMDI[23][22] = MijMDI[23][26] =
      MijMDI[23][28] = MijMDI[23][29] = 1.0 / 6;
  MijMDI[23][30] = MijMDI[23][33] = deletionRate / 2;

  genotypeArrMDI[23][0] = 0;
  genotypeArrMDI[23][1] = 3;
  genotypeArrMDI[23][2] = INVALID;

  // 24, CC
  MijMDI[24][1] = insertionRate;
  MijMDI[24][21] = MijMDI[24][25] = MijMDI[24][26] = 1.0 / 3;
  MijMDI[24][31] = deletionRate;

  genotypeArrMDI[24][0] = 1;
  genotypeArrMDI[24][1] = 1;
  genotypeArrMDI[24][2] = INVALID;

  // 25, CG
  MijMDI[25][8] = MijMDI[25][11] = insertionRate / 2;
  MijMDI[25][21] = MijMDI[25][22] = MijMDI[25][24] = MijMDI[25][26] =
      MijMDI[25][27] = MijMDI[25][28] = 1.0 / 6;
  MijMDI[25][31] = MijMDI[25][32] = deletionRate / 2;

  genotypeArrMDI[25][0] = 1;
  genotypeArrMDI[25][1] = 2;
  genotypeArrMDI[25][2] = INVALID;

  // 26, CT
  MijMDI[26][9] = MijMDI[26][14] = insertionRate / 2;
  MijMDI[26][21] = MijMDI[26][23] = MijMDI[26][24] = MijMDI[26][25] =
      MijMDI[26][28] = MijMDI[26][29] = 1.0 / 6;
  MijMDI[26][31] = MijMDI[26][33] = deletionRate / 2;

  genotypeArrMDI[26][0] = 1;
  genotypeArrMDI[26][1] = 3;
  genotypeArrMDI[26][2] = INVALID;

  // 27, GG
  MijMDI[27][2] = insertionRate;
  MijMDI[27][22] = MijMDI[27][25] = MijMDI[27][28] = 1.0 / 3;
  MijMDI[27][32] = deletionRate;

  genotypeArrMDI[27][0] = 2;
  genotypeArrMDI[27][1] = 2;
  genotypeArrMDI[27][2] = INVALID;

  // 28, GT
  MijMDI[28][12] = MijMDI[28][15] = insertionRate / 2;
  MijMDI[28][22] = MijMDI[28][23] = MijMDI[28][25] = MijMDI[28][26] =
      MijMDI[28][27] = MijMDI[28][29] = 1.0 / 6;
  MijMDI[28][32] = MijMDI[28][33] = deletionRate / 2;

  genotypeArrMDI[28][0] = 2;
  genotypeArrMDI[28][1] = 3;
  genotypeArrMDI[28][2] = INVALID;

  // 29, TT
  MijMDI[29][3] = insertionRate;
  MijMDI[29][23] = MijMDI[29][26] = MijMDI[29][28] = 1.0 / 3;
  MijMDI[29][33] = deletionRate;

  genotypeArrMDI[29][0] = 3;
  genotypeArrMDI[29][1] = 3;
  genotypeArrMDI[29][2] = INVALID;

  // 30, A-
  MijMDI[30][20] = insertionRate / 2;
  MijMDI[30][31] = MijMDI[30][32] = MijMDI[30][33] = 1.0 / 6;
  if (constrained == NO) {
    MijMDI[30][34] = deletionRate / 2;
  }

  genotypeArrMDI[30][0] = 0;
  genotypeArrMDI[30][1] = DELETED;
  genotypeArrMDI[30][2] = INVALID;

  // 31, C-
  MijMDI[31][24] = insertionRate / 2;
  MijMDI[31][30] = MijMDI[31][32] = MijMDI[31][33] = 1.0 / 6;
  if (constrained == NO) {
    MijMDI[31][34] = deletionRate / 2;
  }

  genotypeArrMDI[31][0] = 1;
  genotypeArrMDI[31][1] = DELETED;
  genotypeArrMDI[31][2] = INVALID;

  // 32, G-
  MijMDI[32][27] = insertionRate / 2;
  MijMDI[32][30] = MijMDI[32][31] = MijMDI[32][33] = 1.0 / 6;
  if (constrained == NO) {
    MijMDI[32][34] = deletionRate / 2;
  }

  genotypeArrMDI[32][0] = 2;
  genotypeArrMDI[32][1] = DELETED;
  genotypeArrMDI[32][2] = INVALID;

  // 33, T-
  MijMDI[33][29] = insertionRate / 2;
  MijMDI[33][30] = MijMDI[33][31] = MijMDI[33][32] = 1.0 / 6;
  if (constrained == NO) {
    MijMDI[33][34] = deletionRate / 2;
  }

  genotypeArrMDI[33][0] = 3;
  genotypeArrMDI[33][1] = DELETED;
  genotypeArrMDI[33][2] = INVALID;

  // 34, -
  if (constrained == NO) {
    genotypeArrMDI[34][0] = genotypeArrMDI[34][1] = DELETED;
  } else {
    genotypeArrMDI[34][0] = genotypeArrMDI[34][1] = INVALID;
  }
  genotypeArrMDI[34][2] = INVALID;

  for (int i = 0; i < 35; i++) {
    for (int j = 0; j < 35; j++) {
      if (i != j) {
        MijMDI[i][i] -= MijMDI[i][j];
      }
    }
  }

  InitGenotypeHashMapMDI();
}

/**
 * Get the base name of a folder.
 * E.g., .../path/to/me/, where 'me' is a folder, then 'me' will be returned.
 *
 * @param path a path
 * @return base name of a folder
 */
char *GetFolderBaseName(char *path) {
  int separatorIndex = -1, i = 0;
  while (1) {
    if (i > MAX_NAME) {
      fprintf(stderr, "Error! The path doesn't end with '%c'!", '\0');
      exit(-1);
    } else {
      if (path[i] == '\0')
        break;
      else if (path[i] == '/')
        separatorIndex = i;
    }
    i++;
  }

  if (separatorIndex == -1)
    return path;
  else
    return path + separatorIndex + 1;
}

int GetMaxInt(const int *arr, size_t size) {
  int max = arr[0];

  for (size_t i = 1; i < size; i++)
    if (max < arr[i])
      max = arr[i];

  return max;
}

int GetSumInt(const int *arr, size_t size) {
  int sum = 0;

  for (size_t i = 0; i < size; i++)
    sum += arr[i];

  return sum;
}

/* determine whether to make a coalescence tree or not */
int DoCoalescenceTree(int numDataSets, int replicate, int fixedTree,
                      int noisy) {
  if (numDataSets == 1) {
    return 1;
  } else {
    if (replicate == 0) {
      return 1;
    } else {
      if (fixedTree == 1) {
        if (noisy > 2) {
          fprintf(stderr,
                  "\n Use the fixed tree which has been coalesced before.");
        }

        return 0;
      } else {
        return 1;
      }
    }
  }
}

int GetExistingAllelesNrMD(int code) {
  if (code >= 0 && code < 10) {
    return 2;
  } else if (code < 14 && code >= 10) {
    return 1;
  } else if (code == 14) {
    return 0;
  } else {
    return -1;
  }
}

int GetExistingAllelesNrMDI(int code) {

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // AAA   CCC   GGG   TTT   AAC   AAG   AAT   ACC   CCG   CCT   AGG   CGG

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // GGT   ATT   CTT   GTT   ACG   ACT   AGT   CGT   AA    AC    AG    AT

  // 24    25    26    27    28    29    30    31    32    33    34
  // CC    CG    CT    GG    GT    TT    A-    C-    G-    T-    -

  if (code >= 0 && code < 20) {
    return 3;
  } else if (code < 30 && code >= 20) {
    return 2;
  } else if (code < 34 && code >= 30) {
    return 1;
  } else if (code == 34) {
    return 0;
  } else {
    return -1;
  }
}

int WhichRefGenotype(int nucleotide, int MuInDelModel) {
  if (nucleotide == A) {
    if (MuInDelModel == YES) {
      return 20;
    } else {
      return 0;
    }
  } else if (nucleotide == C) {
    if (MuInDelModel == YES) {
      return 24;
    } else {
      return 4;
    }
  } else if (nucleotide == G) {
    if (MuInDelModel == YES) {
      return 27;
    } else {
      return 7;
    }
  } else if (nucleotide == T) {
    if (MuInDelModel == YES) {
      return 29;
    } else {
      return 9;
    }
  } else {
    fprintf(stderr, "\nERROR in WhichRefGenotype: nucleotide = %c\n",
            nucleotide);
    exit(-1);
  }
}

int NumOverlap(const int *paGenotype, const int *chGenotype, int maxCopies,
               int paStart, int chStart, int *overlap) {
  int num = 0, flag_1, flag_2;
  int localOverlap[maxCopies];
  int *tmp;

  if (overlap == NULL) {
    tmp = localOverlap;
  } else {
    tmp = overlap;
  }

  for (int i = 0; i < maxCopies; i++) {
    tmp[i] = -1;
  }

  for (int i = chStart; i < maxCopies; ++i) {
    flag_1 = 0;
    if (chGenotype[i] >= 0) {
      for (int j = paStart; j < maxCopies; ++j) {
        flag_2 = 0;
        if (paGenotype[j] >= 0) {
          if (chGenotype[i] == paGenotype[j]) {
            for (int k = 0; k < maxCopies; ++k) {
              if (tmp[k] == -1) {
                tmp[k] = j;
                num++;
                flag_1 = 1;
                break;
              } else {
                if (tmp[k] == j) {
                  flag_2 = 1;
                  break;
                } else {
                  continue;
                }
              }
            }
            if (flag_1 == 1)
              break;
            if (flag_2 == 1)
              continue;
          }
        }
      }
    }
  }

  return num;
}

void AllocateThreads(int threadsNr) {
  threads = (pthread_t *)malloc(threadsNr * sizeof(pthread_t));
  if (!threads) {
    fprintf(stderr, "Could not allocate the threads structure\n");
    exit(-1);
  }
}

void InitThreadsStructure(int cellNr, int siteNr, int threadsNr) {
  threadStr = (ThreadStr *)malloc((threadsNr + 1) * sizeof(ThreadStr));
  if (!threadStr) {
    fprintf(stderr, "Could not allocate the threadStr structure\n");
    exit(-1);
  }
  sitesPointsForThreads = (int *)calloc(threadsNr + 2, sizeof(int));
  if (!sitesPointsForThreads) {
    fprintf(stderr, "Could not allocate the sitesPointsForThreads structure\n");
    exit(-1);
  }

  int nrSitesPerThread = siteNr / threadsNr;
  int start = 0, end = nrSitesPerThread;
  for (long i = 0; i < threadsNr; i++) {
    threadStr[i].id = i;
    threadStr[i].start = start;
    threadStr[i].end = i == threadsNr - 1 ? siteNr : end;
    threadStr[i].numSubSites = threadStr[i].end - threadStr[i].start;
    threadStr[i].rd = GetRandDevElement(randDevs, i);

    /* allocate genotype data to be stored in data[cell][site] */
    /* add 1 more size when allocating to make sure free() correctly works */
    /* only allocate space for the first dataset, and release it after all
     * dataset processed to save space and time */
    threadStr[i].data = (DataStr **)calloc(2 * cellNr + 2, sizeof(DataStr *));
    if (!threadStr[i].data) {
      fprintf(stderr, "Could not allocate the data structure in thread %ld\n",
              i);
      exit(-1);
    }

    for (int j = 0; j < 2 * cellNr + 1; j++) {
      threadStr[i].data[j] =
          (DataStr *)calloc(threadStr[i].numSubSites + 1, sizeof(DataStr));
      if (!threadStr[i].data[j]) {
        fprintf(stderr,
                "Could not allocate the data[] structure in thread %ld for "
                "cell %d\n",
                i, j);
        exit(-1);
      }

      for (int k = 0; k < threadStr[i].numSubSites; k++) {
        threadStr[i].data[j][k].pastEvolutionaryEvents =
            (int *)calloc(5, sizeof(unsigned int));
        if (!threadStr[i].data[j][k].pastEvolutionaryEvents) {
          fprintf(stderr,
                  "Could not allocate the data[][].pastEvolutionaryEvents"
                  " structure in thread %ld for cell %d and site %d\n",
                  i, j, k);
          exit(-1);
        }
      }
    }

    threadStr[i].hasParallelMutations =
        (int *)calloc(threadStr[i].numSubSites + 1, sizeof(int));
    threadStr[i].parallelMutations =
        (unsigned **)calloc(threadStr[i].numSubSites + 1, sizeof(unsigned *));
    threadStr[i].summary =
        (int **)calloc(threadStr[i].numSubSites + 1, sizeof(int *));
    for (int j = 0; j < threadStr[i].numSubSites; j++) {
      threadStr[i].summary[j] = (int *)calloc(5, sizeof(int));
      threadStr[i].parallelMutations[j] =
          (unsigned *)calloc(cellNr + 2, sizeof(unsigned));
    }

    threadStr[i].isSNV =
        (int *)calloc(threadStr[i].numSubSites + 1, sizeof(int));
    threadStr[i].isDeleted =
        (int *)calloc(threadStr[i].numSubSites + 1, sizeof(int));
    threadStr[i].isInserted =
        (int *)calloc(threadStr[i].numSubSites + 1, sizeof(int));

    sitesPointsForThreads[i + 1] = threadStr[i].end;

    start = end;
    end += nrSitesPerThread;
  }
}

void AllocateRawDataPrinters(int threadsNr)
{
  rawDataPrinters = (ThreadedRawDataPrinter *)malloc(
      (threadsNr + 1) * sizeof(ThreadedRawDataPrinter));

  for (long i = 0; i < threadsNr; i++)
  {
    rawDataPrinters[i].id = i;
    rawDataPrinters[i].rd = GetRandDevElement(randDevs, i);
  }
}

void InitRawDataPrinters(int tumCellNum, int threadsNr) {
  int numCellsPerThread = 0;

  if (tumCellNum + 1 > threadsNr)
    numCellsPerThread = (tumCellNum + 1) / threadsNr;
  else
    numCellsPerThread = 1;

  int numCellsLeft = tumCellNum + 1;
  int cellIndex = 0;
  for (long i = 0; i < threadsNr; i++) {

    if (i < threadsNr - 1) {
      if (numCellsLeft > 0 && numCellsLeft >= numCellsPerThread)
        rawDataPrinters[i].numCells = numCellsPerThread;
      else if (numCellsLeft > 0 && numCellsLeft < numCellsPerThread)
        rawDataPrinters[i].numCells = numCellsLeft;
      else
        rawDataPrinters[i].numCells = 0;
    } else if (i == threadsNr - 1) {
      if (numCellsLeft > 0)
        rawDataPrinters[i].numCells = numCellsLeft;
      else
        rawDataPrinters[i].numCells = 0;
    }

    if (rawDataPrinters[i].numCells > 0) {

      rawDataPrinters[i].cellIndices = (unsigned int *)malloc(
          (rawDataPrinters[i].numCells + 1) * sizeof(unsigned int));
      rawDataPrinters[i].files = (RawDataFile *)malloc(
          (rawDataPrinters[i].numCells + 1) * sizeof(RawDataFile));

      for (int j = 0; j < rawDataPrinters[i].numCells; j++) {
        rawDataPrinters[i].cellIndices[j] = cellIndex + j;

        RawDataFile tmp = rawDataFiles[rawDataPrinters[i].cellIndices[j] + 1];
        rawDataPrinters[i].files[j].fp = tmp.fp;
        strcpy(rawDataPrinters[i].files[j].fileName, tmp.fileName);
      }

    } else {
      rawDataPrinters[i].cellIndices = NULL;
      rawDataPrinters[i].files = NULL;
    }

    numCellsLeft -= rawDataPrinters[i].numCells;
    if (numCellsLeft < 0)
      numCellsLeft = 0;

    cellIndex += rawDataPrinters[i].numCells;
  }
}

int GetThreadIndex(int siteIndex) {
  for (int i = 0; i < numThreads; i++) {
    if (sitesPointsForThreads[i] <= siteIndex &&
        siteIndex < sitesPointsForThreads[i + 1]) {
      return i;
    }
  }

  fprintf(stderr,
          "Could not find the index of thread (from 0 - %d) according "
          "to the given site %d \n",
          sitesPointsForThreads[numThreads], siteIndex);
  exit(-1);
}

void AllocateCellNamesExtd(int cellNr) {
  cellNamesExtd = (char **)calloc(2 * cellNr + 2, sizeof(char *));
  if (!cellNamesExtd) {
    fprintf(stderr, "Could not allocate tipNames (%ld)\n",
            (2 * cellNr + 2) * sizeof(char *));
    exit(-1);
  }
  for (int i = 0; i < 2 * cellNr + 1; i++) {
    cellNamesExtd[i] = (char *)calloc(MAX_NAME, sizeof(char));
    if (!cellNamesExtd[i]) {
      fprintf(stderr, "Could not allocate tipNames[%d] (%ld)\n", i,
              (MAX_NAME) * sizeof(char));
      exit(-1);
    }

    if (i > cellNr)
      sprintf(cellNamesExtd[i], "tumcell%04d", i);
  }
}

void AllocateCellStructureExtd(int cellNr, int siteNr) {
  cellExtd = (CellStrExtd *)calloc(2 * cellNr + 2, sizeof(CellStrExtd));
  if (!cellExtd) {
    fprintf(stderr, "Could not allocate cellExtd (%ld)\n",
            (2 * cellNr + 1) * sizeof(CellStrExtd));
    exit(-1);
  }

  for (int i = 0; i < 2 * cellNr + 1; ++i) {
    cellExtd[i].site =
        (CellSiteStrExtd *)calloc(siteNr + 1, sizeof(CellSiteStrExtd));
    if (!cellExtd[i].site) {
      fprintf(stderr, "Could not allocate the cellExtd[%d].site structure\n",
              i);
      exit(-1);
    }

    cellExtd[i].isDoublet = NO;
    cellExtd[i].mixeTumorCellIdx = -1;
    cellExtd[i].mixedTumorCellNameIdx = -1;

    for (int j = 0; j < siteNr; ++j) {
      cellExtd[i].site[j].hasADO = NO;
      cellExtd[i].site[j].numADO = 0;
      cellExtd[i].site[j].maxNumExistingAllelesPerStrand = 1;

      cellExtd[i].site[j].readCount = (int *)calloc(5, sizeof(int));
      if (!cellExtd[i].site[j].readCount) {
        fprintf(stderr,
                "Could not allocate the cellExtd[%d].site[%d].readCount "
                "structure\n",
                i, j);
        exit(-1);
      }

      if (i < cellNr + 1)
      {
        cellExtd[i].site[j].observedAlleles = (int *)calloc(5, sizeof(int));
        if (!cellExtd[i].site[j].observedAlleles) {
          fprintf(stderr,
                  "Could not allocate the cellExtd[%d].site[%d].observedAlleles "
                  "structure\n",
                  i, j);
          exit(-1);
        }
      }
    }
  }
}

int GetReadBlockIndex(int site) {
  for (int i = 0; i < numReadBlocks; i++) {
    if (site < startSiteBlocks[i + 1])
      return i;
  }

  fprintf(stderr, "Error!");
  exit(1);
}

void GetCellName(const char *fileName, char *out, int len) {
  for (int i = 0; i < len; i++) {
    if (fileName[i] == '\0' || fileName[i] == '.') {
      out[i] = '\0';
      break;
    } else
      out[i] = fileName[i];
  }
}

int ValueInArray(int value, const int *array, int size) {
  for (int i = 0; i < size; i++)
    if (array[i] == value)
      return 1;

  return 0;
}

int AddValueToArray(int value, int *array, int size, int marker) {
  for (int i = 0; i < size; i++)
    if (array[i] == marker) {
      array[i] = value;
      return 0;
    }

  return 1;
}

void CountMutationNum(const int *chAlleles, int startPos, int maxCopies,
                      int referenceNuc, const int *paAlleles,
                      const int *overlap, int *mutationCounter) {
  int exOverlapPos[maxCopies], compareFlag;

  for (int i = 0; i < maxCopies; i++)
    exOverlapPos[i] = -1;

  for (int i = startPos; i < maxCopies; i++) {
    compareFlag = 1;

    if (chAlleles[i] < 0)
      continue;

    for (int j = 0; j < maxCopies && overlap[j] >= 0; j++) {
      if (chAlleles[i] == paAlleles[overlap[j]]) {
        if (ValueInArray(j, exOverlapPos, maxCopies) == 0) {
          if (AddValueToArray(j, exOverlapPos, maxCopies, -1) == 0) {
            compareFlag = 0;
            break;
          } else {
            fprintf(stderr,
                    "ERROR! Cannot add %d in an array of length "
                    "%d as it is already full: [",
                    j, maxCopies);
            for (int k = 0; k < maxCopies; k++) {
              fprintf(stderr, "%d", exOverlapPos[k]);

              if (k < maxCopies - 1)
                fprintf(stderr, ",");
              else
                fprintf(stderr, "]\n");
            }
            exit(1);
          }
        } else
          continue;
      }
    }

    if (compareFlag == 1) {
      mutationCounter[0]++;

      if (chAlleles[i] == referenceNuc)
        mutationCounter[1]++;
    }
  }
}

/*
void GetMuDelInsNr(ThreadStr *thStr, int anccell, int cell, int site,
                   int maxCopies, int referenceNuc) {
  int overlapNr, overlap[maxCopies];

  if (thStr->data[cell][site].numExistingAlleles !=
      thStr->data[anccell][site].numExistingAlleles) {
    if (thStr->data[anccell][site].numExistingAlleles == 3) {
      // allSites[site].numDeletions++;
      // numDEL++;
      thStr->summary[site][3]++;
      thStr->data[cell][site].pastEvolutionaryEvents[3]++;

      if (thStr->data[cell][site].numExistingAlleles == 2) {
        if (thStr->data[cell][site].alleles[0] < 0) {
          // e.g., A/AA -> -/AA
          overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                                 thStr->data[cell][site].alleles, maxCopies, 1,
                                 1, overlap);
          // allSites[site].numMutations +=
          // (thStr->data[cell][site].numExistingAlleles - overlapNr); numMU +=
          // (thStr->data[cell][site].numExistingAlleles - overlapNr);
          thStr->summary[site][0] +=
              (thStr->data[cell][site].numExistingAlleles - overlapNr);
          CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                           referenceNuc, thStr->data[anccell][site].alleles,
                           overlap,
                           thStr->data[cell][site].pastEvolutionaryEvents);
        } else {
          overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                                 thStr->data[cell][site].alleles, maxCopies, 0,
                                 0, overlap);
          if (overlapNr == 2) {
            if (thStr->data[anccell][site].alleles[0] !=
                    thStr->data[cell][site].alleles[0] &&
                thStr->data[anccell][site].alleles[0] !=
                    thStr->data[cell][site].alleles[1]) {
              // e.g., A/CG -> C/G
              // Ambiguity will raise.
              // allSites[site].numMutations++;
              // numMU++;
              thStr->summary[site][0]++;
              thStr->data[cell][site].pastEvolutionaryEvents[0]++;
              if (ValueInArray(referenceNuc, thStr->data[cell][site].alleles,
                               2) == 1)
                thStr->data[cell][site].pastEvolutionaryEvents[2]++;
            }
          } else {
            // allSites[site].numMutations +=
            // (thStr->data[cell][site].numExistingAlleles - overlapNr); numMU
            // += (thStr->data[cell][site].numExistingAlleles - overlapNr);
            thStr->summary[site][0] +=
                (thStr->data[cell][site].numExistingAlleles - overlapNr);
            CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        }
      } else if (thStr->data[cell][site].numExistingAlleles == 1) {
        // allSites[site].numDeletions++;
        // numDEL++;
        thStr->summary[site][3]++;
        thStr->data[cell][site].pastEvolutionaryEvents[3]++;

        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 0,
                               overlap);
        if (overlapNr == 0) {
          // allSites[site].numMutations++;
          // numMU++;
          thStr->summary[site][0]++;
          CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                           referenceNuc, thStr->data[anccell][site].alleles,
                           overlap,
                           thStr->data[cell][site].pastEvolutionaryEvents);
        }
      } else {
        // allSites[site].numDeletions += 2;
        // numDEL += 2;
        thStr->summary[site][3] += 2;
        thStr->data[cell][site].pastEvolutionaryEvents[3] += 2;
      }
    } else if (thStr->data[anccell][site].numExistingAlleles == 2) {
      if (thStr->data[cell][site].numExistingAlleles == 3) {
        if (thStr->data[anccell][site].alleles[0] >= 0) {
          // allSites[site].numInsertions++;
          // numINS++;
          thStr->summary[site][4]++;
          thStr->data[cell][site].pastEvolutionaryEvents[4]++;

          overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                                 thStr->data[cell][site].alleles, maxCopies, 0,
                                 0, overlap);
          if (thStr->data[anccell][site].alleles[0] ==
              thStr->data[anccell][site].alleles[1]) {
            // e.g., G/G
            if (overlapNr == 2) {
              if (thStr->data[cell][site].alleles[0] !=
                      thStr->data[cell][site].alleles[1] ||
                  thStr->data[cell][site].alleles[0] !=
                      thStr->data[cell][site].alleles[2] ||
                  thStr->data[cell][site].alleles[1] !=
                      thStr->data[cell][site].alleles[2]) {
                // e.g., G/G -> C/GG or G/CG
                // allSites[site].numMutations++;
                // numMU++;
                thStr->summary[site][0]++;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            } else if (overlapNr == 1) {
              if (thStr->data[cell][site].alleles[1] ==
                  thStr->data[cell][site].alleles[2]) {
                // e.g., G/G -> G/AA
                // allSites[site].numMutations++;
                // numMU++;
                thStr->summary[site][0]++;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, 2, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              } else {
                // e.g., G/G -> A/CG
                // allSites[site].numMutations += 2;
                // numMU += 2;
                thStr->summary[site][0] += 2;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            } else {
              if (thStr->data[cell][site].alleles[1] ==
                  thStr->data[cell][site].alleles[2]) {
                // e.g., G/G -> A/CC
                // allSites[site].numMutations += 2;
                // numMU += 2;
                thStr->summary[site][0] += 2;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, 2, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              } else {
                // e.g., G/G -> A/CT
                // allSites[site].numMutations += 3;
                // numMU += 3;
                thStr->summary[site][0] += 3;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            }
          } else {
            // e.g., C/T
            if (overlapNr == 2) {
              if (thStr->data[cell][site].alleles[1] !=
                  thStr->data[cell][site].alleles[2]) {
                if ((thStr->data[cell][site].alleles[0] !=
                     thStr->data[anccell][site].alleles[0]) &&
                    (thStr->data[cell][site].alleles[0] !=
                     thStr->data[anccell][site].alleles[1])) {
                  // e.g., C/T -> A/CT
                  // Ambiguity could raise.
                  // allSites[site].numMutations += 2;
                  // numMU += 2;
                  thStr->summary[site][0] += 2;
                  thStr->data[cell][site].pastEvolutionaryEvents[0] += 2;
                  if (referenceNuc == thStr->data[cell][site].alleles[0])
                    thStr->data[cell][site].pastEvolutionaryEvents[1] += 1;
                  if (referenceNuc == thStr->data[cell][site].alleles[1] ||
                      referenceNuc == thStr->data[cell][site].alleles[2])
                    thStr->data[cell][site].pastEvolutionaryEvents[2] += 1;
                } else {
                  // e.g., C/T -> C/AT
                  // allSites[site].numMutations++;
                  // numMU++;
                  thStr->summary[site][0]++;
                  CountMutationNum(
                      thStr->data[cell][site].alleles, 0, maxCopies,
                      referenceNuc, thStr->data[anccell][site].alleles, overlap,
                      thStr->data[cell][site].pastEvolutionaryEvents);
                }
              }
            } else if (overlapNr == 1) {
              if (thStr->data[cell][site].alleles[1] !=
                  thStr->data[cell][site].alleles[2]) {
                // e.g., C/T -> C/AG
                // allSites[site].numMutations += 2;
                // numMU += 2;
                thStr->summary[site][0] += 2;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              } else {
                // e.g., C/T -> A/TT
                // allSites[site].numMutations++;
                // numMU++;
                thStr->summary[site][0]++;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, 2, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            } else {
              if (thStr->data[cell][site].alleles[1] !=
                  thStr->data[cell][site].alleles[2]) {
                // e.g., C/T -> G/AG
                // allSites[site].numMutations += 3;
                // numMU += 3;
                thStr->summary[site][0] += 3;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              } else {
                // e.g., C/T -> G/AA
                // allSites[site].numMutations += 2;
                // numMU += 2;
                thStr->summary[site][0] += 2;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, 2, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            }
          }
        }
      } else if (thStr->data[cell][site].numExistingAlleles == 1) {
        // allSites[site].numDeletions++;
        // numDEL++;
        thStr->summary[site][3]++;
        thStr->data[cell][site].pastEvolutionaryEvents[3]++;

        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 0,
                               overlap);
        if (overlapNr == 0) {
          // e.g., C/T -> A/-
          // allSites[site].numMutations++;
          // numMU++;
          thStr->summary[site][0]++;
          CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                           referenceNuc, thStr->data[anccell][site].alleles,
                           overlap,
                           thStr->data[cell][site].pastEvolutionaryEvents);
        }
      } else {
        // allSites[site].numDeletions += 2;
        // numDEL += 2;
        thStr->summary[site][3] += 2;
        thStr->data[cell][site].pastEvolutionaryEvents[3] += 2;
      }
    } else if (thStr->data[anccell][site].numExistingAlleles == 1) {
      // e.g., C/-
      if (thStr->data[cell][site].numExistingAlleles == 2 &&
          thStr->data[cell][site].alleles[0] == DELETED) {
        // allSites[site].numInsertions++;
        // numINS++;
        thStr->summary[site][4]++;
        thStr->data[cell][site].pastEvolutionaryEvents[4]++;

        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 0,
                               overlap);
        if (overlapNr == 0) {
          if (thStr->data[cell][site].alleles[1] ==
              thStr->data[cell][site].alleles[2]) {
            // e.g., C/- -> -/AA
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            CountMutationNum(thStr->data[cell][site].alleles, 1, 2,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          } else {
            // e.g., C/- -> -/GT
            // allSites[site].numMutations += 2;
            // numMU += 2;
            thStr->summary[site][0] += 2;
            CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        } else {
          // e.g., C/- -> -/CT
          if (thStr->data[cell][site].alleles[1] !=
              thStr->data[cell][site].alleles[2]) {
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        }
      } else if (thStr->data[cell][site].numExistingAlleles == 0) {
        // allSites[site].numDeletions++;
        // numDEL++;
        thStr->summary[site][3]++;
        thStr->data[cell][site].pastEvolutionaryEvents[3]++;
      }
    }
  } else {
    if (thStr->data[anccell][site].numExistingAlleles == 3) {
      if (thStr->data[anccell][site].alleles[0] !=
          thStr->data[cell][site].alleles[0]) {
        // e.g., A/CC -> G/??
        // allSites[site].numMutations++;
        // numMU++;
        thStr->summary[site][0]++;
        thStr->data[cell][site].pastEvolutionaryEvents[0]++;
        if (referenceNuc == thStr->data[cell][site].alleles[0])
          thStr->data[cell][site].pastEvolutionaryEvents[1]++;
      }

      overlapNr =
          NumOverlap(thStr->data[anccell][site].alleles,
                     thStr->data[cell][site].alleles, maxCopies, 1, 1, overlap);
      // allSites[site].numMutations +=
      // (thStr->data[anccell][site].numExistingAlleles - 1 - overlapNr); numMU
      // += (thStr->data[anccell][site].numExistingAlleles - 1 - overlapNr);
      thStr->summary[site][0] +=
          (thStr->data[anccell][site].numExistingAlleles - 1 - overlapNr);
      CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                       referenceNuc, thStr->data[anccell][site].alleles,
                       overlap, thStr->data[cell][site].pastEvolutionaryEvents);
    } else if (thStr->data[anccell][site].numExistingAlleles == 1) {
      // allSites[site].numMutations++;
      // numMU++;
      thStr->summary[site][0]++;
      CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                       referenceNuc, thStr->data[anccell][site].alleles,
                       overlap, thStr->data[cell][site].pastEvolutionaryEvents);
    } else {
      if (thStr->data[anccell][site].alleles[0] != DELETED &&
          thStr->data[cell][site].alleles[0] == DELETED) {
        // allSites[site].numDeletions++;
        // numDEL++;
        thStr->summary[site][3]++;
        thStr->data[cell][site].pastEvolutionaryEvents[3]++;
        // allSites[site].numInsertions++;
        // numINS++;
        thStr->summary[site][4]++;
        thStr->data[cell][site].pastEvolutionaryEvents[4]++;

        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 1,
                               overlap);
        if (overlapNr == 2) {
          if (thStr->data[anccell][site].alleles[0] !=
              thStr->data[anccell][site].alleles[1]) {
            // e.g., G/T -> -/GT
            // Ambiguity might raise.
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            thStr->data[cell][site].pastEvolutionaryEvents[0]++;
            if (ValueInArray(referenceNuc, thStr->data[cell][site].alleles,
                             3) == 1)
              thStr->data[cell][site].pastEvolutionaryEvents[2]++;
          }
        } else if (overlapNr == 1) {
          if (thStr->data[cell][site].alleles[1] !=
              thStr->data[cell][site].alleles[2]) {
            // e.g., G/T -> -/AT
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        } else {
          if (thStr->data[cell][site].alleles[1] ==
              thStr->data[cell][site].alleles[2]) {
            // e.g., G/T -> -/AA
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            CountMutationNum(thStr->data[cell][site].alleles, 1, 2,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          } else {
            // e.g., G/T -> -/AC
            // allSites[site].numMutations += 2;
            // numMU += 2;
            thStr->summary[site][0] += 2;
            CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        }
      } else if (!(thStr->data[anccell][site].alleles[0] == DELETED &&
                   thStr->data[cell][site].alleles[0] != DELETED)) {
        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 0,
                               overlap);
        // allSites[site].numMutations +=
        // (thStr->data[anccell][site].numExistingAlleles - overlapNr); numMU +=
        // (thStr->data[anccell][site].numExistingAlleles - overlapNr);
        thStr->summary[site][0] +=
            (thStr->data[anccell][site].numExistingAlleles - overlapNr);
        CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                         referenceNuc, thStr->data[anccell][site].alleles,
                         overlap,
                         thStr->data[cell][site].pastEvolutionaryEvents);
      }
    }
  }
}
*/

void GetMuDelInsNr(ThreadStr *thStr, int anccell, int cell, int site,
                   int maxCopies, int referenceNuc) {
  int overlap[maxCopies], cntDiffNr, muCnt[2] = {0};

  cntDiffNr = thStr->data[anccell][site].numExistingAlleles -
              thStr->data[cell][site].numExistingAlleles;

  if (cntDiffNr > 0) {
    thStr->summary[site][2] += cntDiffNr;
    thStr->data[cell][site].pastEvolutionaryEvents[2] += cntDiffNr;
  } else if (cntDiffNr < 0) {
    thStr->summary[site][3] += abs(cntDiffNr);
    thStr->data[cell][site].pastEvolutionaryEvents[3] += abs(cntDiffNr);
  }

  NumOverlap(thStr->data[anccell][site].alleles,
             thStr->data[cell][site].alleles, maxCopies, 0, 0, overlap);
  CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                   thStr->data[anccell][site].alleles, overlap, muCnt);
  thStr->summary[site][0] += muCnt[0];
  thStr->summary[site][1] += muCnt[1];
  thStr->data[cell][site].pastEvolutionaryEvents[0] += muCnt[0];
  thStr->data[cell][site].pastEvolutionaryEvents[1] += muCnt[1];
}

int GetAllele(double randomNr, int num) {
  int alleleNr = -1, intArray[num];

  for (int i = 0; i < num; ++i) {
    intArray[i] = i + 1;
  }

  for (int i = 0; i < num; ++i) {
    if ((randomNr * num) < intArray[i]) {
      alleleNr = intArray[i] - 1;
      break;
    }
  }

  return alleleNr;
}

void GetGenotypeStructure(const int *alleles, int *out) {
  for (int i = 0; i < 4; ++i)
    out[i] = 0;

  for (int i = 0; i < 3; ++i)
    if (alleles[i] >= 0)
      out[alleles[i]]++;
}

/* ******************** GetShapeController ******************** */
int GetShapeController(const int *alleles, int numExistingAlleles) {
  if (numExistingAlleles == 0) {
    fprintf(stderr,
            "Could not get shape controller value for 0 existing alleles.\n");
    exit(-1);
  }

  int result = -1, maxAlleles, genotypeStructure[4];

  GetGenotypeStructure(alleles, genotypeStructure);
  maxAlleles = GetMaxInt(genotypeStructure, 4);

  if (maxAlleles == numExistingAlleles) {
    result = 0;
  } else if (maxAlleles * 2 == numExistingAlleles) {
    result = 1;
  } else {
    result = 1;
  }

  if (result == -1) {
    fprintf(stderr, "Error!");
    exit(-1);
  }

  return result;
}

/******************** ProcessZeroCoverageAndZeroProportion *******************/
void ProcessZeroCoverageAndZeroProportion(int cell, int site,
                                          int MLVariantNucleotide,
                                          int numDataSets) {
  int coverageFlag;

  /* analyze the reason for zero coverage in negative binomial distribution */
  if (cellExtd[cell].site[site].numReads == 0) /* if there is zero coverage */
  {
    coverageFlag = 0;
    zeroCoverage++;
    if (fixedTree == YES && numDataSets > 1) {
      zeroCoverageAll++;
    }

    if (cellExtd[cell].site[site].trueNumExistingAlleles ==
        0) /* if there is complete deletions */
    {
      coverageFlag = 1;
      completeDeletions++;
      if (fixedTree == YES && numDataSets > 1) {
        completeDeletionsAll++;
      }
    } else /* if there are alleles after evolving along the tree */
    {
      if (cellExtd[cell].site[site].observedNumExistingAlleles ==
          0) /* if there is complete ADO */
      {
        coverageFlag = 1;
        completeADOs++;
        if (fixedTree == YES && numDataSets > 1) {
          completeADOsAll++;
        }
      } else /* if there are alleles after evolution and ADO */
      {
        coverageFlag = 1;
        noCoverage++;
        allelicZeroCoverageProportion
            [cellExtd[cell].site[site].observedNumExistingAlleles][1]++;
        if (fixedTree == YES && numDataSets > 1) {
          noCoverageAll++;
        }
      }
    }

    if (coverageFlag == 0) {
      fprintf(
          stderr,
          "WARNING! Unknown reason for zero coverage at site %d in cell %d!",
          site, cell);
    }
  }

  /* analyze the reason for zero proportion in beta binomial distribution */
  if (cellExtd[cell].site[site].readCount[MLVariantNucleotide] ==
      0) /* if there is zero reads for the site-wise variant nucleotide */
  {
    zeroProportionVAF++;
    if (fixedTree == YES && numDataSets > 1)
      zeroProportionVAFAll++;

    if (cellExtd[cell].site[site].numReads >
        0) /* if there is positive coverage */
    {
      noReadsPositiveCoverage++;
      if (fixedTree == YES && numDataSets > 1)
        noReadsPositiveCoverageAll++;
    }
  }
}

char WhichNucExtd(int nucleotide, char missing) {
  if (nucleotide >= 0) {
    if (nucleotide == A)
      return 'A';
    else if (nucleotide == C)
      return 'C';
    else if (nucleotide == G)
      return 'G';
    else if (nucleotide == T)
      return 'T';
    else
      return ('N');
  }

  return missing;
}

/******************** ComputeMeanVarianceDoubleArr *******************/
void ComputeMeanVarianceDoubleArr(const double *arr, int size, double *meanAddr,
                                  double *varAddr, int unbiasedVar) {
  double sum, distSquared, mean, var;

  sum = distSquared = 0.0;

  for (int i = 0; i < size; ++i) {
    sum += arr[i];
  }

  mean = sum / size;

  for (int i = 0; i < size; ++i) {
    distSquared += pow(arr[i] - mean, 2);
  }

  if (unbiasedVar) {
    var = distSquared / (size - 1);
  } else {
    var = distSquared / size;
  }

  *meanAddr = mean;
  *varAddr = var;
}

void PrintGenotypeInLettersToFile(FILE *fp, int numExistingAlleles,
                                  const int *alleles, int phased, char missing,
                                  int withSpaceFlag) {
  if (numExistingAlleles == 3) {
    if (fp) {
      if (withSpaceFlag)
        fprintf(fp, " ");

      fprintf(fp, "%c%c%c", WhichNucExtd(alleles[0], MISSING_CHAR_ALT),
              WhichNucExtd(alleles[1], MISSING_CHAR_ALT),
              WhichNucExtd(alleles[2], MISSING_CHAR_ALT));
    }
  } else if (numExistingAlleles == 2 || numExistingAlleles == 1) {
    if (fp) {
      if (withSpaceFlag)
        fprintf(fp, " ");

      fprintf(fp, "%c%c%c", WhichNucExtd(alleles[0], MISSING_CHAR_ALT),
              phased > 0 ? '|' : '/',
              WhichNucExtd(alleles[1], MISSING_CHAR_ALT));
    }
  } else if (numExistingAlleles == 0) {
    if (fp) {
      if (withSpaceFlag)
        fprintf(fp, " ");

      fprintf(fp, "%c", missing);
    }
  }
}

void GetGenotypeInNumbersToFileGivenExistingAlleles(int referenceAllele,
                                                    int stopPosition,
                                                    const int *alleles,
                                                    char missing, int *altNucs,
                                                    int *nucCounts, int *out) {
  /* sanity check */
  if (stopPosition < 2 || stopPosition > 3) {
    fprintf(stderr,
            "ARGUMENT ERROR! Stop position should be either 2 or 3! "
            "(%d is given)\n",
            stopPosition);
    exit(-1);
  }

  for (int i = 0; i < stopPosition; ++i) {
    if (alleles[i] < 0) {
      nucCounts[4]++;
      out[i] = (int)missing;
    } else if (alleles[i] >= 0 && alleles[i] == referenceAllele) {
      nucCounts[0]++;
      out[i] = '0';
    } else {
      for (int j = 0; j < 3; j++) {
        if (altNucs[j] < 0) {
          altNucs[j] = alleles[i];
          nucCounts[j + 1]++;
          out[i] = j + 1 + '0';
          break;
        } else if (altNucs[j] == alleles[i]) {
          nucCounts[j + 1]++;
          out[i] = j + 1 + '0';
          break;
        }
      }
    }
  }
}

void PrintGenotypeInNumbers(FILE *fp, int stopPosition, int useSeparator,
                            char separator, const int *results) {
  for (int i = 0; i < stopPosition; i++) {
    if (useSeparator && i == 1)
      fprintf(fp, "%c", separator);

    fprintf(fp, "%c", results[i]);
  }
}

/*
 * -3: -
 * -2: 1/-
 * -1: 0/-
 * 0: 0/0
 * 1: 0/1
 * 2: 1/1
 * 3: 1/2 (namely 1/1')
 * 4: 000
 * 5: 001
 * 6: 011
 * 7: 012
 * 8: 111
 * 9: 112
 * 10: 123
 */
int GetGenotypeInTernaryMatrix(const int *nucCounts, const int *results) {
  int code;

  if (nucCounts[4] >= 2)
    code = -3;
  else if (nucCounts[4] == 1) {
    if (nucCounts[0] > 0)
      code = -1;
    else
      code = -2;
  } else {
    const int sum = GetSumInt(nucCounts, 4);
    const int maxCnt = GetMaxInt(nucCounts, 4);

    if (sum == 2) {
      if (nucCounts[0] == 2)
        code = 0;
      else if (nucCounts[0] == 1)
        code = 1;
      else {
        if (maxCnt == 2)
          code = 2;
        else if (maxCnt == 1)
          code = 3;
        else {
          fprintf(stderr, "Unrecognized genotype: %d%d%d", results[0],
                  results[1], results[2]);
          exit(-1);
        }
      }
    } else if (sum == 3) {
      if (nucCounts[0] == 3)
        code = 4;
      else if (nucCounts[0] == 2)
        code = 5;
      else if (nucCounts[0] == 1) {
        if (maxCnt == 1)
          code = 7;
        else if (maxCnt == 2)
          code = 6;
        else {
          fprintf(stderr, "Unrecognized genotype: %d%d%d", results[0],
                  results[1], results[2]);
          exit(-1);
        }
      } else {
        if (maxCnt == 3)
          code = 8;
        else if (maxCnt == 2)
          code = 9;
        else if (maxCnt == 1)
          code = 10;
        else {
          fprintf(stderr, "Unrecognized genotype: %d%d%d", results[0],
                  results[1], results[2]);
          exit(-1);
        }
      }
    } else {
      fprintf(stderr, "Unrecognized genotype: %d%d%d", results[0], results[1],
              results[2]);
      exit(-1);
    }
  }

  return code;
}

void PrintGenotypeInNumbersToFile(FILE *fpNr, FILE *fpMt, int referenceAllele,
                                  int numExistingAlleles, const int *alleles,
                                  int phased, char missing,
                                  int withWhiteSpaceFlag, int useWhiteSpace) {
  if (withWhiteSpaceFlag) {
    if (fpNr)
      fprintf(fpNr, "%c", useWhiteSpace == 1 ? ' ' : '\t');

    if (fpMt)
      fprintf(fpMt, "%c", useWhiteSpace == 1 ? ' ' : '\t');
  }

  int altNucs[] = {-1, -1, -1};
  // for ref, three alts in the order of appearance, and "missing" character
  int nucCounts[] = {0, 0, 0, 0, 0};
  int results[] = {'0', '0', '0'};

  if (numExistingAlleles == 3) {
    GetGenotypeInNumbersToFileGivenExistingAlleles(
        referenceAllele, 3, alleles, missing, altNucs, nucCounts, results);

    if (!phased)
      qsort(results, 3, sizeof(int), GreaterThanInt);

    if (fpNr)
      PrintGenotypeInNumbers(fpNr, 3, NO, phased > 0 ? '|' : '/', results);

    if (fpMt)
      fprintf(fpMt, "%d", GetGenotypeInTernaryMatrix(nucCounts, results));

  } else if (numExistingAlleles == 2 || numExistingAlleles == 1) {
    GetGenotypeInNumbersToFileGivenExistingAlleles(
        referenceAllele, 2, alleles, missing, altNucs, nucCounts, results);

    if (!phased)
      qsort(results, 2, sizeof(int),
            numExistingAlleles == 2 ? GreaterThanInt : LessThanInt);

    if (fpNr)
      PrintGenotypeInNumbers(fpNr, 2, YES, phased > 0 ? '|' : '/', results);

    if (fpMt)
      fprintf(fpMt, "%d", GetGenotypeInTernaryMatrix(nucCounts, results));

  } else if (numExistingAlleles == 0) {
    GetGenotypeInNumbersToFileGivenExistingAlleles(
        referenceAllele, 2, alleles, missing, altNucs, nucCounts, results);

    if (fpNr)
      fprintf(fpNr, "%c", missing);

    if (fpMt)
      fprintf(fpMt, "%d", GetGenotypeInTernaryMatrix(nucCounts, results));
  }
}

void FreeFiniteMuInDelSimulator(int size) {
  for (int i = 0; i < size; ++i) {
    struct NodeIntToIntArr *tmp = NULL;

    while (genotypeHashMapIntToIntArr->list[i]) {
      tmp = genotypeHashMapIntToIntArr->list[i];
      genotypeHashMapIntToIntArr->list[i] =
          genotypeHashMapIntToIntArr->list[i]->next;
      free(tmp);
      tmp = NULL;
    }
  }
  free(genotypeHashMapIntToIntArr->list);
  genotypeHashMapIntToIntArr->list = NULL;
  free(genotypeHashMapIntToIntArr);
  genotypeHashMapIntToIntArr = NULL;

  for (int i = 0; i < size; ++i) {
    struct NodeIntArrToInt *tmp = NULL;

    while (genotypeHashMapIntArrToInt->list[i]) {
      tmp = genotypeHashMapIntArrToInt->list[i];
      genotypeHashMapIntArrToInt->list[i] =
          genotypeHashMapIntArrToInt->list[i]->next;
      free(tmp);
      tmp = NULL;
    }
  }
  free(genotypeHashMapIntArrToInt->list);
  genotypeHashMapIntArrToInt->list = NULL;
  free(genotypeHashMapIntArrToInt);
  genotypeHashMapIntArrToInt = NULL;
}

int SumInt(const int *arr, unsigned int length) {
  int ret = 0;
  for (int i = 0; i < length; i++) {
    ret += arr[i];
  }
  return ret;
}

double SumDouble(const double *arr, unsigned int length) {
  double ret = 0;
  for (int i = 0; i < length; i++) {
    ret += arr[i];
  }
  return ret;
}

/* Ascending order */
int GreaterThanInt(const void *a, const void *b) {
  return (*(int *)a - *(int *)b);
}

/* Descending order */
int LessThanInt(const void *a, const void *b) {
  return (*(int *)b - *(int *)a);
}

int BinarySearch(const int *arr, int start, int end, int val) {
  if (start == end)
    return -1;

  int midPoint = start + (end - start) / 2;

  if (arr[midPoint] == val)
    return midPoint;
  else if (arr[midPoint] < val)
    return BinarySearch(arr, midPoint + 1, end, val);
  else
    return BinarySearch(arr, start, midPoint, val);
}

int GetTumCellsTotalNum(int cellNr)
{
  int num = cellNr;

  for (int i = 0; i < cellNr + 1; ++i)
    if (cellExtd[i].isDoublet == YES)
      num++;

  return num;
}

void FillMpileupEntry(
  const int refNuc,
  const CellSiteStrExtd *csite1, const CellSiteStrExtd *csite2,
  const int cov, char *reads, char *quals
  )
{
  if (csite1 == NULL && csite2 == NULL)
  {
    fprintf(stderr, "Invalid arguments.\n");
    return;
  }

  if (cov == 0) {
    reads[0] = '*';
    quals[0] = '*';
    reads[1] = '\0';
    quals[1] = '\0';
  } else if (cov > 0) {
    int l = 0;
    for (int k = 0; k < 4; k++) {
      int agg_read_count = 0;
      if (csite1 != NULL) agg_read_count += csite1->readCount[k];
      if (csite2 != NULL) agg_read_count += csite2->readCount[k];

      for (int m = 0; m < agg_read_count; m++) {
        if (k == refNuc)
          reads[l] = '.';
        else
          reads[l] = WhichNucExtd(k, MISSING_CHAR);

        quals[l] = '~';
        l++;
      }
    }

    reads[l] = '\0';
    quals[l] = '\0';
  }
}

#endif // MYSIMULATOR_FINITEMUINDELSIMULATOR_H
