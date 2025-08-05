/******************************************************************************
        GENERAL PUBLIC LICENSE
        This program is free software; you can redistribute it and/or
        modify it under the terms of the GNU General Public License
        as published by the Free Software Foundation; either version 2
        of the License, or (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program; if not, write to the Free Software
        Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
 USA.
 ******************************************************************************/

/*  cellcoal
//  CellCoal
//
//  Copyright (c) 2016-2019 David Posada. All rights reserved.

        Name:		CellCoal
        Purpose:	To simulate the clonal neutral evolution of single cells
(for example within a tumor) Code:		David Posada (DP) Started:
01 June 2016

 Version 0.1
- started with my code from SNPsim 0.9
- thinnned coalescent without recombination
- attached an additional node to the coalescent MRCA node. This healthy node
will be the root for the simulation.
- attached a contemporaneous healthy cell as outgroup
- new data structures
- evolve maternal then paternal genomes over coalescent tree
- redefined output files
- branch length healthy tip and root user defined
- reimplemented ISM diploid (standard) and implemented ISM haploid model
- checked mutation distribution (made a fix for snpsim at the simISMsite
function)
- reimplemented fixed SNP model (distributed at random in maternal/paternal
genomes)
- reimplemented JC for 0/1 data (=MK2 or CFN model) (made a fix for snpsim at
the simJCsite function)
- introduced rate variation among sites
- added rate variation among lineages
- added DNA ISM model
- added rate matrix
- produce read counts
- output format CTAG
- added errors to the NGS reads
- calculate genotype likelihoods (checked to work fine with ANGSD)
- added ADO
- generate reference genome (referenceAllele in struct site)
- implemented HKY+G and nested models (JC69, K80, F81)
- printf FASTA haplotypes with IUPAC codes
- print replicates in different files
- specify results folder name
- add germline variation
- automatically clean previous results directory
- different error probabilities for different bases (Eij)
- can add error directly to genotypes (no reads produced in such a case)
- overdispersed coverage heterogeneity following a negative binomial
- implemented two models for amplification errors. The frequency of wrong
templates (max 1; or all possible 3) follows a Beta distribution.
- when ADO we assume there can be a proportion, controlled by a parameter, in
the number of reads
- added GTR and GTnR (avoided old redundant calculations)
- print consensus haps for binary data
- simulate point deletions
- added a specific function for JC
- changed letters/symbols for command line arguments
- prints true variants to a file
- can read and use a user tree
- genotype error and ADO now parameterized at the genotype level
- simulate reads for and print all sites in the VCF (changed)
- implemented trinucleotide mutational signatures (alone or mixture)
- implemented user genome
- implemented copy-neutral LOH
- implemented doublets at the NGS reads level
- added option for allelic imbalance
- NGS data: alternate alleles are defined by the reads
- NGS data: SNVs are defined according to the ML genotypes
- NGS data: genotypes and haplotypes show ML genotypes
- VCF: GL is now G10; GL now follows VCF standards
- variable ADO rates among sites and/or cells
- transforming/healthytip branch lengths are not a function of the tumor MRCA /
total depth
- option for different taxa names for tumor/non-tumor scenarios
- fixed memory leaks and better printing for VCF (thanks to Alexey Kozlov)
- fixed ISM binary models so germline sites can also accumulate somatic
mutations
- new implementation of the 2-template model that corrects the genotype
likelihoods
- can use custom name for parameter file

[TO-DOs]
- add new signatures
- move data structure into cell structure (optional)
- dated tips (optional)
- population structure (optional)
- change random number generator (optional)

[MODELS implemented]
0:ISM_01
1:ISMhap_01
2:M2_01
3:ISM_DNA
4:ISMhap_DNA
5:HKY_JC_DNA
6:GTR_DNA
7:GTRnrDNA
*/

#include "cellcoal.h"
#include "definitions.h"

int main(int argc, char **argv) {
  FILE *fp;
  int i, j, k;
  int dataSetNum;
  long int seed;
  double a, a2;
  double start, secs;

  /* Default settings */
  numThreads = 1;   /* the number of threads */
  numDataSets = 10; /* the number of samples to simulate */
  numCells = 8;     /* number of ingrou cells in each data set */
  ploidy = 2;       /* we assume diploid genomes */
  numSites = 1000;  /* number of sites (markers, loci) to simulate = length of
                       the genome */
  N = 1000;         /* effective population size */
  numPeriods = 0;   /* number of distinct demographic periods */
  doDemographics = NO;   /* whether to implement demographics */
  doExponential = NO;    /* whether to do exponential growth */
  growthRate = 1.0e-5;   /* rate for the exponential population growth */
  mutationRate = 1.0e-7; /* nucleotide mutation rate per site per generation */
  SNPrate = 0.0;         /* germline variation rate for rooth healthy genome */
  rateVarAmongLineages =
      NO;         /* modify rate variation among branches  (to non-clock) */
  alphabet = DNA; /* alphabet 0/1 or DNA */
  altModel = 0;   /* by default the alternative model will be ISM haploid */
  propAltModelSites = 0; /* proportion of sites that will mutate according to
                            alternative model */
  nonISMRelMutRate =
      1.0;             /* relative rate alternative/default model for sites */
  equalBaseFreq = YES; /* DNA base frequencies */
  freq[0] = freq[1] = freq[2] = freq[3] = 0.25;
  titv = 0.5;      /* transition/transversion rate ratio */
  thereIsMij = NO; /* mutation rate matrix*/
  Mij[0][0] = Mij[1][1] = Mij[2][2] = Mij[3][3] =
      0; /* mutation probabilities */
  Mij[0][1] = Mij[0][2] = Mij[0][3] = 1.0;
  Mij[1][0] = Mij[1][2] = Mij[1][3] = 1.0;
  Mij[2][0] = Mij[2][1] = Mij[2][3] = 1.0;
  Mij[3][0] = Mij[3][1] = Mij[3][2] = 1.0;
  thereIsEij = NO; /* error rate matrix*/
  Eij[0][0] = Eij[1][1] = Eij[2][2] = Eij[3][3] =
      0; /* sequencing error probabilities */
  Eij[0][1] = Eij[0][2] = Eij[0][3] = 1.0;
  Eij[1][0] = Eij[1][2] = Eij[1][3] = 1.0;
  Eij[2][0] = Eij[2][1] = Eij[2][3] = 1.0;
  Eij[3][0] = Eij[3][1] = Eij[3][2] = 1.0;
  doJC = YES;
  doHKY = NO;
  doGTR = NO;
  doGTnR = NO;
  rateVarAmongSites =
      NO; /* rate variation among different sites along the genome */
  alphaSites = infinity;    /* alpha shape of the gamma distribution for rate
                               variation among sites */
  alphaBranches = infinity; /* alpha shape of the gamma distribution for rate
                               variation among lineages */
  alphaCoverage =
      infinity; /* alpha shape of the gamma distribution for coverage */
  doSimulateFixedNumMutations =
      NO; /* whether to simulate a fixed number of mutations */
  doUserTree =
      NO; /* whether to assume a user tree instead od making the coalescent */
  doUserGenome =
      NO; /* whether to use a user genome instead of a simulated one */
  doPrintSNVgenotypes = NO;   /* whether to print SNVs */
  doPrintSNVhaplotypes = NO;  /* whether to print haplotypes */
  doPrintTrueHaplotypes = NO; /* whether to print the true haplotypes s */
  doPrintFullHaplotypes = NO; /* whether to print sequences */
  doPrintFullGenotypes =
      NO; /* whether to print all genotypes (variable + invariable) */
  doPrintTree = YES;     /* whether to print the coalescent tree */
  doPrintTimes = NO;     /* whether to print coalescent times */
  doPrintAncestors = NO; /* whether to print data for ancestral cells */
  doNGS = NO;            /* do not produce reads by default */
  doTumorNames = NO;     /* use specific name for taxa in the tumor scenario */
  stringPrecision = 12;  /* precision for taxa names */
  doPrintCATG = NO; /* whether to print read counts for SNVs in CATG format*/
  doSimulateData = YES; /* whether to simulate any data (or just look at the
                           expectations; useful for debugging) */
  doPrintSeparateReplicates =
      NO; /* whether to put every replica in its own file */
  doPrintIUPAChaplotypes = NO; /* whether to print IUPAC halotypes */
  doGeneticSignatures = NO;    /* whether to use a genetic signature to model
                                  trinucleotide mutations */
  numUserSignatures = 0; /* by default we do not use a genetic signature */
  transformingBranchLengthRatio =
      0.5; /* ratio of the transforming branch length compare to the tumor MRCA
              length */
  healthyTipBranchLengthRatio =
      1.0; /* ratio the branch leading to the healthy cell to the (tumor MRCA +
              transfoming branch) length*/
  coverage = 0;         /* NGS  depth for read counts */
  rateVarCoverage = NO; /* there is coverage dispersion */
  ADOrate = 0;          /* allelic dropout */
  ADOvarAmongSites =
      NO; /* ADO variation among different sites along the genome */
  ADOvarAmongCells = NO; /* ADO variation among different cells */
  alphaADOsites =
      infinity; /* alpha shape allelic dropout heterogeneity among sites */
  alphaADOcells =
      infinity; /* alpha shape allelic dropout heterogeneity among cells */
  sequencingError = 0;        /* NGS error rate */
  genotypingError = 0;        /* add errors directly in the genotypes */
  meanAmplificationError = 0; /* mean of beta distribution for WGA errors */
  varAmplificationError =
      0.000001; /* variance of beta distribution for WGA errors */
  simulateOnlyTwoTemplates =
      NO; /* whether simulating maximum of two templates after single-cell
             amplification, or there can be all four */
  haploidCoverageReduction =
      0.5; /* proportion of reads produced when a single allele is present */
  allelicImbalance = 0.5; /* proportion of maternal/ paternal reads */
  meanDoubletRate = 0.0;      /* no doublets by default */
  varDoubletRate = 1.0e-06;
  numNodes = 1000; /* initial number of nodes allocated to build the coalescent
                      trees */
  seed = time(NULL); /* seed for random numbers */
  doSpecificParameterFile =
      NO;       /* by default do not use a custom parameter file */
  userSeed = 0; /* seed entered by the user */
  noisy = 1; /* level of information to be printed in the screen (see below) */

  /*
   noisy = 0: does not print anything
   noisy = 1: + simulation summary,
   noisy = 2: + replicate information,
   noisy = 3: + calculation status and event information
  */

  strcpy(parameterFile, "parameters");
  ReadParametersFromCommandLine(argc, argv);
  if (argc < 2 || doSpecificParameterFile == YES) {
    readingParameterFile = YES;
    if ((fp = freopen(parameterFile, "r", stdin)) != NULL)
      ReadParametersFromFile();
    else {
      fprintf(stderr, "\nERROR: Could not open parameter file: %s\n",
              parameterFile);
      PrintUsage(stderr);
      exit(-1);
    }
  } else
    readingParameterFile = NO;

  start = omp_get_wtime();

  if (noisy > 0) {
    PrintHeader(stderr);
    PrintDate(stderr);
  }

  if (doSimulateData == NO) {
    doPrintSNVgenotypes = NO;
    doPrintSNVhaplotypes = NO;
    doPrintTrueHaplotypes = NO;
    doPrintFullHaplotypes = NO;
    doPrintFullGenotypes = NO;
    doPrintAncestors = NO;
    doNGS = NO;
    doPrintCATG = NO;
  }

  /* initialize a few variables */
  numNodes = 2 * numCells + 1;
  AllocateTreeStructure(numNodes);
  a = 0;
  a2 = 0;
  for (i = 1; i < numCells; i++) /* ignoring outgroup (healthy cell) */
  {
    a += 1.0 / (double)i;
    a2 += 1.0 / (double)pow(i, 2);
  }
  theta = 4.0 * N * mutationRate * 2.0 * numSites *
          ((1.0 - propAltModelSites) + nonISMRelMutRate * propAltModelSites);
  expNumMU = theta * a;
  expVarNumMU = theta * a + pow(theta, 2) * a2;
  expTMRCA = 2 * (1 - 1.0 / numCells) /* 2.0 * N */;
  expVarTMRCA = 0;
  for (i = 2; i <= numCells; i++)
    expVarTMRCA += /* 4.0 * pow(N,2) */ 4.0 / (pow(i, 2) * pow(i - 1, 2));
  cumNumCA = cumNumMU = cumNumMUSq = cumNumSNVs = cumNumSNVsSq =
      cumCountMLgenotypeErrors = 0;
  zeroSNVs = cumTMRCA = cumTMRCASq = 0;

  /* added code */
  cumNumDEL = cumNumINS = 0;
  /* end added code */

  altModelMutationRate = mutationRate * nonISMRelMutRate;

  if (doSimulateData == YES && alphabet == DNA) {
    /* initialize cumfreq */
    cumfreq[0] = freq[0];
    for (i = 1; i < 4; i++)
      cumfreq[i] = cumfreq[i - 1] + freq[i];

    /* initialize cumMij */
    for (i = 0; i < 4; i++) {
      cumMij[i][0] = Mij[i][0];
      for (j = 1; j < 4; j++)
        cumMij[i][j] = cumMij[i][j - 1] + Mij[i][j];
    }

    /* initialize cumEij */
    for (i = 0; i < 4; i++) {
      cumEij[i][0] = Eij[i][0];
      for (j = 1; j < 4; j++)
        cumEij[i][j] = cumEij[i][j - 1] + Eij[i][j];
    }
  }

  /* Set seed and spin wheels of pseudorandom number generator */
  /* seed = (unsigned int) clock();*/
  if (userSeed > 0)
    seed = userSeed;

  InitRandomDevices(seed);

  if (noisy > 0) {
    PrintCommandLine(stderr, argc, argv);
  }

  /* Check the distribution of the mutations */
#ifdef CHECK_MUT_DISTRIBUTION
  MutCount = (int *)calloc(2 * (numCells + 1), sizeof(int));
  dataSetsWithSNVs = 0;
  meansumPos = 0;
#endif

  /* set names for taxa */
  if (doTumorNames == YES) {
    strcpy(inCellName, "tumcell");
    strcpy(outCellName, "healthycell");
    strcpy(inRootCellName, "tumoralroot");
    strcpy(outRootCellName, "healthyroot");
    stringPrecision = 12;
  } else {
    strcpy(inCellName, "cell");
    strcpy(outCellName, "outgcell");
    strcpy(inRootCellName, "ingrroot");
    strcpy(outRootCellName, "outgroot");
    stringPrecision = 9;
  }

  if (doPrintIUPAChaplotypes == YES)
    stringPrecision += 1;

  /* set file dirs and names */
  if (strlen(resultsDir) == 0)
    strcpy(resultsDir, "Results");
  // clean previous results directory
  RemoveDir(resultsDir);
  // create new results directory
  mkdir(resultsDir, S_IRWXU);

  strcpy(treeDir, "trees_dir");
  strcpy(timesDir, "times_dir");
  strcpy(SNVgenotypesDir, "snv_genotypes_dir");
  strcpy(SNVhaplotypesDir, "snv_haplotypes_dir");
  strcpy(trueHaplotypesDir, "true_haplotypes_dir");
  strcpy(fullHaplotypesDir, "full_haplotypes_dir");
  strcpy(MLhaplotypesDir, "ML_haplotypes_dir");
  strcpy(fullGenotypesDir, "full_genotypes_dir");
  strcpy(CATGdir, "catg_dir");
  strcpy(VCFdir, "vcf_dir");
  /* added code */
  strcpy(treeWithTrunkDir, "with_trunk_without_normal_dir");
  strcpy(treeWithoutTrunkDir, "without_trunk_without_normal_dir");
  strcpy(treeWithHealthyTipDir, "without_trunk_with_normal_dir");
  strcpy(cellNamesDir, "cell_names_dir");
  strcpy(sizeFactorsDir, "size_factors_dir");
  strcpy(SNVSitesNamesDir, "true_evo_events_sites_names_dir");
  strcpy(SNVAllelicSeqInfoDir, "true_evo_events_allelic_seq_info_dir");
  strcpy(trueSNVDir, "true_evo_events_dir");
  strcpy(parallelMutDir, "parallel_mutation_dir");
  strcpy(rawSeqDataDir, "raw_seq_data_dir");
  strcpy(dummyBamDataDir, "dummy_bam_data_dir");
  strcpy(mpileupDataDir, "mpileup_data_dir");
  strcpy(adoStatesDir, "ado_state_dir");
  /* end added code */

  strcpy(SNVgenotypesFile, "snv_gen");
  strcpy(SNVhaplotypesFile, "snv_hap");
  strcpy(trueHaplotypesFile, "true_hap");
  strcpy(fullHaplotypesFile, "full_hap");
  strcpy(MLhaplotypesFile, "ML_hap");
  strcpy(fullGenotypesFile, "full_gen");
  // strcpy(treeFile, "trees");
  strcpy(timesFile, "times");
  if (strlen(userTreeFile) == 0)
    strcpy(userTreeFile, "usertree");
  if (strlen(userGenomeFile) == 0)
    strcpy(userGenomeFile, "usergenome");
  strcpy(CATGfile, "catg");
  strcpy(VCFfile, "vcf");
  strcpy(logFile, "log");
  strcpy(settingsFile, "log");
  /* added code */
  strcpy(sizeFactorsFile, "size_factors");
  strcpy(cellNamesFile, "cell_names");
  strcpy(parallelMutFile, "parallel_mutations");
  strcpy(SNVSitesNamesFile, "true_evo_events_sites_names");
  strcpy(SNVAllelicSeqInfoFile, "true_evo_events_allelic_seq_info");
  if (altModel == FiniteMu || altModel == FiniteMuDel ||
      altModel == FiniteMuInDel)
    strcpy(treeFile, "tree");
  else
    strcpy(treeFile, "trees");
  strcpy(treeWithoutTrunkFile, "tree");
  strcpy(treeWithHealthyTipFile, "tree");
  strcpy(trueSNVGenotypesAlFile, "true_evo_events_genotypes_alphabetic");
  strcpy(trueSNVGenotypesNrFile, "true_evo_events_genotypes_numeric");
  strcpy(trueSNVTernaryMatrixFile, "true_evo_events_ternary_matrix");
  strcpy(rawSeqDataRefFile, "reference");
  strcpy(trueSNVCovFile, "true_evo_events_coverage");
  strcpy(evolutionaryEventsFile, "historical_evo_events_summary");
  strcpy(mpileupDataNMFile, "normal_merged");
  strcpy(mpileupDataMFile, "merged");
  strcpy(mpileupDataNUFile, "normal_unmerged");
  strcpy(mpileupDataUFile, "unmerged");
  strcpy(adoStatesFile, "ado_states");
/* end added code */
#ifdef MYDEBUG
  strcpy(mutationsFile, "mutations");
#endif

  /* added code */
  if (doPrintSeparateReplicates == NO) {
    if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu)) {
      PrepareGlobalFiles(argc, argv);
    } else {
      PrepareGlobalFilesExtd(argc, argv);
    }
  }
  /* end added code */

  /* original code */
  /*
  if (doPrintSeparateReplicates == NO)
    PrepareGlobalFiles(argc, argv);
    */

  /* Start printing info to log */
  sprintf(File, "%s/%s", resultsDir, settingsFile);
  if ((fpLog = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  PrintHeader(fpLog);
  PrintDate(fpLog);
  sprintf(File, "%s/%s", resultsDir, settingsFile);
  PrintCommandLine(fpLog, argc, argv);

  /* Read user treefile */
  if (doUserTree == YES) {
    if (noisy > 2)
      fprintf(stderr, "\n>> Reading user tree ...");
    ReadUserTree(fpUserTree);
  }

  /* Read user genomefile  */
  if (doUserGenome == YES) {
    if (noisy > 2)
      fprintf(stderr, "\n>> Reading user genome ...");
    ReadUserGenome(fpUserGenome);
  }

  OUTGROUP_ROOT = 2 * numCells;
  INGROUP_ROOT = (2 * numCells) - 1;

  if (doGeneticSignatures == YES && doSimulateData == YES)
    PrepareGeneticSignatures();

  /* added code */
  if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
      altModel == FiniteMu) {
    InitThreadsStructure(numCells, numSites, numThreads);
    AllocateRawDataPrinters(numThreads);
  }

  allelicCov = (double *)malloc((numSites + 1) * sizeof(double));
  allelicRawVar = (double *)malloc((numSites + 1) * sizeof(double));

  sizeFactors = (double *)malloc((numCells + 2) * sizeof(double));
  /* end added code */

  for (dataSetNum = 0; dataSetNum < numDataSets; dataSetNum++) {
    numCA = numMU = numDEL = numCNLOH = numProposedMU = TMRCA = numSNVmaternal =
        0;

    /* added code */
    numINS = 0;

    if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
        altModel == FiniteMu) {
      if (dataSetNum == 0)
        rawDataFiles =
            (RawDataFile *)malloc((2 * numCells + 3) * sizeof(RawDataFile));

      completeDeletions = completeADOs = noCoverage = noReadsPositiveCoverage =
          zeroCoverage = zeroProportionVAF = 0;

      for (i = 0; i < 4; ++i) {
        for (j = 0; j < 2; ++j) {
          allelicZeroCoverageProportion[i][j] = 0;
        }
      }

      if (fixedTree == YES && numDataSets > 1 && dataSetNum == 0) {
        completeDeletionsAll = completeADOsAll = noCoverageAll =
            noReadsPositiveCoverageAll = zeroCoverageAll =
                zeroProportionVAFAll = 0;

        for (i = 0; i < 4; ++i) {
          for (j = 0; j < 2; ++j) {
            allelicZeroCoverageProportionAll[i][j] = 0;
          }
        }
      }
    }

    if (doPrintSeparateReplicates == YES) {
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        PrepareSeparateFiles(dataSetNum);
      } else {
        PrepareSeparateFilesExtd(dataSetNum);
      }
    } else {
      if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu) {
        PrepareSeparateFoldersAndFiles(dataSetNum);
      }
    }
    /* end added code */

    /* original code */
    /*
    if (doPrintSeparateReplicates == YES)
      PrepareSeparateFiles(dataSetNum);
      */

    if (noisy == 1) {
      fprintf(stdout, "\rReplicate #%3d/%d ", dataSetNum + 1, numDataSets);
      fflush(stdout);

      /* if stdout is redirected elsewhere, print also to stderr */
      if (!isatty(fileno(stdout))) {
        fprintf(stderr, "\rreplicate #%3d/%d", dataSetNum + 1, numDataSets);
        fflush(stderr);
      }
    }

    /* Generate coalescence tree */
    if (doUserTree == NO) {
      if (noisy > 2) {
        fprintf(
            stderr,
            "\n______________________________________________________________"
            "________\n");
        fprintf(stderr, "\n\rReplicate #%3d/%d ", dataSetNum + 1, numDataSets);
        fprintf(stderr, "\n\n>> Starting coalescent tree ...");
      }

      /* added code */
      if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu) {
        doCoalescenceTree =
            DoCoalescenceTree(numDataSets, dataSetNum, fixedTree, noisy);
      }
      /* end added code */

      if (doCoalescenceTree == YES) {
        MakeCoalescenceTree(numCells, N, mainThreadRandDev);
      }
      if (noisy > 2)
        fprintf(stderr, "\n>> Finishing coalescent tree ... DONE\n");

      /* Make cell tree non-clock if needed */
      if (rateVarAmongLineages == YES)
        MakeTreeNonClock(coalTreeMRCA, mainThreadRandDev);

      if (doCoalescenceTree == YES) {
        cumNumCA += numCA;
        cumTMRCA += TMRCA;
        cumTMRCASq += pow(TMRCA, 2);
      }
    }

    /* Put the tree structure in a matrix */
    InitTreeStructure(numNodes);
    //    PrintTreeStructure(numNodes);
    SetTreeStructure(healthyRoot, numNodes);
    //    PrintTreeStructure(numNodes);

    if (doCoalescenceTree == YES) {
      if (dataSetNum == 0)
        AllocateCellNamesExtd(numCells);

      if (doPrintTree == YES)
        PrintTree(healthyRoot, fpTrees);

      PrintTree(healthyRoot->left, fpTreeWithoutTrunk);

      {
        int tmp = doPrintHealthyTip;
        doPrintHealthyTip = YES;
        PrintTree(healthyRoot, fpTreeWithHealthyTip);
        doPrintHealthyTip = tmp;
      }

      if (dataSetNum == 0 && !(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu)) {
        PrintTumCellNames(fpCellNames, ',', 0);
      }
    }

    if (doPrintTimes == YES && doCoalescenceTree == YES)
      PrintTimes(0, fpTimes);

    totalTreeLength = SumBranches(healthyRoot);

    /* added code */
    /* allocate genotype data to be stored in data[genome][cell][site] */
    if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu)) {
      data = (int ***)calloc(ploidy, sizeof(int **));
      if (!data) {
        fprintf(stderr, "Could not allocate the data structure\n");
        exit(-1);
      }
      for (i = 0; i < ploidy; i++) {
        data[i] = (int **)calloc(2 * numCells + 1, sizeof(int *));
        if (!data[i]) {
          fprintf(stderr, "Could not allocate the data[] structure\n");
          exit(-1);
        }
        for (j = 0; j < 2 * numCells + 1; j++) {
          data[i][j] = (int *)calloc(numSites, sizeof(int));
          if (!data[i][j]) {
            fprintf(stderr, "Could not allocate the data[][] structure\n");
            exit(-1);
          }
        }
      }
    }
    /* end added code */

    /* allocate genotype data to be stored in data[genome][cell][site] */
    /* original code
    data = (int ***)calloc(ploidy, sizeof(int **));
    if (!data) {
      fprintf(stderr, "Could not allocate the data structure\n");
      exit(-1);
    }
    for (i = 0; i < ploidy; i++) {
      data[i] = (int **)calloc(2 * numCells + 1, sizeof(int *));
      if (!data[i]) {
        fprintf(stderr, "Could not allocate the data[] structure\n");
        exit(-1);
      }
      for (j = 0; j < 2 * numCells + 1; j++) {
        data[i][j] = (int *)calloc(numSites, sizeof(int));
        if (!data[i][j]) {
          fprintf(stderr, "Could not allocate the data[][] structure\n");
          exit(-1);
        }
      }
    }
    */

    if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu)) {
      AllocateCellStructure();
    } else {
      AllocateCellStructureExtd(numCells, numSites);
    }
    /* allocate cell structure */
    // AllocateCellStructure();

    /* allocate memory for site information (equal for maternal and paternal) */
    allSites = (SiteStr *)calloc(numSites + 1, sizeof(SiteStr));
    if (!allSites) {
      fprintf(stderr, "Could not allocate the allSites structure\n");
      exit(-1);
    }
    for (i = 0; i < numSites; i++) {
      allSites[i].alternateAlleles = (int *)calloc(4, sizeof(int));
      if (!allSites[i].alternateAlleles) {
        fprintf(stderr, "Could not allocate the allSites[i].alternateAlleles[] "
                        "structure\n");
        exit(-1);
      }

      /* added code */
      allSites[i].hasEventsOnLeaves = NO;
      /* end added code */
    }

    /* the arrays below keep the index for different types of sites */
    SNVsites = (int *)calloc(numSites, sizeof(int));
    if (!SNVsites) {
      fprintf(stderr, "Could not allocate the SNVsites structure\n");
      exit(-1);
    }

    /* the arrays below keep the index for different types of sites */
    variantSites = (int *)calloc(numSites, sizeof(int));
    if (!variantSites) {
      fprintf(stderr, "Could not allocate the variantSites structure\n");
      exit(-1);
    }

    /* added code */
    if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
        altModel == FiniteMu) {
      mutatedSites = (int *)calloc(numSites, sizeof(int));
      if (!mutatedSites) {
        fprintf(stderr, "Could not allocate the mutatedSites structure\n");
        exit(-1);
      }

      insertedSites = (int *)calloc(numSites, sizeof(int));
      if (!insertedSites) {
        fprintf(stderr, "Could not allocate the insertedSites structure\n");
        exit(-1);
      }

      deletedSites = (int *)calloc(numSites, sizeof(int));
      if (!deletedSites) {
        fprintf(stderr, "Could not allocate the deletedSites structure\n");
        exit(-1);
      }
    }
    /* end added code */

    DefaultModelSites = (int *)calloc(numSites, sizeof(int));
    if (!DefaultModelSites) {
      fprintf(stderr, "Could not allocate the DefaultModelSites structure\n");
      exit(-1);
    }

    AltModelSites = (int *)calloc(numSites, sizeof(int));
    if (!AltModelSites) {
      fprintf(stderr, "Could not allocate the AltModelSites structure\n");
      exit(-1);
    }

    if (doSimulateData == YES) {
      numISMmutations = 0;
      numISMdeletions = 0;
      numISMCNLOH = 0;

      if (noisy > 2)
        fprintf(stderr, "\n>> Starting evolving genome along the tree ... ");

      /* assign ancestral states to all genomes */
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        InitializeGenomes(healthyRoot, mainThreadRandDev);
      } else {
        if (noisy > 2)
          fprintf(stderr, "\n>> Initializing the reference genome ... ");
        InitializeGenomesExtd(healthyRoot, mainThreadRandDev);
        if (noisy > 2)
          fprintf(stderr, "DONE");
      }

      /* add germline variation if needed */
      if (SNPrate > 0)
        AddGermlineVariation(mainThreadRandDev);

      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        /* evolve maternal genome */
        if (noisy > 2)
          fprintf(stderr, "\n>> Evolving maternal genome ... ");
        EvolveSitesOnTree(healthyRoot, MATERNAL, mainThreadRandDev);
        if (noisy > 2)
          fprintf(stderr, "DONE");

        /* evolve paternal genome */
        if (noisy > 2)
          fprintf(stderr, "\n>> Evolving paternal genome ... ");
        EvolveSitesOnTree(healthyRoot, PATERNAL, mainThreadRandDev);
        if (noisy > 2)
          fprintf(stderr, "DONE\n\n");
      } else {
        /* pick finiteMuInDel model */
        if (noisy > 2)
          fprintf(stderr, "\n>> Evolving genome along the tree ... \n");
        EvolveSitesOnTreeExtd(healthyRoot, mainThreadRandDev);
        if (noisy > 2) {
          fprintf(stderr, "DONE");
          if (numSites >= 1000) {
            fprintf(stderr, " @ ");
            PrintDate(stderr);
          } else {
            fprintf(stderr, "\n");
          }
        }
      }

      cumNumMU += numMU;
      cumNumMUSq += pow(numMU, 2);

      /* added code */
      cumNumDEL += numDEL;
      cumNumDELSq += pow(numMU, 2);

      cumNumINS += numINS;
      cumNumINSSq += pow(numMU, 2);
      /* end added code */

      /* added code */
      if (CNLOHrate > 0 && !(altModel == FiniteMuInDel ||
                             altModel == FiniteMuDel || altModel == FiniteMu)) {
        /* end added code */
        // if (CNLOHrate > 0) {
        /* evolve maternal CN_LOH */
        if (noisy > 2)
          fprintf(stderr, "\n>> Evolving maternal CN_LOH ... ");
        EvolveCNLOHonTree(healthyRoot, MATERNAL, mainThreadRandDev);
        if (noisy > 2)
          fprintf(stderr, "DONE");

        /* evolve paternal CN_LOH  */
        if (noisy > 2)
          fprintf(stderr, "\n>> Evolving paternal CN_LOH ... ");
        EvolveCNLOHonTree(healthyRoot, PATERNAL, mainThreadRandDev);
        if (noisy > 2)
          fprintf(stderr, "DONE\n\n");

        cumNumCNLOH += numCNLOH;
        cumNumCNLOHSq += pow(numCNLOH, 2);
      }

      /* added code */
      if (deletionRate > 0 &&
          !(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        /* end added code */
        // if (deletionRate > 0) {
        /* evolve maternal deletions */
        if (noisy > 2)
          fprintf(stderr, "\n>> Evolving maternal deletions ... ");
        EvolveDeletionsOnTree(healthyRoot, MATERNAL, mainThreadRandDev);
        if (noisy > 2)
          fprintf(stderr, "DONE");

        /* evolve paternal deletions  */
        if (noisy > 2)
          fprintf(stderr, "\n>> Evolving paternal deletions ... ");
        EvolveDeletionsOnTree(healthyRoot, PATERNAL, mainThreadRandDev);
        if (noisy > 2)
          fprintf(stderr, "DONE\n\n");

        cumNumDEL += numDEL;
        cumNumDELSq += pow(numDEL, 2);
      }

      /* added code */
      if (noisy > 2)
        fprintf(stderr, "\n>> Finishing evolving the genome ... DONE\n");
      /* end added code */

      /* count how many SNVs we observe before ADO and genotype errors, but
       * after CNLOH and deletion */
      // numSNVs = CountTrueVariants();

      /* added code */
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        for (i = 0; i < numCells + 1; i++)
          for (j = 0; j < numSites; j++) {
            cell[i].site[j].trueMaternalAllele = data[MATERNAL][i][j];
            cell[i].site[j].truePaternalAllele = data[PATERNAL][i][j];
          }

        numSNVs = CountTrueVariants();
      } else {
        for (j = 0; j < numSites; j++) {
          int threadIndex = GetThreadIndex(j);
          int relaSite = j - threadStr[threadIndex].start;
          for (i = 0; i < numCells + 1; i++) {
            cellExtd[i].site[j].trueCode = cellExtd[i].site[j].observedCode =
                threadStr[threadIndex].data[i][relaSite].code;
            cellExtd[i].site[j].trueNumExistingAlleles =
                cellExtd[i].site[j].observedNumExistingAlleles =
                    threadStr[threadIndex].data[i][relaSite].numExistingAlleles;
            cellExtd[i].site[j].trueAlleles =
                threadStr[threadIndex].data[i][relaSite].alleles;
            cellExtd[i].site[j].pastEvolutionaryEvents =
                threadStr[threadIndex].data[i][relaSite].pastEvolutionaryEvents;

            if (cellExtd[i].site[j].trueNumExistingAlleles == 0) {
              cellExtd[i].site[j].maxNumExistingAllelesPerStrand = 0;
            } else if (cellExtd[i].site[j].trueNumExistingAlleles > 2) {
              cellExtd[i].site[j].maxNumExistingAllelesPerStrand = 2;
            } else {
              cellExtd[i].site[j].maxNumExistingAllelesPerStrand = 1;
            }

            int *pCell = cellExtd[i].site[j].observedAlleles;
            int *pData = threadStr[threadIndex].data[i][relaSite].alleles;
            for (k = 0; k < 3; ++k) {
              *(pCell++) = *(pData++);
            }

            threadStr[threadIndex].data[i][relaSite].alleles =
                cellExtd[i].site[j].observedAlleles;
          }
        }

        numSNVs = CountTrueEvents();
      }
      /* end added code */

      /* copy true genotypes before doing ADO but after CNLOH and deletion */
      /* original code
      for (i = 0; i < numCells + 1; i++)
        for (j = 0; j < numSites; j++) {
          cell[i].site[j].trueMaternalAllele = data[MATERNAL][i][j];
          cell[i].site[j].truePaternalAllele = data[PATERNAL][i][j];
        }
        */

      /* print reference haplotypes without errors */
      if (doPrintTrueHaplotypes == YES) {
        if (doPrintSeparateReplicates == NO)
          fprintf(fpTrueHaplotypes, "[#%d]\n", dataSetNum + 1);
        PrintTrueFullHaplotypes(fpTrueHaplotypes);
      }

      /* do alellic dropout */
      if (ADOrate > 0) {
        /* added code */
        if (noisy > 2)
          fprintf(stderr, "\n>> Performing allelic dropout (%s) ... ",
                  (ADOtype == SINGLE_ADO) ? "single ADO" : "locus dropout");
        /* end added code */
        AllelicDropout(mainThreadRandDev);
        /* added code */
        if (noisy > 2)
          fprintf(stderr, "DONE\n");
        /* end added code */
      }

      /* introduce errors directly in the genotypes */
      if (genotypingError > 0)
        GenotypeError(mainThreadRandDev);

      /* generate NGS read counts */
      if (doNGS == YES && !(altModel == FiniteMuInDel ||
                            altModel == FiniteMuDel || altModel == FiniteMu))
        GenerateReadCounts(mainThreadRandDev);

      if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu) {
        if (noisy > 2) {
          fprintf(stderr, "\n>> Generating read counts ... ");
          if (numSites >= 1000) {
            fprintf(stderr, " @ ");
            PrintDate(stderr);
          }
        }
        GenerateReadCountsExtd(dataSetNum, mainThreadRandDev);
        if (noisy > 2) {
          fprintf(stderr, "DONE");
          if (numSites >= 1000) {
            fprintf(stderr, " @ ");
            PrintDate(stderr);
          } else {
            fprintf(stderr, "\n");
          }
        }
      }

      /* added code */
      /* count how many alleles are at each site and how many SNVs we observe
       * now, after ADO and genotype errors/read generation */
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        if (doNGS == YES)
          numSNVs = CountAllelesInMLGenotypes();
        else
          numSNVs = CountAllelesInObservedGenotypes();
      }

      if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu) {
        PrintTumCellNames(fpCellNames, ',', 0);
        PrintSNVSitesNames(fpSNVSitesNames);
        PrintSNVAllelicInformation(fpSNVAllelicSeqInfo);
        PrintAdoStates(fpAdoStates);
        PrintSizeFactors(fpSizeFactors);
        PrintParallelMut(fpParallelMut);

        for (i = 0; i < numCells + 1; ++i) {
          for (j = 0; j < numSNVs; ++j) {
            /* reasons for zero VAF:
             * zero coverage
             * zero reads for variant nucleotide with positive coverage */
            allelicZeroCoverageProportion
                [cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles][0]++;

            ProcessZeroCoverageAndZeroProportion(
                i, SNVsites[j], allSites[SNVsites[j]].MLVariantNucleotide,
                numDataSets);
          }
        }

        ProcessCoverageProb(dataSetNum);
      }
      /* end added code */

      /* original code */
      /*
      if (doNGS == YES)
        numSNVs = CountAllelesInMLGenotypes();
      else
        numSNVs = CountAllelesInObservedGenotypes();
      */

      cumNumSNVs += numSNVs;
      cumNumSNVsSq += pow(numSNVs, 2);

      if (numSNVs == 0)
        zeroSNVs++;

      /* added code */
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        if (noisy > 1) {
          fprintf(stderr, "\nData set %d", dataSetNum + 1);
          fprintf(stderr, "\nNumber of coalescence events   =   %d", numCA);
          fprintf(stderr, "\nNumber of mutational events    =   %d", numMU);
          fprintf(stderr, "\nNumber of SNVs                 =   %d", numSNVs);
          fprintf(stderr, "\nNumber of deletion events      =   %d", numDEL);
          fprintf(stderr, "\nNumber of cnLOH events         =   %d\n\n",
                  numCNLOH);
        }
      } else {
        PrintSummaryToLog(fpLog, dataSetNum, noisy > 1);
      }
      /* end added code */

      /* original code */
      /*
      if (noisy > 1) {
        fprintf(stderr, "\nData set %d", dataSetNum + 1);
        fprintf(stderr, "\nNumber of coalescence events   =   %d", numCA);
        fprintf(stderr, "\nNumber of mutational events    =   %d", numMU);
        fprintf(stderr, "\nNumber of SNVs                 =   %d", numSNVs);
        fprintf(stderr, "\nNumber of deletion events      =   %d", numDEL);
        fprintf(stderr, "\nNumber of cnLOH events         =   %d\n\n",
                numCNLOH);
      }
       */

      if (doPrintSNVgenotypes == YES &&
          numSNVs > 0) /* we only print replicates with variation */
      {
        /* added code */
        if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
              altModel == FiniteMu)) {
          PrintSNVGenotypes(fpSNVgenotypes);
        } else {
          PrintSNVGenotypesExtd(fpSNVgenotypes, noisy > 3);
        }
        /* end added code */

        /* original code */
        // PrintSNVGenotypes(fpSNVgenotypes);
      }
      if (doPrintSNVhaplotypes == YES &&
          numSNVs > 0) /* we only print replicates with variation */
      {
        if (doPrintSeparateReplicates == NO)
          fprintf(fpSNVhaplotypes, "[#%d]\n", dataSetNum + 1);
        PrintSNVHaplotypes(fpSNVhaplotypes, NO);
      }
      if (doPrintFullGenotypes == YES) {
        if (doPrintSeparateReplicates == NO)
          fprintf(fpFullGenotypes, "[#%d]\n", dataSetNum + 1);

        /* added code */
        if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
              altModel == FiniteMu))
          PrintFullGenotypes(fpFullGenotypes);
        else
          PrintFullGenotypesExtd(fpFullGenotypes, noisy > 3);
        /* end added code */

        /* original code */
        // PrintFullGenotypes(fpFullGenotypes);
      }
      if (doPrintFullHaplotypes == YES) {
        if (doPrintSeparateReplicates == NO)
          fprintf(fpFullHaplotypes, "[#%d]\n", dataSetNum + 1);
        PrintFullHaplotypes(fpFullHaplotypes);
      }

      /* added code */
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        if (alphabet == DNA) {
          if (doNGS == YES) {
            if (doPrintCATG == YES) {
              PrintCATG(fpCATG);
              if (doPrintSeparateReplicates == NO)
                fprintf(fpCATG, "[#%d]\n", dataSetNum + 1);
            }
          }
        }
      } else {
        PrintTrueSNV(fpTrueSNVGenotypesAl, fpTrueSNVGenotypesNr,
                     fpTrueSNVTernaryMatrix);

        PrintEvolutionaryEventsSummary(fpEvolutionaryEvents);

        PrintTrueSNVCov(fpTrueSNVCov);

        const int tumCellRawFileNum = GetTumCellsTotalNum(numCells);
        PrepareRawReadSeqDataFiles(dataSetNum, tumCellRawFileNum);

        InitRawDataPrinters(tumCellRawFileNum, numThreads);

        PrintRawData(rawDataFiles);
        PrintMPileUp(fpMpileupDataNM,
          meanDoubletRate > 0 ? fpMpileupDataNU : NULL, fpMpileupDataM,
          meanDoubletRate > 0 ? fpMpileupDataU : NULL);
      }
      /* end added code */

      /* original code */
      /*
      if (alphabet == DNA) {
        if (doNGS == YES) {
          if (doPrintSeparateReplicates == NO)
            fprintf(fpVCF, "[#%d]\n", dataSetNum + 1);
          PrintVCF(fpVCF);

          if (doPrintCATG == YES) {
            PrintCATG(fpCATG);
            if (doPrintSeparateReplicates == NO)
              fprintf(fpCATG, "[#%d]\n", dataSetNum + 1);
          }
        }
      }
       */
    } // simdata

#ifdef CHECK_MUT_DISTRIBUTION
    /* Check expectations regarding the position of the mutations in the tree */
    /* Find out how many mutations are per site; we expect theta/i number of
       sites with i mutations This approximation might not be good for the
       haploid ISM model We will count how many sites have 1, 2, 3 ...
       numCells-1 mutant alleles */
    SiteMut = (int *)calloc(numSites, sizeof(int));
    for (j = 0; j < numSites; j++)
      for (i = 0; i < numCells; i++) {
        if (data[0][i][j] == 1)
          SiteMut[j]++;
        if (data[1][i][j] == 1)
          SiteMut[j]++;
      }
    for (j = 0; j < numSites; j++)
      for (i = 0; i < numCells; i++)
        if (SiteMut[j] == i)
          MutCount[i]++;

    /* Check expectations regarding the distribution of the mutations along the
     * genome */
    /* the average sum of the positions of the segregating sites has to be L/2
     */
    sumPos = 0;
    if (numSNVs > 0) {
      dataSetsWithSNVs++;
      for (j = 0; j < numSites; j++)
        if (allSites[j].numMutations > 0)
          sumPos += (j + 1);
      meansumPos += (sumPos / (double)numSNVs);
    }
#endif

    /* added code */
    if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu)) {
      if (noisy > 2) {
        for (i = 0; i < numSites; i++) {
          PrintSiteInfo(stderr, i);
        }
      }
    } else {
      if (noisy > 3) {
        for (i = 0; i < numSites; i++) {
          PrintSiteInfoExtd(stderr, i);
        }
      }
    }
    /* end added code */

    /* added code */
    /*
    if (noisy > 2)
      for (i = 0; i < numSites; i++)
        PrintSiteInfo(stderr, i);
        */
    /* end added code */

    // free all the necessary stuff for this replicate
    if (doUserTree == NO) {
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu) ||
          ((altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu) &&
           fixedTree == NO)) {
        for (i = 0; i < numNodes; i++) {
          free(treeNodes[i].name);
          treeNodes[i].name = NULL;
        }
        free(treeNodes);
        treeNodes = NULL;
      }
    }

    /* added code */
    if (doSimulateData == YES &&
        !(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu))
    /* end added code */
    /* original code */
    // if (doSimulateData == YES)
    {
      /* close file pointers for separate replicate file */
      if (doPrintSeparateReplicates == YES) {
        if (doPrintSNVgenotypes == YES)
          fclose(fpSNVgenotypes);
        if (doPrintFullGenotypes == YES)
          fclose(fpFullGenotypes);
        if (doPrintSNVhaplotypes == YES)
          fclose(fpSNVhaplotypes);
        if (doPrintFullHaplotypes == YES)
          fclose(fpFullHaplotypes);
        if (doPrintTrueHaplotypes == YES)
          fclose(fpTrueHaplotypes);
        if (doNGS == YES)
          fclose(fpVCF);
        if (doPrintCATG == YES)
          fclose(fpCATG);
      }

      /* free data structures */
      for (i = 0; i < ploidy; i++) {
        for (j = 0; j < 2 * numCells + 1; j++)
          free(data[i][j]);
        data[i][j] = NULL;
        free(data[i]);
        data[i] = NULL;
      }
      free(data);
      data = NULL;

      for (i = 0; i < numSites; i++) {
        free(allSites[i].alternateAlleles);
        allSites[i].alternateAlleles = NULL;
      }
      free(allSites);
      allSites = NULL;
      free(SNVsites);
      SNVsites = NULL;
      free(variantSites);
      variantSites = NULL;
      free(DefaultModelSites);
      DefaultModelSites = NULL;
      free(AltModelSites);
      AltModelSites = NULL;

      if (doGeneticSignatures == YES) {
        for (i = 0; i < 64; i++) {
          free(triNucleotideMaternal[i].position);
          triNucleotideMaternal[i].position = NULL;
          free(triNucleotidePaternal[i].position);
          triNucleotidePaternal[i].position = NULL;
        }
        free(triNucleotideMaternal);
        triNucleotideMaternal = NULL;
        free(triNucleotidePaternal);
        triNucleotidePaternal = NULL;
      }
    }

    /* original code */
    // if (doNGS == YES) {
    /* added code */
    if (doNGS == YES && !(altModel == FiniteMuInDel ||
                          altModel == FiniteMuDel || altModel == FiniteMu))
    /* end added code */
    {
      for (i = 0; i < numCells + 1; i++) {
        for (j = 0; j < numSites; j++) {
          for (k = 0; k < 4; k++) {
            free(cell[i].site[j].genLike[k]);
            cell[i].site[j].genLike[k] = NULL;
            free(cell[i].site[j].scaledGenLike[k]);
            cell[i].site[j].scaledGenLike[k] = NULL;
            free(cell[i].site[j].genLikeDoublet[k]);
            cell[i].site[j].genLikeDoublet[k] = NULL;
            free(cell[i].site[j].scaledGenLikeDoublet[k]);
            cell[i].site[j].scaledGenLikeDoublet[k] = NULL;
          }
          free(cell[i].site[j].readCount);
          cell[i].site[j].readCount = NULL;
          free(cell[i].site[j].readCountDoublet);
          cell[i].site[j].readCountDoublet = NULL;
          free(cell[i].site[j].genLike);
          cell[i].site[j].genLike = NULL;
          free(cell[i].site[j].genLikeDoublet);
          cell[i].site[j].genLikeDoublet = NULL;
          free(cell[i].site[j].scaledGenLike);
          cell[i].site[j].scaledGenLike = NULL;
          free(cell[i].site[j].scaledGenLikeDoublet);
          cell[i].site[j].scaledGenLikeDoublet = NULL;
        }
        free(cell[i].site);
        cell[i].site = NULL;
      }
      free(cell);
      cell = NULL;
    }

    /* added code */
    if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
        altModel == FiniteMu) {
      if (doPrintTree == YES)
        fclose(fpTrees);

      fclose(fpTreeWithoutTrunk);

      fclose(fpTreeWithHealthyTip);

      if (dataSetNum == 0)
        fclose(fpCellNames);

      fclose(fpParallelMut);

      fclose(fpTrueSNVGenotypesAl);
      fclose(fpTrueSNVGenotypesNr);
      fclose(fpTrueSNVTernaryMatrix);

      for (i = 0; i < GetTumCellsTotalNum(numCells) + 2; i++) {
        fclose(rawDataFiles[i].fp);
        rawDataFiles[i].fp = NULL;
      }

      fclose(fpTrueSNVCov);

      fclose(fpEvolutionaryEvents);

      fclose(fpMpileupDataNM);
      fclose(fpMpileupDataM);
      if (meanDoubletRate > 0)
      {
        fclose(fpMpileupDataNU);
        fclose(fpMpileupDataU);
      }

      /* close file pointers for separate replicate file */
      if (doPrintSeparateReplicates == YES && doPrintFullGenotypes == YES)
        fclose(fpFullGenotypes);

      if (doPrintSNVgenotypes == YES)
        fclose(fpSNVgenotypes);

      fclose(fpSNVSitesNames);
      fclose(fpSNVAllelicSeqInfo);
      fclose(fpAdoStates);
      fclose(fpSizeFactors);

      for (i = 0; i < 2 * numCells + 1; ++i) {
        for (j = 0; j < numSites; ++j) {
          cellExtd[i].site[j].pastEvolutionaryEvents = NULL;

          if (i < numCells + 1)
          {
            free(cellExtd[i].site[j].observedAlleles);
          }
          cellExtd[i].site[j].observedAlleles = NULL;

          free(cellExtd[i].site[j].readCount);
          cellExtd[i].site[j].readCount = NULL;
        }
        free(cellExtd[i].site);
        cellExtd[i].site = NULL;
      }
      free(cellExtd);
      cellExtd = NULL;

      for (i = 0; i < numSites; i++) {
        free(allSites[i].alternateAlleles);
        allSites[i].alternateAlleles = NULL;
      }
      free(allSites);
      allSites = NULL;

      free(SNVsites);
      SNVsites = NULL;
      free(variantSites);
      variantSites = NULL;
      free(DefaultModelSites);
      DefaultModelSites = NULL;
      free(AltModelSites);
      AltModelSites = NULL;

      free(mutatedSites);
      mutatedSites = NULL;
      free(deletedSites);
      deletedSites = NULL;
      free(insertedSites);
      insertedSites = NULL;

      for (i = 0; i < 4; i++) {
        free(probNBDetails[i]);
        probNBDetails[i] = NULL;
      }

      for (i = 0; i < numThreads; i++)
      {
        if (rawDataPrinters[i].cellIndices != NULL)
        {
          free(rawDataPrinters[i].cellIndices);
          rawDataPrinters[i].cellIndices = NULL;
        }
        if (rawDataPrinters[i].files != NULL)
        {
          free(rawDataPrinters[i].files);
          rawDataPrinters[i].files = NULL;
        }
      }
    } else {
      if (doPrintSeparateReplicates == YES && doCoalescenceTree == YES) {
        if (doPrintTree == YES)
          fclose(fpTrees);

        if (doPrintTimes == YES)
          fclose(fpTimes);
      }
    }
    /* end added code */

    /* if(((dataSetNum*10)%numDataSets) == 0)
     PrintMemory(stderr); */
  } // replicate

  FreeTreeStructure(numNodes);

  if (cellNamesExtd != NULL) {
    for (i = 0; i < 2 * numCells + 1; i++) {
      free(cellNamesExtd[i]);
      cellNamesExtd[i] = NULL;
    }
    free(cellNamesExtd);
    cellNamesExtd = NULL;
  }

  if (!(numDataSets > 1 && fixedTree == YES)) {
    meanNumCA = cumNumCA / (double)numDataSets;
  }
  meanNumMU = cumNumMU / (double)numDataSets;
  meanNumSNVs = cumNumSNVs / (double)numDataSets;
  meanNumDEL = cumNumDEL / (double)numDataSets;

  /* added code */
  if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
        altModel == FiniteMu)) {
    meanNumCNLOH = cumNumCNLOH / (double)numDataSets;
  } else {
    meanNumINS = cumNumINS / (double)numDataSets;
  }
  /* end added code */
  /* original code */
  // meanNumCNLOH = cumNumCNLOH / (double)numDataSets;

  if (!(numDataSets > 1 && fixedTree == YES)) {
    meanTMRCA = cumTMRCA / (double)numDataSets;
  }

  varNumMU = (1.0 / (double)(numDataSets - 1)) *
             (cumNumMUSq - pow(cumNumMU, 2) / (double)numDataSets);
  varNumSNVs = (1.0 / (double)(numDataSets - 1)) *
               (cumNumSNVsSq - pow(cumNumSNVs, 2) / (double)numDataSets);
  varNumDEL = (1.0 / (double)(numDataSets - 1)) *
              (cumNumDELSq - pow(cumNumDEL, 2) / (double)numDataSets);
  /* added code */
  varNumINS = (1.0 / (double)(numDataSets - 1)) *
              (cumNumINSSq - pow(cumNumINS, 2) / (double)numDataSets);
  /* end added code */
  varTMRCA = (1.0 / (double)(numDataSets - 1)) *
             (cumTMRCASq - pow(cumTMRCA, 2) / (double)numDataSets);

  /* close file pointers */
  /* added code */
  if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
        altModel == FiniteMu)) {
    if (doPrintSeparateReplicates == NO) {
      if (doPrintSNVgenotypes == YES)
        fclose(fpSNVgenotypes);
      if (doPrintFullGenotypes == YES)
        fclose(fpFullGenotypes);
      if (doPrintSNVhaplotypes == YES)
        fclose(fpSNVhaplotypes);
      if (doPrintFullHaplotypes == YES)
        fclose(fpFullHaplotypes);
      if (doPrintTrueHaplotypes == YES)
        fclose(fpTrueHaplotypes);
      if (doNGS == YES)
        fclose(fpVCF);
      if (doPrintCATG == YES)
        fclose(fpCATG);
      if (doPrintTree == YES)
        fclose(fpTrees);
      if (doPrintTimes == YES)
        fclose(fpTimes);
    }
  } else {
    free(rawDataFiles);
    rawDataFiles = NULL;

    if (doPrintSeparateReplicates == NO) {
      if (doPrintFullGenotypes == YES)
        fclose(fpFullGenotypes);

      if (doPrintTimes == YES)
        fclose(fpTimes);
    }
  }
  /* end added code */
  /* original code */
  /*
  if (doPrintSeparateReplicates == NO) {
    if (doPrintSNVgenotypes == YES)
      fclose(fpSNVgenotypes);
    if (doPrintFullGenotypes == YES)
      fclose(fpFullGenotypes);
    if (doPrintSNVhaplotypes == YES)
      fclose(fpSNVhaplotypes);
    if (doPrintFullHaplotypes == YES)
      fclose(fpFullHaplotypes);
    if (doPrintTrueHaplotypes == YES)
      fclose(fpTrueHaplotypes);
    if (doNGS == YES)
      fclose(fpVCF);
    if (doPrintCATG == YES)
      fclose(fpCATG);
    if (doPrintTree == YES)
      fclose(fpTrees);
    if (doPrintTimes == YES)
      fclose(fpTimes);
  }
   */

  if (noisy > 0) {
    /* if (noisy < 2)
     PrintHeader (stdout);*/

    fprintf(stdout, "\n\n** Simulations finished at ");
    PrintDate(stdout);

    /* added code */
    if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu)) {
      PrintRunInformation(stdout);
      /* print run settings to file */
      PrintRunInformation(fpLog);
    } else {
      PrintRunInformationExtd(fpLog);
      if (noisy > 2) {
        PrintRunInformationExtd(stderr);
      }
    }
    fclose(fpLog);
    /* end added code */
    /*original code*/
    /*
    PrintRunInformation(stdout);
     */
    /*original code*/
    /* print run settings to file */
    /*
    PrintRunInformation(fpLog);
     */

    fprintf(stdout, "\n\nOutput files are in folder \"%s\":", resultsDir);

    if (doPrintTree == YES) {
      fprintf(stdout, "\n Tree printed to file \"%s\"", treeFile);
      if (doPrintSeparateReplicates == YES)
        fprintf(stdout, " in folder \"%s\"", treeDir);

      fprintf(stdout,
              "\n Tree without healthy root and healthy tip printed "
              "to file \"%s\"",
              treeWithoutTrunkFile);
      if (doPrintSeparateReplicates == YES)
        fprintf(stdout, " in folder \"%s\"", treeDir);
    }

    if (doPrintTimes == YES) {
      fprintf(stdout, "\n Times printed to file \"%s\"", timesFile);
      if (doPrintSeparateReplicates == YES)
        fprintf(stdout, " in folder \"%s\"", timesDir);
    }

    /* added code */
    if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu)) {
      if (doSimulateData == YES) {
        if (doPrintSNVgenotypes == YES) {
          fprintf(stdout, "\n SNV genotypes printed to file \"%s\"",
                  SNVgenotypesFile);
          if (doPrintSeparateReplicates == YES)
            fprintf(stdout, " in folder \"%s\"", SNVgenotypesDir);
        }

        if (doPrintFullGenotypes == YES) {
          fprintf(stdout, "\n Full genotypes printed to file \"%s\"",
                  fullGenotypesFile);
          if (doPrintSeparateReplicates == YES)
            fprintf(stdout, " in folder \"%s\"", fullGenotypesDir);
        }

        if (doPrintSNVhaplotypes == YES) {
          if (doPrintIUPAChaplotypes == YES)
            fprintf(stdout,
                    "\n SNV haplotypes (IUPAC codes) printed to file \"%s\"",
                    SNVhaplotypesFile);
          else
            fprintf(stdout, "\n SNV haplotypes printed to file \"%s\"",
                    SNVhaplotypesFile);
          if (doPrintSeparateReplicates == YES)
            fprintf(stdout, " in folder \"%s\"", SNVhaplotypesDir);
        }

        if (doPrintFullHaplotypes == YES) {
          if (doPrintIUPAChaplotypes == YES)
            fprintf(stdout,
                    "\n Full haplotypes (IUPAC codes) printed to file \"%s\"",
                    fullHaplotypesFile);
          else
            fprintf(stdout, "\n Full haplotypes printed to file \"%s\"",
                    fullHaplotypesFile);
          if (doPrintSeparateReplicates == YES)
            fprintf(stdout, " in folder \"%s\"", fullHaplotypesDir);
        }

        if (doPrintTrueHaplotypes == YES) {
          if (doPrintIUPAChaplotypes == YES)
            fprintf(stdout,
                    "\n True haplotypes (IUPAC codes) printed to file \"%s\"",
                    trueHaplotypesFile);
          else
            fprintf(stdout, "\n True haplotypes printed to file \"%s\"",
                    trueHaplotypesFile);
          if (doPrintSeparateReplicates == YES)
            fprintf(stdout, " in folder \"%s\"", trueHaplotypesDir);
        }

        if (doNGS == YES) {
          fprintf(stdout, "\n Genotype likelihoods printed to file \"%s\"",
                  VCFfile);
          if (doPrintSeparateReplicates == YES)
            fprintf(stdout, " in folder \"%s\"", VCFdir);
        }

        if (doPrintCATG == YES) {
          fprintf(stdout, "\n Read counts printed to file \"%s\"", CATGfile);
          if (doPrintSeparateReplicates == YES)
            fprintf(stdout, " in folder \"%s\"", CATGdir);
        }
      } else
        fprintf(stdout, "\n\n\nNo output data files");
    } else {
      if (doPrintSNVgenotypes == YES) {
        fprintf(stdout, "\n SNV genotypes printed to file \"%s\"",
                SNVgenotypesFile);
        if (doPrintSeparateReplicates == YES)
          fprintf(stdout, " in folder \"%s\"", SNVgenotypesDir);
      }

      if (doPrintFullGenotypes == YES) {
        fprintf(stdout, "\n Full genotypes printed to file \"%s\"",
                fullGenotypesFile);
        if (doPrintSeparateReplicates == YES)
          fprintf(stdout, " in folder \"%s\"", fullGenotypesDir);
      }

      fprintf(stdout,
              "\n Genotype likelihoods printed to file \"%s\" in folder \"%s\"",
              VCFfile, VCFdir);
    }
    /* end added code */
    /* original code */
    /*
    if (doSimulateData == YES) {
      if (doPrintSNVgenotypes == YES) {
        fprintf(stdout, "\n SNV genotypes printed to file \"%s\"",
                SNVgenotypesFile);
        if (doPrintSeparateReplicates == YES)
          fprintf(stdout, " in folder \"%s\"", SNVgenotypesDir);
      }

      if (doPrintFullGenotypes == YES) {
        fprintf(stdout, "\n Full genotypes printed to file \"%s\"",
                fullGenotypesFile);
        if (doPrintSeparateReplicates == YES)
          fprintf(stdout, " in folder \"%s\"", fullGenotypesDir);
      }

      if (doPrintSNVhaplotypes == YES) {
        if (doPrintIUPAChaplotypes == YES)
          fprintf(stdout,
                  "\n SNV haplotypes (IUPAC codes) printed to file \"%s\"",
                  SNVhaplotypesFile);
        else
          fprintf(stdout, "\n SNV haplotypes printed to file \"%s\"",
                  SNVhaplotypesFile);
        if (doPrintSeparateReplicates == YES)
          fprintf(stdout, " in folder \"%s\"", SNVhaplotypesDir);
      }

      if (doPrintFullHaplotypes == YES) {
        if (doPrintIUPAChaplotypes == YES)
          fprintf(stdout,
                  "\n Full haplotypes (IUPAC codes) printed to file \"%s\"",
                  fullHaplotypesFile);
        else
          fprintf(stdout, "\n Full haplotypes printed to file \"%s\"",
                  fullHaplotypesFile);
        if (doPrintSeparateReplicates == YES)
          fprintf(stdout, " in folder \"%s\"", fullHaplotypesDir);
      }

      if (doPrintTrueHaplotypes == YES) {
        if (doPrintIUPAChaplotypes == YES)
          fprintf(stdout,
                  "\n True haplotypes (IUPAC codes) printed to file \"%s\"",
                  trueHaplotypesFile);
        else
          fprintf(stdout, "\n True haplotypes printed to file \"%s\"",
                  trueHaplotypesFile);
        if (doPrintSeparateReplicates == YES)
          fprintf(stdout, " in folder \"%s\"", trueHaplotypesDir);
      }

      if (doNGS == YES) {
        fprintf(stdout, "\n Genotype likelihoods printed to file \"%s\"",
                VCFfile);
        if (doPrintSeparateReplicates == YES)
          fprintf(stdout, " in folder \"%s\"", VCFdir);
      }

      if (doPrintCATG == YES) {
        fprintf(stdout, "\n Read counts printed to file \"%s\"", CATGfile);
        if (doPrintSeparateReplicates == YES)
          fprintf(stdout, " in folder \"%s\"", CATGdir);
      }
    } else
      fprintf(stdout, "\n\n\nNo output data files");
      */
  }

#ifdef CHECK_MUT_DISTRIBUTION
  fprintf(stdout, "\n\nChecking distribution of mutations");
  fprintf(stdout,
          "\n Avg # sites with  0 mutant alelles = %6.2f  (Exp = %6.2f)",
          (double)MutCount[0] / numDataSets, numSites - expNumMU);
  for (i = 1; i < numCells; i++)
    fprintf(stdout,
            "\n Avg # sites with %2d mutant alelles = %6.2f  (Exp = %6.2f)", i,
            (double)MutCount[i] / numDataSets, theta / i);
  fprintf(stdout, "\n MeansumPos = %4.2f  (expected (L/2) = %4.2f)\n\n",
          meansumPos / (double)dataSetsWithSNVs, numSites / 2.0);
#endif

  /* free pending stuff */
  free(CommandLine);
  CommandLine = NULL;
  free(cellNames);
  cellNames = NULL;
  free(triMutationsCounter);
  triMutationsCounter = NULL;
  free(signatureWeight);
  signatureWeight = NULL;
  free(signatureID);
  signatureID = NULL;
  if (doUserGenome == YES) {
    free(paternalUserGenome);
    paternalUserGenome = NULL;
    free(maternalUserGenome);
    maternalUserGenome = NULL;
  }

  /* added code */
  if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
      altModel == FiniteMu) {
    free(threads);
    free(sitesPointsForThreads);
    sitesPointsForThreads = NULL;

    for (i = 0; i < numThreads; i++) {
      for (j = 0; j < 2 * numCells + 1; j++) {
        for (k = 0; k < threadStr[i].numSubSites; k++) {
          free(threadStr[i].data[j][k].pastEvolutionaryEvents);
          threadStr[i].data[j][k].pastEvolutionaryEvents = NULL;
        }

        free(threadStr[i].data[j]);
        threadStr[i].data[j] = NULL;
      }

      rawDataPrinters[i].rd = NULL;
      threadStr[i].rd = NULL;

      free(threadStr[i].hasParallelMutations);
      threadStr[i].hasParallelMutations = NULL;

      free(threadStr[i].data);
      threadStr[i].data = NULL;

      for (k = 0; k < threadStr[i].numSubSites; k++) {
        free(threadStr[i].summary[k]);
        threadStr[i].summary[k] = NULL;
        free(threadStr[i].parallelMutations[k]);
        threadStr[i].parallelMutations[k] = NULL;
      }
      free(threadStr[i].summary);
      threadStr[i].summary = NULL;
      free(threadStr[i].parallelMutations);
      threadStr[i].parallelMutations = NULL;

      free(threadStr[i].isSNV);
      threadStr[i].isSNV = NULL;
      free(threadStr[i].isDeleted);
      threadStr[i].isDeleted = NULL;
      free(threadStr[i].isInserted);
      threadStr[i].isInserted = NULL;
    }
    free(threadStr);
    threadStr = NULL;

    free(rawDataPrinters);
    rawDataPrinters = NULL;

    mainThreadRandDev = NULL;
    FreeMem(randDevs);

    free(allelicCov);
    allelicCov = NULL;

    free(allelicRawVar);
    allelicRawVar = NULL;

    free(sizeFactors);
    sizeFactors = NULL;

    if (fixedTree == YES) {
      for (i = 0; i < numNodes; i++) {
        free(treeNodes[i].name);
        treeNodes[i].name = NULL;
      }
      free(treeNodes);
      treeNodes = NULL;
    }

    if (numDataSets > 1 && fixedTree == YES) {
      for (i = 0; i < 4; i++) {
        free(probNBDetailsAll[i]);
        probNBDetailsAll[i] = NULL;
      }
    }

    FreeFiniteMuInDelSimulator(altModel == FiniteMuInDel ? 7 : 3);
  }
  /* end added code */

  secs = omp_get_wtime() - start;

  if (noisy > 0) {
    printf("\n\n_______________________________________________________________"
           "__");
    printf("\nTime processing: %G seconds", secs);
    printf("\nFinishing @ ");
    PrintDate(stdout);
    printf("\n");
  }

  //  pthread_exit(NULL);
  return 0;
}

/************************* MakeCoalescenceTree  ************************/
/*	Builds a genealogy  under the coalescent with exponential growth without
 * recombination */

void MakeCoalescenceTree(int numCells, int N, RandomDevice *rd) {
  int i, *activeNodes, numActiveNodes, nextAvailableNode, firstNode, secondNode,
      newNode, eventNum, period;
  double rateCA, timeCA, currentTime, eventTime;
  TreeNode *p, *q, *r;

  p = q = r = NULL;

  /* allocate space for tree */
  treeNodes = (TreeNode *)calloc(numNodes + 1, sizeof(TreeNode));
  if (!treeNodes) {
    fprintf(stderr, "Could not allocate nodes (%ld)\n",
            numNodes * sizeof(TreeNode));
    exit(-1);
  }
  activeNodes = (int *)calloc(numNodes + 1, sizeof(int));
  if (!activeNodes) {
    fprintf(stderr, "Could not allocate activeNodes (%ld)\n",
            numNodes * sizeof(int));
    exit(-1);
  }

  /* set everything to null */
  for (i = 0; i < numNodes; i++) {
    treeNodes[i].left = NULL;
    treeNodes[i].right = NULL;
    treeNodes[i].anc = NULL;
    treeNodes[i].time = 0.0;
    treeNodes[i].length = 0.0;
    treeNodes[i].index = 0;
    treeNodes[i].label = 0;
    treeNodes[i].isHealthyRoot = NO;
    treeNodes[i].isHealthyTip = NO;
    treeNodes[i].name = (char *)calloc(MAX_NAME, sizeof(char));
    if (!treeNodes[i].name) {
      fprintf(stderr, "Could not allocate the treeNodes[i].name structure\n");
      exit(-1);
    }
  }

  /* set up initial nodes */
  numActiveNodes = 0;
  for (i = 0; i < numCells; i++) {
    p = treeNodes + i;
    activeNodes[numActiveNodes] = i;
    p->index = i;
    numActiveNodes++;
  }
  nextAvailableNode = numActiveNodes;

  /* make coalescence tree */
  eventNum = 0;
  period = 1;
  currentTime = 0.0;
  eventTime = 0.0;
  while (numActiveNodes > 1) {
    rateCA = numActiveNodes * (numActiveNodes - 1) / 2.0;

    if (doDemographics == YES) {
      periodGrowth[period] =
          -log(Nend[period] / (double)Nbegin[period]) /
          (double)(cumDuration[period] - cumDuration[period - 1]);
      if (isnan(periodGrowth[period]) == YES) {
        fprintf(stderr, "\nERROR: period growth (%f) is NaN",
                periodGrowth[period]);
        fprintf(
            stderr,
            "\n       This might suggest that the growth rate is too negative");
        fprintf(stderr,
                "\n       and the coalescent time is therefore infinite.");
        fprintf(stderr, "\n       Try a smaller value");
        exit(-1);
      }

      timeCA = log(1 + RandExponential(rateCA, rd) * periodGrowth[period] *
                           2.0 * Nbegin[period] *
                           exp(-periodGrowth[period] *
                               (currentTime - cumDuration[period - 1]))) /
               periodGrowth[period]; /* time is kept in generations (scaled by
                                        2N)  */

      /*	When growth rate is very negative, coalescent time may be
       infinite this results in log (-x) => timCA = NaN. If this not the last
       period just jump to the next. If this is the last period, we have to exit
       the program
       */
      if (isnan(timeCA) == YES) {
        if (period < numPeriods) {
          currentTime = cumDuration[period];
          period++;
          periodGrowth[period] = -log(Nend[period] / (double)Nbegin[period]) /
                                 (double)cumDuration[period];
          continue;
        } else {
          fprintf(stderr, "\nERROR: Coalescent time (%f) is infinite ", timeCA);
          fprintf(stderr, "\n       This might suggest that the exponential "
                          "growth rate is too negative");
          fprintf(stderr,
                  "\n       and the coalescent time is therefore infinite.");
          fprintf(stderr, "\n       Try a smaller value");
          exit(-1);
        }
      }
    } else {
      timeCA = RandExponential(rateCA, rd) * 2.0 *
               N; /* time is kept in generations */
      if (doExponential == YES) {
        timeCA = log(exp(growthRate * currentTime) + growthRate * timeCA) /
                     growthRate -
                 currentTime;

        /*	When growth rate is very negative, coalescent time may be
         infinite this results in log (-x) => timCA = NaN. We have to exit the
         program
         */
        if (isnan(timeCA) == YES) {
          fprintf(stderr, "\nERROR: Coalescent time (%f) is infinite ", timeCA);
          fprintf(stderr, "\n       This might suggest that the exponential "
                          "growth rate is too negative");
          fprintf(stderr,
                  "\n       and the coalescent time is therefore infinite.");
          fprintf(stderr, "\n       Try a smaller value");
          exit(-1);
        }
      }
    }

    eventTime = timeCA;

    /*	if this period is not the last one and if the event time is outside the
       current interval, update period and start again */
    if (doDemographics == YES && period < numPeriods &&
        (currentTime + eventTime) > cumDuration[period]) {
      currentTime = cumDuration[period];
      period++;
      periodGrowth[period] =
          -log(Nend[period] - Nbegin[period]) / (double)cumDuration[period];
      continue;
    }
    numCA++;

    /* set time */
    currentTime += eventTime;

    if (noisy > 3)
      fprintf(stderr, "\nEvent %3d :  rateCA = %lf, currentTime = %lf\n",
              ++eventNum, rateCA, currentTime);

    /* figure out which two nodes will coalesce */
    firstNode = numActiveNodes * RandUniform(rd);
    if (firstNode >= numActiveNodes) {
      fprintf(stderr, "\n\nERROR: firstNode out of range!\n");
      exit(-1);
    }
    do {
      secondNode = numActiveNodes * RandUniform(rd);
    } while (firstNode == secondNode);

    newNode = nextAvailableNode;
    if (noisy > 2)
      fprintf(stderr, "\n  Coalescence involving %d and %d to create node %d",
              activeNodes[firstNode] + 1, activeNodes[secondNode] + 1,
              newNode + 1);

    /* readjust MRCA, pointers and active sites */
    p = treeNodes + activeNodes[firstNode];
    q = treeNodes + activeNodes[secondNode];
    r = treeNodes + newNode;
    r->index = nextAvailableNode;
    r->left = p;
    r->right = q;
    p->anc = r;
    q->anc = r;
    r->time = currentTime;
    p->length = r->time - p->time;
    q->length = r->time - q->time;
    p->branchLength = p->length * mutationRate;
    q->branchLength = q->length * mutationRate;

    /* readjust active nodes */
    activeNodes[firstNode] = newNode;
    activeNodes[secondNode] = activeNodes[numActiveNodes - 1];
    numActiveNodes--;
    nextAvailableNode++;
    if (nextAvailableNode >= numNodes) {
      fprintf(stderr, "\n\nERROR: Too many nodes!\n");
      exit(-1);
    }
  } /* coalescent cell tree finished */

  TMRCA = currentTime /
          (2.0 * N); /* TMRCA is now is back in 2N coalescent generations */
  coalTreeMRCA = r;

  /* Connect the coalescent cell MRCA node with the healthy ancestral cell */
  if (noisy > 2)
    fprintf(stderr, "\n>> Adding healthy root ... ");
  healthyRoot = treeNodes + nextAvailableNode;
  healthyRoot->index = nextAvailableNode;
  healthyRoot->left = coalTreeMRCA;
  coalTreeMRCA->anc = healthyRoot;
  healthyRoot->time =
      currentTime + (currentTime * transformingBranchLengthRatio);
  healthyRoot->length = 0;
  healthyRoot->isHealthyRoot = YES;
  coalTreeMRCA->length = healthyRoot->time - coalTreeMRCA->time;
  coalTreeMRCA->branchLength = coalTreeMRCA->length * mutationRate;
  if (noisy > 2)
    fprintf(stderr, "DONE");

  nextAvailableNode++;

  /* Connect the healthy ancestral cell with the tip healthy cell*/
  if (noisy > 2)
    fprintf(stderr, "\n>> Adding healthy tip ... ");
  healthyTip = treeNodes + nextAvailableNode;
  healthyTip->index = nextAvailableNode;
  healthyTip->left = NULL;
  healthyTip->right = NULL;
  healthyTip->anc = healthyRoot;
  healthyRoot->right = healthyTip;
  healthyTip->time = 0;
  healthyTip->length =
      (healthyRoot->time - healthyTip->time) * healthyTipBranchLengthRatio;
  healthyTip->branchLength = healthyTip->length * mutationRate;
  healthyTip->isHealthyTip = YES;
  if (noisy > 2)
    fprintf(stderr, "DONE");

  /* Relabel nodes on tree (label tips with indexes) */
  if (noisy > 2)
    fprintf(stderr, "\n>> Relabeling nodes on tree ... ");
  tipLabel = numCells;
  intLabel = numCells + 1;
  RelabelNodes(healthyRoot);
  if (noisy > 2)
    fprintf(stderr, "DONE");

  free(activeNodes);
  activeNodes = NULL;
}

/********************************** ReadUserGenome
 * ***********************************/
/* Reads a FASTA with a diploid reference genome */

void ReadUserGenome(FILE *fp) {
  char *header;

  /* allocate space for reading */
  header = (char *)calloc(MAX_LINE, sizeof(char));
  if (!header) {
    fprintf(stderr, "Could not allocate the header string\n");
    exit(-1);
  }

  maternalUserGenome = (char *)calloc(MAX_GENOME, sizeof(char));
  if (!maternalUserGenome) {
    fprintf(stderr, "Could not allocate the maternalUserGenome structure\n");
    exit(-1);
  }

  paternalUserGenome = (char *)calloc(MAX_GENOME, sizeof(char));
  if (!paternalUserGenome) {
    fprintf(stderr, "Could not allocate the paternalUserGenome structure\n");
    exit(-1);
  }

  /* open user genome */
  if ((fpUserGenome = fopen(userGenomeFile, "r")) == NULL) {
    fprintf(stderr, "\nERROR: Can't open user genome file \"%s\"\n",
            userGenomeFile);
    PrintUsage(stderr);
  }

  /* read the FASTA file */
  fgets(header, MAX_LINE, fpUserGenome);
  fgets(maternalUserGenome, MAX_GENOME, fpUserGenome);
  if (maternalUserGenome == NULL) {
    fprintf(stderr, "\nERROR: Can't read the user maternal genome");
    PrintUsage(stderr);
  }

  fgets(header, MAX_LINE, fpUserGenome);
  fgets(paternalUserGenome, MAX_GENOME, fpUserGenome);
  if (paternalUserGenome == NULL) {
    fprintf(stderr, "\nERROR: Can't read the user paternal genome");
    PrintUsage(stderr);
  }

  /* check lengths */
  if (strlen(maternalUserGenome) != strlen(paternalUserGenome)) {
    fprintf(
        stderr,
        "\nERROR: maternal and paternal user genomes have different lengths");
    PrintUsage(stderr);
  }

  numSites = (int)strlen(maternalUserGenome) - 1;

  /* realloc memory */
  maternalUserGenome =
      (char *)realloc(maternalUserGenome, numSites * sizeof(char));
  paternalUserGenome =
      (char *)realloc(paternalUserGenome, numSites * sizeof(char));

  fprintf(stderr, "\nnumsites = %d", numSites);
  fprintf(stderr, "\nuser maternal genome = %s", maternalUserGenome);
  fprintf(stderr, "\nuser paternal genome = %s", paternalUserGenome);

  free(header);
  header = NULL;
  fclose(fpUserGenome);
}

/********************************** ReadUserTree
 * ***********************************/
/* Reads a rooted phyogenetic tree with branch lengths in Newick format */
static void ReadUserTree(FILE *fp) {
  int taxonNumber, numDigits;
  int nextAvailableNode;
  int i, j;
  int digit, digit1, digit2, digit3;
  char *tempString;
  TreeNode *p, *q;

  p = q = NULL;

  /* allocate space for reading the tree string */
  treeString = (char *)calloc(MAX_LINE, sizeof(char));
  if (!treeString) {
    fprintf(stderr, "Could not allocate the treeString structure\n");
    exit(-1);
  }

  tempString = (char *)calloc(MAX_NAME, sizeof(char));
  if (!tempString) {
    fprintf(stderr, "Could not allocate the tempString structure\n");
    exit(-1);
  }

  /* allocate space for tree */
  treeNodes = (TreeNode *)calloc(numNodes + 1, sizeof(TreeNode));
  if (!treeNodes) {
    fprintf(stderr, "Could not allocate nodes (%ld)\n",
            numNodes * sizeof(TreeNode));
    exit(-1);
  }

  /* allocate space for cellNames */
  cellNames = (char **)calloc(numCells + 1, sizeof(char *));
  if (!cellNames) {
    fprintf(stderr, "Could not allocate cellNames (%ld)\n",
            numCells * sizeof(char *));
    exit(-1);
  }
  for (i = 0; i < numCells; i++) {
    cellNames[i] = (char *)calloc(MAX_NAME, sizeof(char));
    if (!cellNames[i]) {
      fprintf(stderr, "Could not allocate cellNames[%d] (%ld)\n", i,
              MAX_NAME * sizeof(char));
      exit(-1);
    }
  }

  /* initialize tree nodes */
  for (i = 0; i < numNodes; i++) {
    treeNodes[i].left = NULL;
    treeNodes[i].right = NULL;
    treeNodes[i].anc = NULL;
    treeNodes[i].time = 0.0;
    treeNodes[i].length = 0.0;
    treeNodes[i].index = -1;
    treeNodes[i].label = -1;
    treeNodes[i].isHealthyRoot = NO;
    treeNodes[i].isHealthyTip = NO;
    treeNodes[i].name = (char *)calloc(MAX_NAME, sizeof(char));
    if (!treeNodes[i].name) {
      fprintf(stderr, "Could not allocate the treeNodes[i].name structure\n");
      PrintUsage(stderr);
    }
  }

  /* open treefile */
  if ((fpUserTree = fopen(userTreeFile, "r")) == NULL) {
    fprintf(stderr, "\nERROR: Can't open user treefile \"%s\"\n", userTreeFile);
    PrintUsage(stderr);
  }

  /* read the tree string */
  fgets(treeString, MAX_LINE, fpUserTree);
  if (treeString == NULL) {
    fprintf(stderr, "\nERROR: Can't read tree string");
    PrintUsage(stderr);
  }
  fclose(fpUserTree);

  CheckTree(treeString);

  /* build the tree linked list */
  numCells = 0;
  nextAvailableNode = 0;
  i = 0;
  do {
    if (treeString[i] == '(') {
      if (nextAvailableNode == 0) {
        p = treeNodes + nextAvailableNode;
        coalTreeMRCA = p;
        p->index = nextAvailableNode;
        nextAvailableNode++;
      } else {
        q = treeNodes + nextAvailableNode;
        if (p->left == NULL) {
          p->left = q;
          q->anc = p;
        } else {
          p->right = q;
          q->anc = p;
        }
        p = q;
        p->index = nextAvailableNode;
        nextAvailableNode++;
      }
    } else if (treeString[i] == ')') {
      q = p->anc;
      p = q;
    } else if (treeString[i] == ',') {
      q = p->anc;
      p = q;
    } else if (treeString[i] == ':') {
      i++;
      for (j = 0; j < 10; j++)
        tempString[j] = ' ';
      j = 0;
      while (treeString[i] != ',' && treeString[i] != ')' &&
             treeString[i] != ':' && treeString[i] != ';' &&
             treeString[i] != '(')
        tempString[j++] = treeString[i++];
      i--;
      sscanf(tempString, "%lf", &p->branchLength);

      if (p->left == NULL)
        p->numLeftDescendants = 1;

      if (p->right == NULL)
        p->numRightDescendants = 1;

      if (p->left != NULL && p->right != NULL) {
        p->numLeftDescendants +=
            p->left->numLeftDescendants + p->left->numRightDescendants;
        p->numRightDescendants +=
            p->right->numLeftDescendants + p->right->numRightDescendants;

        p->meanDepth = ((p->left->branchLength + p->left->meanDepth) *
                            p->numLeftDescendants +
                        (p->right->branchLength + p->right->meanDepth) *
                            p->numRightDescendants) /
                       (p->numLeftDescendants + p->numRightDescendants);
      }
    } else
    /* now read taxon names*/
    {
      if (/* DISABLES CODE */ (0)) /* let's read every name as a word */
      /* if (isdigit(treeString[i]) == YES)*/ /* taxon names are numbers */
      {
        numDigits = 0;
        digit1 = digit2 = digit3 = 0;
        while (treeString[i] != ',' && treeString[i] != ')' &&
               treeString[i] != ':') {
          if (treeString[i] == '0')
            digit = 0;
          else if (treeString[i] == '1')
            digit = 1;
          else if (treeString[i] == '2')
            digit = 2;
          else if (treeString[i] == '3')
            digit = 3;
          else if (treeString[i] == '4')
            digit = 4;
          else if (treeString[i] == '5')
            digit = 5;
          else if (treeString[i] == '6')
            digit = 6;
          else if (treeString[i] == '7')
            digit = 7;
          else if (treeString[i] == '8')
            digit = 8;
          else if (treeString[i] == '9')
            digit = 9;
          if (numDigits == 0)
            digit1 = digit;
          else if (numDigits == 1)
            digit2 = digit;
          else if (numDigits == 2)
            digit3 = digit;
          numDigits++;
          i++;
        }
        i--;
        if (numDigits == 1)
          taxonNumber = digit1;
        else if (numDigits == 2)
          taxonNumber = 10 * digit1 + digit2;
        else if (numDigits == 3)
          taxonNumber = 100 * digit1 + 10 * digit2 + digit3;
      } else /* taxon names are letters */
      {
        taxonName = (char *)calloc(MAX_NAME, sizeof(char));
        if (!taxonName) {
          fprintf(stderr, "Could not allocate the taxonName structure\n");
          exit(-1);
        }
        taxonNamesAreChars = YES;
        j = 0;
        while (treeString[i] != ',' && treeString[i] != ')' &&
               treeString[i] != ':')
          taxonName[j++] = treeString[i++];
        i--;
      }
      q = treeNodes + nextAvailableNode;
      if (p->left == NULL) {
        p->left = q;
        q->anc = p;
      } else {
        p->right = q;
        q->anc = p;
      }
      strcpy(q->name, taxonName);
      strcpy(cellNames[numCells], taxonName);
      p = q;
      p->index = nextAvailableNode;
      nextAvailableNode++;
      numCells++;
      free(taxonName);
      taxonName = NULL;
    }
    i++;
  } while (treeString[i] != ';');

  /* Connect the coalescent cell MRCA node with the outgroup ancestral cell */
  p = coalTreeMRCA;
  if (p->left != NULL && p->right != NULL) {
    p->numLeftDescendants +=
        p->left->numLeftDescendants + p->left->numRightDescendants;
    p->numRightDescendants +=
        p->right->numLeftDescendants + p->right->numRightDescendants;

    p->meanDepth =
        ((p->left->branchLength + p->left->meanDepth) * p->numLeftDescendants +
         (p->right->branchLength + p->right->meanDepth) *
             p->numRightDescendants) /
        (p->numLeftDescendants + p->numRightDescendants);
  }

  if (noisy > 2)
    fprintf(stderr, "\n>> Adding healthy root to user tree... ");
  healthyRoot = treeNodes + nextAvailableNode;
  healthyRoot->index = nextAvailableNode;
  healthyRoot->left = coalTreeMRCA;
  healthyRoot->length = 0;
  healthyRoot->isHealthyRoot = YES;
  healthyRoot->name = "healthyroot";
  coalTreeMRCA->anc = healthyRoot;
  coalTreeMRCA->branchLength =
      coalTreeMRCA->meanDepth * transformingBranchLengthRatio;
  if (noisy > 2)
    fprintf(stderr, "DONE");

  nextAvailableNode++;

  /* Connect the healthy ancestral cell with the tip healthy cell*/
  if (noisy > 2)
    fprintf(stderr, "\n>> Adding healthy tip to user tree... ");
  healthyTip = treeNodes + nextAvailableNode;
  healthyTip->left = NULL;
  healthyTip->right = NULL;
  healthyTip->anc = healthyRoot;
  healthyRoot->right = healthyTip;
  healthyTip->time = 0;
  healthyTip->branchLength =
      (coalTreeMRCA->meanDepth + coalTreeMRCA->branchLength) *
      healthyTipBranchLengthRatio;
  healthyTip->isHealthyTip = YES;
  healthyTip->name = "healthy";
  healthyTip->index = nextAvailableNode;
  if (noisy > 2)
    fprintf(stderr, "DONE");

  /* Relabel nodes on tree (label tips with consecutive indexes) */
  if (noisy > 2)
    fprintf(stderr, "\n>> Relabeling nodes on tree ... ");

  TipNodeNum = 0;
  IntNodeNum = numCells + 1;
  RelabelUserTree(healthyRoot);

  if (noisy > 2)
    fprintf(stderr, "DONE");

  numNodes = nextAvailableNode + 1;

  /* reallocate tree nodes */
  treeNodes = (TreeNode *)realloc(treeNodes, numNodes * sizeof(TreeNode));
  if (treeNodes == NULL) {
    fprintf(stderr, "Could not reallocate treeNodes (%ld bytes)\n",
            numNodes * sizeof(TreeNode));
    PrintUsage(stderr);
  }

  fprintf(fpLog, "\n\nUser-defined tree: ");
  PrintTree(healthyRoot, fpLog);

  free(tempString);
  tempString = NULL;
}

/***************** CheckTree *******************/
/*  Checks for an error in the input tree string
        Prints message when there is an error in the
        tree and prints correct part of the tree */
static void CheckTree(char *treeStr) {
  int i, k;
  int thereIsError, numLeftPar, numRightPar;

  thereIsError = NO;
  numLeftPar = 1;
  numRightPar = 0;
  i = 1;

  do {
    if (treeStr[i] == '(') {
      if (treeStr[i - 1] != '(' && treeStr[i - 1] != ',')
        thereIsError = YES;
      if (treeStr[i + 1] != '(' && isalnum(treeStr[i + 1]) == NO)
        thereIsError = YES;
      numLeftPar++;
    } else if (treeStr[i] == ')') {
      if (treeStr[i - 1] != ')' && isalnum(treeStr[i - 1]) == NO)
        thereIsError = YES;
      if (treeStr[i + 1] != ')' && treeStr[i + 1] != ',' &&
          treeStr[i + 1] != ':')
        thereIsError = YES;
      if (treeStr[i + 1] == ';' && i == strlen(treeStr) - 3)
        thereIsError = NO;
      numRightPar++;
    } else if (treeStr[i] == ',') {
      if (treeStr[i - 1] != ')' && isalnum(treeStr[i - 1]) == NO)
        thereIsError = YES;
      if (treeStr[i + 1] != '(' && isalnum(treeStr[i + 1]) == NO)
        thereIsError = YES;
    } else if (treeStr[i] == ':') {
      if (treeStr[i - 1] != ')' && isalnum(treeStr[i - 1]) == NO)
        thereIsError = YES;
      if (isdigit(treeStr[i + 1]) == NO)
        thereIsError = YES;
    }

    if (thereIsError == YES) {
      fprintf(stderr, "\nERROR: There is something wrong in the tree\n");
      for (k = 0; k <= i; k++)
        fprintf(stderr, "%c", treeStr[k]);
      fprintf(stderr, " <- HERE");
      PrintUsage(stderr);
    }

    i++;

  } while (treeStr[i] != ';');

  if (numLeftPar != numRightPar) {
    fprintf(stderr, "Tree seems unbalanced (%d left and %d right parentheses)",
            numLeftPar, numRightPar);
    PrintUsage(stderr);
  }
  if (strrchr(treeString, ':') == NULL) {
    fprintf(stderr, "Tree does not have branch lengths");
    PrintUsage(stderr);
  }
}

/****************** RelabelUserTree **************************/
/*	Reindex nodes, count nodes to tips,
        and initialize components in the user tree */

void RelabelUserTree(TreeNode *p) {
  if (p != NULL) {
    RelabelUserTree(p->left);
    RelabelUserTree(p->right);
    if (p->left == NULL && p->right == NULL)
      p->index = TipNodeNum++;
    else
      p->index = IntNodeNum++;
    p->label = p->index;
  }
}

/**************** RelabelNodes **************/
/*	Relabel nodes so all tips are consecutive */

void RelabelNodes(TreeNode *p) {
  if (p != NULL) {
    RelabelNodes(p->left);
    RelabelNodes(p->right);
    if (p->left == NULL && p->right == NULL &&
        p->isHealthyTip == NO) /* is cell tip */
      p->label = p->index;
    else if (p->left == NULL && p->right == NULL &&
             p->isHealthyTip == YES) /* is outgroup tip */
      p->label = tipLabel++;
    else
      p->label = intLabel++; /* is internal */
  }
}

/**************** MakeTreeNonClock **************/
/*	Introduce rate variation among lineages using a gamma rate  */

void MakeTreeNonClock(TreeNode *p, RandomDevice *rd) {
  if (p != NULL) {
    p->branchLength *= RandGamma(alphaBranches, rd) / alphaBranches;
    MakeTreeNonClock(p->left, rd);
    MakeTreeNonClock(p->right, rd);
  }
}

/********************************** InitializeGenomes
 * ***********************************/
/* Initialize all genomes with the reference states  */
void InitializeGenomes(TreeNode *p, RandomDevice *rd) {
  int i, cell, anccell, site;
  double ran;

  if (p != NULL) {
    cell = p->label;
    if (alphabet == DNA) {
      if (p->isHealthyRoot == YES) {
        if (doUserGenome == YES) {
          for (site = 0; site < numSites;
               site++) /* initialize genome with the user provided sequences */
          {
            data[MATERNAL][cell][site] = WhichNucChar(maternalUserGenome[site]);
            data[PATERNAL][cell][site] = WhichNucChar(paternalUserGenome[site]);
            allSites[site].referenceAllele =
                data[MATERNAL][cell]
                    [site]; // then allSites[site].referenceAllele hosts the
                            // reference genome
          }
        } else if (doGeneticSignatures ==
                   NO) /* initialize genome with mononucleotide frequencies */
        {
          for (site = 0; site < numSites; site++) {
            ran = RandUniform(rd);
            for (i = 0; i < 4; i++) {
              if (ran <= cumfreq[i]) {
                allSites[site].referenceAllele = data[MATERNAL][cell][site] =
                    data[PATERNAL][cell][site] = i;
                // then allSites[site].referenceAllele hosts the reference
                // genome
                break;
              }
            }
          }
        } else /* initialize genome with trinucleotide frequencies */
        {
          SimulateTriNucFreqGenome(cell, rd);
        }
      } else {
        anccell = p->anc->label;
        for (site = 0; site < numSites; site++) {
          data[MATERNAL][cell][site] = data[MATERNAL][anccell][site];
          data[PATERNAL][cell][site] = data[PATERNAL][anccell][site];
        }
      }
    } else /* binary alphabet */
      for (site = 0; site < numSites; site++)
        data[MATERNAL][cell][site] = data[PATERNAL][cell][site] = 0;

    InitializeGenomes(p->left, rd);
    InitializeGenomes(p->right, rd);
  }
}

/* added code */
/* Initialize all genomes with the reference states,
 * namely A_A, C_C, G_G or T_T. */
void InitializeGenomesExtd(TreeNode *p, RandomDevice *rd) {
  int i, cell, anccell, site, threadIndex, relaSite;
  double ran;

  if (p != NULL) {
    cell = p->label;
    if (p->isHealthyRoot == YES) {
      for (i = 0; i < numThreads; i++) {
        memset(threadStr[i].hasParallelMutations, 0,
               threadStr[i].numSubSites * sizeof(int));
        memset(threadStr[i].isSNV, 0, threadStr[i].numSubSites * sizeof(int));
        memset(threadStr[i].isDeleted, 0, threadStr[i].numSubSites * sizeof(int));
        memset(threadStr[i].isInserted, 0, threadStr[i].numSubSites * sizeof(int));
      }

      if (doUserGenome == YES) {
        /* initialize genome with the user provided sequences */
        for (site = 0; site < numSites; site++) {
          threadIndex = GetThreadIndex(site);
          relaSite = site - threadStr[threadIndex].start;

          if (maternalUserGenome[site] != paternalUserGenome[site]) {
            fprintf(stderr,
                    "\nWARNING! The two nucleotides in cell %d at "
                    "site %d are not the same: %c and %c. Choose %c "
                    "as the reference allele.\n",
                    cell, site, maternalUserGenome[site],
                    paternalUserGenome[site], maternalUserGenome[site]);
          }
          // then allSites[site].referenceAllele hosts the reference genome
          allSites[site].referenceAllele = maternalUserGenome[site];

          memset(threadStr[threadIndex].summary[relaSite], 0, 5 * sizeof(int));
          memset(threadStr[threadIndex].parallelMutations[relaSite], 0,
                 (numCells + 2) * sizeof(unsigned));
          threadStr[threadIndex].data[cell][relaSite].code =
              WhichRefGenotype(allSites[site].referenceAllele,
                               altModel == FiniteMuInDel ? YES : NO);
          threadStr[threadIndex].data[cell][relaSite].numExistingAlleles = 2;
          threadStr[threadIndex].data[cell][relaSite].alleles =
              GetValIntToIntArr(
                  genotypeHashMapIntToIntArr,
                  threadStr[threadIndex].data[cell][relaSite].code);
          memset(threadStr[threadIndex]
                     .data[cell][relaSite]
                     .pastEvolutionaryEvents,
                 0, 5 * sizeof(unsigned int));
        }
      } else {
        for (site = 0; site < numSites; site++) {
          threadIndex = GetThreadIndex(site);
          relaSite = site - threadStr[threadIndex].start;

          ran = RandUniform(rd);
          for (i = 0; i < 4; i++) {
            if (ran <= cumfreq[i]) {
              // then allSites[site].referenceAllele hosts the reference genome
              allSites[site].referenceAllele = i;

              memset(threadStr[threadIndex].summary[relaSite], 0,
                     5 * sizeof(int));
              memset(threadStr[threadIndex].parallelMutations[relaSite], 0,
                     (numCells + 2) * sizeof(unsigned));
              threadStr[threadIndex].data[cell][relaSite].code =
                  WhichRefGenotype(i, altModel == FiniteMuInDel ? YES : NO);
              threadStr[threadIndex].data[cell][relaSite].numExistingAlleles =
                  2;
              threadStr[threadIndex].data[cell][relaSite].alleles =
                  GetValIntToIntArr(
                      genotypeHashMapIntToIntArr,
                      threadStr[threadIndex].data[cell][relaSite].code);
              memset(threadStr[threadIndex]
                         .data[cell][relaSite]
                         .pastEvolutionaryEvents,
                     0, 5 * sizeof(unsigned int));
              break;
            }
          }
        }
      }
    } else {
      anccell = p->anc->label;
      for (site = 0; site < numSites; site++) {
        threadIndex = GetThreadIndex(site);
        relaSite = site - threadStr[threadIndex].start;

        threadStr[threadIndex].data[cell][relaSite].code =
            threadStr[threadIndex].data[anccell][relaSite].code;
        threadStr[threadIndex].data[cell][relaSite].numExistingAlleles =
            threadStr[threadIndex].data[anccell][relaSite].numExistingAlleles;
        threadStr[threadIndex].data[cell][relaSite].alleles =
            threadStr[threadIndex].data[anccell][relaSite].alleles;
        memset(
            threadStr[threadIndex].data[cell][relaSite].pastEvolutionaryEvents,
            0, 5 * sizeof(unsigned int));
      }
    }

    InitializeGenomesExtd(p->left, rd);
    InitializeGenomesExtd(p->right, rd);
  }
}
/* end added code */

/********************************** SimulateTriNucFreqGenome
 * ***********************************/
/*
         Simulate a homozygous diploid genome with the trinucleotide frequencies
   of the human genome

         We simulate first a random trinucleotide.
         Then we go site by site taking into account the last two letters  of
   the previous trinucleotide.  if the first trinucleotide is, for example, CAT
         then we simulate trinucleotide starting in AT_  (i.e. ATA, ATC, ATG or
   ATT), and keep going  taking into account the last two letters of the new
   trinucleotide, and so on.
*/

void SimulateTriNucFreqGenome(int cell, RandomDevice *rd) {
  int chosenTriNucleotide, nextNucleotide;
  int k, n1, n2, n3, rest, site;
  double *prob4, sum;

  /* memory allocations */
  prob4 = (double *)calloc(4, sizeof(double));
  if (!prob4) {
    fprintf(stderr, "Could not allocate the prob4 vector\n");
    exit(-1);
  }

  /* choose first trinucleotide */
  chosenTriNucleotide = ChooseUniformState(triNucFreq, rd);

  /* find bases of the selected trinucleotide */
  n1 = chosenTriNucleotide / 16;
  rest = chosenTriNucleotide % 16;
  n2 = rest / 4;
  n3 = rest % 4;
  // fprintf (stderr, "\n%2d %c%c%c ", chosenTriNucleotide, WhichNuc(n1),
  // WhichNuc(n2), WhichNuc(n3));

  site = 0;
  data[MATERNAL][cell][site] = data[PATERNAL][cell][site] =
      allSites[site].referenceAllele = n1;
  site++;
  data[MATERNAL][cell][site] = data[PATERNAL][cell][site] =
      allSites[site].referenceAllele = n2;
  site++;
  data[MATERNAL][cell][site] = data[PATERNAL][cell][site] =
      allSites[site].referenceAllele = n3;

  /* fill the rest of the genome */
  /* choose next nucleotide given the last two bases of the previous
   * trinucleotide  */
  for (site = 3; site < numSites; site++) {
    /* normalize frequencies given the last two bases */
    sum = 0;
    for (k = 0; k < 4; k++)
      sum += triNucFreq[trinuc(n2, n3, k)];
    for (k = 0; k < 4; k++)
      prob4[k] = triNucFreq[trinuc(n2, n3, k)] / sum;

    nextNucleotide = ChooseUniformState(prob4, rd);
    data[MATERNAL][cell][site] = data[PATERNAL][cell][site] =
        allSites[site].referenceAllele = nextNucleotide;

    /* move downstream one position */
    n1 = n2;
    n2 = n3;
    n3 = nextNucleotide;
  }

  free(prob4);
  prob4 = NULL;
}

/********************************** RecordTriNucObservation
 * ***********************************/
/* add a new occurrence of a trinucleotide along the genome */
void RecordTriNucObservation(TriNucStr *trin, int site) {
  if (trin->numAvailablePositions == trin->tempLength) {
    trin->tempLength *= 2;
    trin->position =
        (int *)realloc(trin->position, trin->tempLength * sizeof(int));
    if (!trin->position) {
      fprintf(
          stderr,
          "Could not reallocate trinucleotide position array when adding\n");
      exit(-1);
    }
  }
  trin->position[trin->numAvailablePositions++] = site;
}

/********************************** AddGermlineVariation
 * ***********************************/
/* Introduce germline SNPs in the healthy root genome for DNA or 0/1 data  */
void AddGermlineVariation(RandomDevice *rd) {
  int j, site, genome, current_state;
  double ran;

  for (site = 0; site < numSites; site++) {
    if (RandUniform(rd) < SNPrate) // there is a SNP at this site
    {
      allSites[site].isSNP = YES;
      if (RandUniform(rd) < 0.5) // in maternal/paternal genome
        genome = MATERNAL;
      else
        genome = PATERNAL;

      if (alphabet == BINARY)
        data[genome][OUTGROUP_ROOT][site] = 1;
      else // DNA mutation according to Mij matrix
      {
        current_state = data[genome][OUTGROUP_ROOT][site];
        ran = RandUniform(rd) * cumMij[current_state][3];
        for (j = 0; j < 4; j++) {
          if (ran <= cumMij[current_state][j]) {
            data[genome][OUTGROUP_ROOT][site] = j;
            break;
          }
        }
      }
    }
  }
}

/********************************** EvolveSitesOnTree
 * ***********************************/
/* Evolves all sites (maternal and paternal genomes) on the given tree
   We assume that a site will be ISM or Mk in both maternal and paternal genome
 */
void EvolveSitesOnTree(TreeNode *treeRoot, int genome, RandomDevice *rd) {
  int i = 0;

  if (rateVarAmongSites == YES)
    for (i = 0; i < numSites; i++)
      allSites[i].rateMultiplier = RandGamma(alphaSites, rd) / alphaSites;

  if (doGeneticSignatures == YES) {
    SimulateSignatureISM(treeRoot, genome,
                         rd); /* we assume diploid ISM for genetic signatures */
  } else {
    if (propAltModelSites == 0) /* only default model (ISM diploid) sites */
    {
      numDefaultModelSites = numSites;
      numAltModelSites = 0;
      for (i = 0; i < numSites; i++)
        DefaultModelSites[i] = i;
      SimulateISM(treeRoot, genome, NO, rd);
    } else if (propAltModelSites == 1.0) /* only alternative model sites */
    {
      numDefaultModelSites = 0;
      numAltModelSites = numSites;
      for (i = 0; i < numSites; i++)
        AltModelSites[i] = i;
      if (altModel == ISMhap)
        SimulateISM(treeRoot, genome, YES, rd);
      else if (altModel == Mk2)
        SimulateMk2(treeRoot, genome, rd);
      else if (altModel == finiteDNA)
        SimulateFiniteDNA(treeRoot, genome, rd);
      else {
        fprintf(stderr,
                "\n\nERROR: Sorry, the specified model is unknown for me");
        exit(0);
      }
    } else /* both ISM and non-ISM sites */
    {
      numDefaultModelSites = 0;
      numAltModelSites = 0;

      for (i = 0; i < numSites; i++) {
        if (RandUniform(rd) < propAltModelSites)
          AltModelSites[numAltModelSites++] = i;
        else
          DefaultModelSites[numDefaultModelSites++] = i;
      }

      /* Evolve ISM sites */
      if (numDefaultModelSites > 0)
        SimulateISM(treeRoot, genome, NO, rd);

      /* Evolve non-SIM sites */
      if (numAltModelSites > 0) {
        if (altModel == ISMhap)
          SimulateISM(treeRoot, genome, YES, rd);
        else if (altModel == Mk2)
          SimulateMk2(treeRoot, genome, rd);
        else if (altModel == finiteDNA)
          SimulateFiniteDNA(treeRoot, genome, rd);
        else {
          fprintf(stderr,
                  "\n\nERROR: Sorry, the specified model is unknown for me");
          exit(-1);
        }
      }
    }
  }
}

/* added code */
void EvolveSitesOnTreeExtd(TreeNode *treeRoot, RandomDevice *rd) {
  if (rateVarAmongSites == YES)
    for (int i = 0; i < numSites; i++)
      allSites[i].rateMultiplier = RandGamma(alphaSites, rd) / alphaSites;

  ThreadedSimulateFiniteMuInDel(
      treeRoot,
      GetGenotypeTotalNum(altModel, constrained)
      );
}

void ThreadedSimulateFiniteMuInDel(TreeNode *treeRoot, int nrOfGenotypes) {
  TreeRootRandomizer treeRootRandomizer[numThreads];
  pthread_attr_t attr;
  double t_start, t_end;

  for (long i = 0; i < numThreads; i++) {
    treeRootRandomizer[i].id = i;
    treeRootRandomizer[i].root = treeRoot;
    treeRootRandomizer[i].numGenotypes = nrOfGenotypes;
  }

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  t_start = omp_get_wtime();
  for (long i = 0; i < numThreads; i++) {
    pthread_create(&threads[i], &attr, ThreadedSimulateFiniteMuInDelForSite,
                   (void *)(treeRootRandomizer + i));
  }
  pthread_attr_destroy(&attr);

  for (long i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }
  t_end = omp_get_wtime();

  fprintf(stderr, "Multi-threading finished in %f seconds.\n", t_end - t_start);

  /* collect the number of mutations, deletions, and insertions for each site
   * it is done here for maximizing the performance of multi-threading
   */
  for (int i = 0; i < numSites; i++) {
    int thIndex = GetThreadIndex(i);
    int relaSite = i - threadStr[thIndex].start;

    allSites[i].numMutations = threadStr[thIndex].summary[relaSite][0];
    allSites[i].numDeletions = threadStr[thIndex].summary[relaSite][3];
    allSites[i].numInsertions = threadStr[thIndex].summary[relaSite][4];

    numMU += allSites[i].numMutations;
    numDEL += allSites[i].numDeletions;
    numINS += allSites[i].numInsertions;
  }
}

void *ThreadedSimulateFiniteMuInDelForSite(void *arg) {
  TreeRootRandomizer *treeRootRandomizer;
  TreeNode *treeRoot;
  long id;
  int nrOfGenotypes;
  double **EPij;
  double t_start, t_end;

  treeRootRandomizer = (TreeRootRandomizer *)arg;
  id = treeRootRandomizer->id;
  treeRoot = treeRootRandomizer->root;
  nrOfGenotypes = treeRootRandomizer->numGenotypes;

  fprintf(stderr, "Thread %ld is starting with %d sites (%d - %d)...\n", id,
          threadStr[id].numSubSites, threadStr[id].start,
          threadStr[id].end - 1);

  t_start = omp_get_wtime();
  for (int i = 0; i < threadStr[id].numSubSites; i++) {
    SimulateFiniteMuInDelForSite(treeRoot, threadStr + id, i, nrOfGenotypes);
  }
  t_end = omp_get_wtime();

  fprintf(stderr, "Thread %ld finished in %f seconds.\n", id, t_end - t_start);
  pthread_exit(NULL);
}

void SimulateFiniteMuInDelForSite(TreeNode *p, ThreadStr *thStr, int site,
                                  int nrOfGenotypes) {
  double branchLength, *EPij = NULL;
  int cell, anccell, ancstate, newstate, nrOfCopies = -1;

  if (p != NULL) {
    if (p->isHealthyRoot == NO) {
      cell = p->label;
      anccell = p->anc->label;
      ancstate = thStr->data[anccell][site].code;

      if (doUserTree == YES)
        branchLength = p->branchLength;
      else {
        if (rateVarAmongSites == YES)
          branchLength = altModelMutationRate * p->length *
                         allSites[site + thStr->start].rateMultiplier;
        else
          branchLength = altModelMutationRate * p->length;
      }

      EPij = FillSubstitutionMatrixExtd(branchLength, nrOfGenotypes);
      newstate = thStr->data[cell][site].code =
          ChooseUniformStateExtd(EPij, ancstate * nrOfGenotypes, thStr->rd);

      free(EPij);
      EPij = NULL;

      if (altModel == FiniteMu) {
        nrOfCopies = 2;
      } else if (altModel == FiniteMuDel) {
        nrOfCopies = GetExistingAllelesNrMD(newstate);
      } else if (altModel == FiniteMuInDel) {
        nrOfCopies = GetExistingAllelesNrMDI(newstate);
      }

      if (nrOfCopies >= 0) {
        thStr->data[cell][site].numExistingAlleles = nrOfCopies;
      } else {
        fprintf(stderr, "ERROR!!! Invalid state (%d)!", newstate);
        exit(-1);
      }

      thStr->data[cell][site].alleles =
          GetValIntToIntArr(genotypeHashMapIntToIntArr, newstate);

      {
        int *tmp1 = thStr->data[anccell][site].pastEvolutionaryEvents;
        int *tmp2 = thStr->data[cell][site].pastEvolutionaryEvents;

        for (int i = 0; i < 5; i++)
          *tmp2++ = *tmp1++;
      }

      if (newstate != ancstate) {
        GetMuDelInsNr(thStr, anccell, cell, site,
                      altModel == FiniteMuInDel ? 3 : 2,
                      allSites[site + thStr->start].referenceAllele);
      }

      if (newstate !=
              WhichRefGenotype(allSites[site + thStr->start].referenceAllele,
                               altModel == FiniteMuInDel ? YES : NO) &&
          p->label < numCells)
        allSites[site + thStr->start].hasEventsOnLeaves = YES;
    }

    SimulateFiniteMuInDelForSite(p->left, thStr, site, nrOfGenotypes);
    SimulateFiniteMuInDelForSite(p->right, thStr, site, nrOfGenotypes);

    if (p->isHealthyRoot == YES) {
      unsigned mutationState[numNodes], numMutatedBranches = 0;
      for (int i = 0; i < numNodes; i++) {
        mutationState[i] =
            thStr->data[i][site].code == thStr->data[p->label][site].code ? 0
                                                                          : 1;
        if (i == healthyTip->label)
          mutationState[i] = 0;
        if (mutationState[i] > 0)
          numMutatedBranches++;
        if (i <= numCells)
          thStr->parallelMutations[site][i] = mutationState[i];
      }

      if (numMutatedBranches > 1) {
        for (int bwd = numNodes - 1; bwd >= 0; bwd--) {
          if (mutationState[bwd] != 0) {
            for (int fwd = 0; fwd < numNodes; fwd++) {
              if (bwd != fwd && mutationState[fwd] != 0 &&
                  treeStructure[bwd][fwd] > 1) {
                mutationState[fwd] = 0;
                numMutatedBranches--;
              }
            }
          }
        }

        if (numMutatedBranches > 1) {
          thStr->hasParallelMutations[site] = YES;
        }
      }
    }
  }
}

double *FillSubstitutionMatrixExtd(double branchLength, int nrOfGenotypes) {
  if (branchLength < 1e-6) {
    double *ch_prob = (double *) malloc((nrOfGenotypes + 1) * (nrOfGenotypes + 1) * sizeof(double));
    for (int i = 0; i < nrOfGenotypes; i++) {
      for (int j = 0; j < nrOfGenotypes; j++) {
        if (i == j)
          ch_prob[i * nrOfGenotypes + j] = 1.0;
        else
          ch_prob[i * nrOfGenotypes + j] = 0.0;
      }
    }
    return ch_prob;
  } else {
    return GetTransitionMtx(
        altModel,
        nrOfGenotypes,
        GetRateMtx(altModel),
        branchLength
        );
  }

  return NULL;
}

int ChooseUniformStateExtd(const double *prob, int start_pos, RandomDevice *rd) {
  int chosenState;
  double ran, cumProb;

  chosenState = 0;
  cumProb = prob[start_pos + chosenState];
  ran = RandUniform(rd);

  while (ran > cumProb)
    cumProb += prob[start_pos + ++chosenState];

  return chosenState;
}
/* end added code */

/*************************** SimulateISM **********************************/
/*	Simulates mutations under an infinite diploid/haploid sites model (ISM).

    Diploid here means that if site 33 mutates in the female genome,
    site 33 in the male genome will not mutate, and viceversa.

    Haploid here means that site 33 in the maternal genome and site 33 in the
   paternal genome will be considered different sites, so both can mutate. The
   rest is as in the diploid, standard ISM

    The number of mutations will be distributed as:

    (1) a Poisson with parameter the sum of
        branch lengths over all sites  Times are already scaled in
        2N. We will force it to get just one mutation per site

        (2) A used-defined number of mutations (segregating sites in the ISM
   diploid)

        Approach 1, which conditions on theta and takes in account
        the branch lengths, is preferred than Hudson's approach 2, which selects
   per each site a random branch and puts a mutation there, and therefore it is
        conditioned on the number of segregating sites. With the coalescent
        you get large or short trees, which imply different number of mutations
        but if you do Hudson's approach you always put the same number of
   mutations.
*/

void SimulateISM(TreeNode *p, int genome, int doISMhaploid, RandomDevice *rd) {
  int i, trials, numMutations, mutationsSoFar;
  double totalBranchSum;
  int *modelSites, numModelSites;

  if (doISMhaploid == NO) {
    modelSites = DefaultModelSites;
    numModelSites = numDefaultModelSites;
  } else {
    modelSites = AltModelSites;
    numModelSites = numAltModelSites;
  }

  /*
      totalBranchSum = 0;
      for (i=0; i<numModelSites; i++)
          {
          allSites[modelSites[i]].branchSum = SumBranches (healthyRoot);
          totalBranchSum += allSites[modelSites[i]].branchSum;
          }
  */

  totalBranchSum = totalTreeLength * numModelSites;

  if (doSimulateFixedNumMutations ==
      YES) /* conditioned on a specific number of mutations */
  {
    if (genome == MATERNAL)
      numMutations = numSNVmaternal = RandBinomial(0.5, numFixedMutations, rd);
    else
      numMutations = numFixedMutations - numSNVmaternal;
  } else
    numMutations = RandPoisson(
        totalBranchSum, rd); /* the number of mutations will be distributed as
                                  a Poisson with parameter totalBranchSum */

  /* if the number of ISM mutations is bigger than the number of ISM sites quit
   * the program and warn the user about violation of ISM */
  if (doISMhaploid == NO && numISMmutations + numMutations > (numModelSites)) {
    fprintf(stderr, "\n\nERROR: The diploid infinite sites model (ISM) has "
                    "been violated. There will be");
    fprintf(stderr,
            "\nmore mutations (%d existing + %d proposed]) than available "
            "sites under this model (%d).",
            numISMmutations, numMutations, numDefaultModelSites);
    fprintf(stderr,
            "\nTry using a smaller mutation rate or effective population size");
    fprintf(stderr,
            "\n(i.e. smaller theta) or check the proportion of JC sites.");
    fprintf(stderr, "\nIf using an user tree, try with smaller branch lengths "
                    "(including the transforming and healthy branches).\n");
    exit(-1);
  } else if (doISMhaploid == YES &&
             numISMmutations + numMutations > (2 * numModelSites)) {
    fprintf(stderr, "\n\nERROR: The haploid infinite sites model (ISM) has "
                    "been violated. There will be");
    fprintf(stderr,
            "\nmore mutations (%d existing + %d proposed]) than available "
            "sites under this model (%d).",
            numISMmutations, numMutations, 2 * numModelSites);
    fprintf(stderr,
            "\nTry using a smaller mutation rate or effective population size");
    fprintf(stderr,
            "\n(i.e. smaller theta) or check the proportion of JC sites.");
    fprintf(stderr, "\nIf using an user tree, try with smaller branch lengths "
                    "(including the transforming and healthy branches).\n");
    exit(-1);
  }

  numISMmutations += numMutations;
  trials = 0;
  mutationsSoFar = 0;
  while (mutationsSoFar < numMutations) {
    i = RandUniformInt(0, numModelSites, rd); /* choose a site at random */
    while (
        ((doISMhaploid == NO) && (allSites[modelSites[i]].numMutations != 0)) ||
        ((doISMhaploid == YES) && (genome == MATERNAL) &&
         (allSites[modelSites[i]].numMutationsMaternal != 0)) ||
        ((doISMhaploid == YES) && (genome == PATERNAL) &&
         (allSites[modelSites[i]].numMutationsPaternal != 0))) {
      i = RandUniformInt(0, numModelSites, rd);
      if (trials++ > 100 * numModelSites) {
        fprintf(stderr,
                "\n\n ERROR: after %d trials cannot find an unmuted site",
                100 * numModelSites);
        fprintf(stderr, "\nmutations = %d   mutations so far = %d\n\n",
                numMutations, mutationsSoFar);
        for (i = 0; i < numModelSites; i++)
          fprintf(stderr, "\nsite %d has %d mutations", modelSites[i] + 1,
                  allSites[modelSites[i]].numMutations);
        exit(-1);
      }
    }
    if (alphabet == DNA)
      SimulateISMDNAforSite(p, genome, modelSites[i], doISMhaploid, rd);
    else
      SimulateISMforSite(p, genome, modelSites[i], doISMhaploid, rd);

    mutationsSoFar++;

#ifdef MYDEBUG
    fprintf(stderr, "\nmutations = %d   mutations so far = %d\n", numMutations,
            mutationsSoFar);
    if (allSites[modelSites[i]].numMutations > 1) {
      fprintf(stderr, "\n\n ERROR: %d mutations in site %d",
              allSites[modelSites[i]].numMutations, modelSites[i] + 1);
      exit(-1);
    }
    for (i = 0; i < numDefaultModelSites; i++)
      fprintf(stderr, "%2d[%d] ", modelSites[i] + 1,
              allSites[modelSites[i]].numMutations);
#endif
  }
}

/********************************** SimulateISMForSite
 * ***********************************/
/*	Simulates a 0/1 mutation under an infinite sites model (ISM) for a given
   site. The branch where this mutation is placed is chosen according to its
   length. 0 is the reference (healthy) allele
 */
void SimulateISMforSite(TreeNode *p, int genome, int site, int doISMhaploid,
                        RandomDevice *rd) {
  static double cumBranchLength, uniform;
  int cell, anccell;

  if (p != NULL) {
    cell = p->label;

    if (p->isHealthyRoot == YES) {
      cumBranchLength = 0;
      uniform = RandUniform(rd) * totalTreeLength;
    } else {
      anccell = p->anc->label;
      cumBranchLength += p->branchLength;

      if ((cumBranchLength < uniform) || /* => there will be no change */
          ((doISMhaploid == NO) && (allSites[site].numMutations > 0)) ||
          ((doISMhaploid == YES) && (genome == MATERNAL) &&
           (allSites[site].numMutationsMaternal > 0)) ||
          ((doISMhaploid == YES) && (genome == PATERNAL) &&
           (allSites[site].numMutationsPaternal > 0))) {
        data[genome][cell][site] = data[genome][anccell][site];
      } else /* => there will be change */
      {
        if (data[genome][anccell][site] ==
            0) /* checking all this might be excessive */
          data[genome][cell][site] = 1;
        else if (data[genome][anccell][site] == 1) {
          if (allSites[site].isSNP == YES) // this is a germline mutation
            data[genome][cell][site] = 0;
          else {
            fprintf(stderr,
                    "\n\nERROR: site %d in genome %d of cell %d cannot mutate "
                    "twice under the ISM model",
                    site, genome, cell);
            exit(-1);
          }
        } else {
          fprintf(stderr,
                  "\n\nERROR: site %d in genome %d of cell %d has an unknow "
                  "state %d under the ISM model",
                  site, genome, anccell, data[genome][anccell][site]);
          exit(-1);
        }

        if (genome == MATERNAL)
          allSites[site].numMutationsMaternal++;
        else if (genome == PATERNAL)
          allSites[site].numMutationsPaternal++;
        allSites[site].numMutations++;
        numMU++;
      }
    }
    SimulateISMforSite(p->left, genome, site, doISMhaploid, rd);
    SimulateISMforSite(p->right, genome, site, doISMhaploid, rd);
  }
}

/********************************** SimulateISMDNAforSite
 * ***********************************/
/*	Simulates a ACGT mutation under an infinite sites model (ISM) for a
   given site. The branch where this mutation is placed is chosen according to
   its length. The reference (healthy) allele for each site will be determined
   by the nucleotide frequencies A mutation matrix will define the probability
   of changing to other nucleotides given the healthy alelle chose
 */
void SimulateISMDNAforSite(TreeNode *p, int genome, int site, int doISMhaploid,
                           RandomDevice *rd) {
  static double cumBranchLength, uniform, ran;
  int j, cell, anccell, ancstate;

  if (p != NULL) {
    cell = p->label;

    if (p->isHealthyRoot == YES) {
      cumBranchLength = 0;
      uniform = RandUniform(rd) * totalTreeLength;
    } else {
      anccell = p->anc->label;
      ancstate = data[genome][anccell][site];
      cumBranchLength += p->branchLength;

      if ((cumBranchLength < uniform) || /* => there will be no change */
          ((doISMhaploid == NO) && (allSites[site].numMutations > 0)) ||
          ((doISMhaploid == YES) && (genome == MATERNAL) &&
           (allSites[site].numMutationsMaternal > 0)) ||
          ((doISMhaploid == YES) && (genome == PATERNAL) &&
           (allSites[site].numMutationsPaternal > 0))) {
        data[genome][cell][site] = ancstate;
      } else /* => there will be change */
      {
        ran = RandUniform(rd) * cumMij[ancstate][3];
        for (j = 0; j < 4; j++) {
          if (ran <= cumMij[ancstate][j]) {
            data[genome][cell][site] = j;
            break;
          }
        }
        if (genome == MATERNAL)
          allSites[site].numMutationsMaternal++;
        else if (genome == PATERNAL)
          allSites[site].numMutationsPaternal++;
        allSites[site].numMutations++;
        numMU++;
      }
    }
    SimulateISMDNAforSite(p->left, genome, site, doISMhaploid, rd);
    SimulateISMDNAforSite(p->right, genome, site, doISMhaploid, rd);
  }
}

/*************************** SimulateSignatureISM
 * **********************************/
/*	Simulates genetic signatures under an infinite diploid trinucleotide
model (ISM).

    Diploid here means that if site 33 mutates in the female genome,
    site 33 in the male genome will not mutate, and viceversa.

    The number of mutations will be distributed as:

    (1) a Poisson with parameter the sum of
        branch lengths over all sites  Times are already scaled in
        2N. We will force it to get just one mutation per site

        (2) A used-defined number of mutations (segregating sites in the ISM
diploid)

=== Importantly, in the current implementation we assume that the context does
not evolve ==

 Also, note that to get the expected counts under each signature it is important
to simulate long genomes, otherwise trinucleotide exhaustion will bias the
counts

*/

void SimulateSignatureISM(TreeNode *p, int genome, RandomDevice *rd) {
  int i, numMutations, mutationsSoFar;
  double totalBranchSum;
  int newState, chosenSite;

  newState = -9;

  /* allocate memory to store the trinucleotide frequencies */
  triNucleotideMaternal = (TriNucStr *)calloc(64, sizeof(TriNucStr));
  if (!triNucleotideMaternal) {
    fprintf(stderr, "Could not allocate the triNucleotideMaternal vector\n");
    exit(-1);
  }
  for (i = 0; i < 64; i++) {
    triNucleotideMaternal[i].tempLength = 100;
    triNucleotideMaternal[i].numAvailablePositions = 0;
    triNucleotideMaternal[i].position = (int *)calloc(100, sizeof(int));
  }

  triNucleotidePaternal = (TriNucStr *)calloc(64, sizeof(TriNucStr));
  if (!triNucleotidePaternal) {
    fprintf(stderr, "Could not allocate the triNucleotidePaternal vector\n");
    exit(-1);
  }
  for (i = 0; i < 64; i++) {
    triNucleotidePaternal[i].tempLength = 100;
    triNucleotidePaternal[i].numAvailablePositions = 0;
    triNucleotidePaternal[i].position = (int *)calloc(100, sizeof(int));
  }

  /* count trinucleotide frequencies in the healthy root genome */
  CountTriNucFrequencies(data[genome][OUTGROUP_ROOT], genome);
  /*
      totalBranchSum = 0;
      for (site=0; site<numSites; site++)
          {
          allSites[site].branchSum = SumBranches (healthyRoot);
          totalBranchSum += allSites[site].branchSum;
          }
  */
  totalBranchSum = totalTreeLength * numSites;

  if (doSimulateFixedNumMutations ==
      YES) /* conditioned on a specified number of mutations */
  {
    if (genome == MATERNAL)
      numMutations = numSNVmaternal = RandBinomial(0.5, numFixedMutations, rd);
    else
      numMutations = numFixedMutations - numSNVmaternal;
  } else
    numMutations = RandPoisson(
        totalBranchSum, rd); /* the number of mutations will be distributed as
                                  a Poisson with parameter totalBranchSum */

  /* if the number of ISM mutations is bigger than the number of ISM sites quit
   * the program and warn the user about violation of ISM */
  if ((numISMmutations + numMutations) > numSites) {
    fprintf(stderr, "\n\nERROR: The diploid infinite sites model (ISM) has "
                    "been violated. There will be");
    fprintf(stderr,
            "\nmore mutations (%d existing + %d proposed]) than available "
            "sites under this model (%d).",
            numISMmutations, numMutations, numDefaultModelSites);
    fprintf(stderr,
            "\nTry using a smaller mutation rate or effective population size");
    fprintf(stderr,
            "\n(i.e. smaller theta) or check the proportion of JC sites.");
    fprintf(stderr, "\nIf using an user tree, try with smaller branch lengths "
                    "(including the transforming and healthy branches).\n");
    exit(-1);
  }

  numISMmutations += numMutations;

#ifdef PRINT_TRIMUTCOUNTER
  if (genome == MATERNAL) {
    triMutationsCounter = (int *)calloc(96, sizeof(int));
    if (!triMutationsCounter) {
      fprintf(stderr, "Could not allocate the triMutationsCounter vector\n");
      exit(-1);
    }
  }
#endif

  mutationsSoFar = 0;
  while (mutationsSoFar < numMutations) {
    chosenSite = ChooseTrinucleotideSite(rd, &newState, genome);

    if (allSites[chosenSite].numMutations == 0) /* force ISM diploid */
    {
      SimulateSignatureISMforSite(p, genome, chosenSite, newState, rd);
      mutationsSoFar++;
    } else // chose another site as this site has already changed in the
           // maternal chromosome
    {
      // fprintf (stderr, "\nsite %d has already changed", chosenSite+1);
    }
  }

  // fprintf (stderr, "\n\nnumISMmutations = %d  numMutations = %d
  // mutationsSoFar = %d (genome = %d)",numISMmutations, numMutations,
  // mutationsSoFar, genome);

#ifdef PRINT_TRIMUTCOUNTER
  int j, k;
  int doDeconstructSigs = YES;
  FILE *fpCounts, *fpRefGenome;

  sprintf(File, "%s/%s", resultsDir, "trimutcounts");
  if ((fpCounts = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  if (doDeconstructSigs == NO) /* counter for SomaticSignatures */
  {
    for (i = 0; i < 96; i++)
      fprintf(fpCounts, "%d\n", triMutationsCounter[i]);
  } else /* counter for deconstructSigs */
  {
    for (j = 0; j < 6; j++)
      for (i = 0; i < 4; i++)
        for (k = 0; k < 4; k++)
          fprintf(fpCounts, "%d\n", triMutationsCounter[trimut(i, j, k)]);
  }
  fclose(fpCounts);

  /* print reference genome as genome context for SomaticSignatures */
  sprintf(File, "%s/%s", resultsDir, "contextGenome");
  if ((fpRefGenome = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }
  fprintf(fpRefGenome, ">1\n");
  for (i = 0; i < numSites; i++)
    fprintf(fpRefGenome, "%c", WhichNuc(allSites[i].referenceAllele));
  fclose(fpRefGenome);
#endif
}

/********************************** CountTriNucFrequencies
 * ***********************************/
/*
         Counts the trinucleotide frequencies in the given genome
*/

void CountTriNucFrequencies(int *genome_array, int genome) {
  int i, n1, n2, n3;
  TriNucStr *triNucleotide;

  if (genome == MATERNAL)
    triNucleotide = triNucleotideMaternal;
  else
    triNucleotide = triNucleotidePaternal;

  for (i = 0; i < numSites - 2; i++) {
    n1 = genome_array[i];
    n2 = genome_array[i + 1];
    n3 = genome_array[i + 2];
    RecordTriNucObservation(&triNucleotide[trinuc(n1, n2, n3)], i);
  }

#ifdef PRINT_TRINUC_GENOME
  if (genome == MATERNAL) {
    fprintf(stderr, "\n\nSimulated trinuc maternal genome: ");
    int rest;
    for (i = 0; i < numSites; i++)
      fprintf(stderr, "%c", WhichNuc(allSites[i].referenceAllele));

    fprintf(stderr, "\n\nObserved trinuc freqs: ");
    for (i = 0; i < 64; i++) {
      n1 = i / 16;
      rest = i % 16;
      n2 = rest / 4;
      n3 = rest % 4;
      fprintf(stderr, "\n%2d %c%c%c [%d] = ", i + 1, WhichNuc(n1), WhichNuc(n2),
              WhichNuc(n3), triNucleotide[i].numAvailablePositions);
      for (int j = 0; j < triNucleotide[i].numAvailablePositions; j++)
        fprintf(stderr, "%d ", triNucleotide[i].position[j] + 1);
    }
    fprintf(stderr, "\n\n");
  }
#endif
}

/********************************** ChooseTrinucleotideSite
 * ***********************************/
/*
Find a valid trinucleotide site to introduce a mutation
The trinucleotide signatures are observed mutations counts, so we should choose
sites according to that. We cannot just choose at random

Given signatures are referred by the pyrimidine of the mutated Watson–Crick base
pair, counts for reverse complements trinucleotide changes are
indistinguishable. For example:

        5-AGG-3 => 5-AAG-3 is indistinguishable from
        3-TCC-5 => 3-TTC-5 so is equivalent to
        5-CCT-3 => 5-CTT-3
*/

int ChooseTrinucleotideSite(RandomDevice *rd, int *newState, int genome) {
  int n1, n2, n3, r1, r2, r3, rest, index, trindex, trials;
  int chosenSignature, chosenTriMutation, chosenPosition;
  int signID;
  TriNucStr *triNucleotide;

  if (genome == MATERNAL)
    triNucleotide = triNucleotideMaternal;
  else
    triNucleotide = triNucleotidePaternal;

  /* 1: Choose a mutation signature */
  signID = ChooseUniformState(signatureWeight, rd);
  chosenSignature = signatureID[signID];

  /* 2: Select a mutation out of the 96 possible ones for the choosen signature
   */
  // chosenTriMutation = ChooseTrinucMutation(geneticSignature[chosenSignature],
  // seed);
  chosenTriMutation =
      ChooseUniformState(signatureProbs[chosenSignature - 1], rd);

  /* 3: Define in which of the 64 trinucleotide classes this mutation will take
   * place */
  r1 = chosenTriMutation / 24;
  rest = chosenTriMutation % 24;
  if (rest / 4 < 3)
    r2 = C;
  else
    r2 = T;
  r3 = rest % 4;

  if (RandUniform(rd) > 0.5) // plus strand
  {
    n1 = r1;
    n2 = r2;
    n3 = r3;
    *newState = targetTriChange[rest / 4];
    // fprintf (stderr, "\n\nplus strand");
  } else // minus strand
  {
    n1 = complementBase[r3];
    n2 = complementBase[r2];
    n3 = complementBase[r1];
    *newState = complementBase[targetTriChange[rest / 4]];
    // fprintf (stderr, "\n\ncomplement reverse");
  }

  trindex = trinuc(n1, n2, n3);

#ifdef PRINT_TRIMUTATIONS
  fprintf(stderr, "\n\nchosen triMutation (id %d) = %c%c%c", chosenTriMutation,
          WhichNuc(r1), WhichNuc(r2), WhichNuc(r3));
  fprintf(stderr, "\n%c%c%c [%d] => ", WhichNuc(n1), WhichNuc(n2), WhichNuc(n3),
          triNucleotide[trindex].numAvailablePositions);
  fprintf(stderr, "%c%c%c", WhichNuc(n1), WhichNuc(*newState), WhichNuc(n3));
  fprintf(stderr, " (%c>%c)", WhichNuc(n2), WhichNuc(*newState));
#endif

  trials = 0;
  while (triNucleotide[trindex].numAvailablePositions == 0) {
    /* 1: Choose a mutation signature */
    signID = ChooseUniformState(signatureWeight, rd);
    chosenSignature = signatureID[signID];
    // fprintf (stderr,"\nchosenSignature = %d", chosenSignature);
    /* 2: Select a mutation out of the 96 possible ones for the choosen
     * signature */
    chosenTriMutation =
        ChooseUniformState(signatureProbs[chosenSignature - 1], rd);

    /* 3: Define the trinucleotide class where this mutation would take place */
    r1 = chosenTriMutation / 24;
    rest = chosenTriMutation % 24;
    if (rest / 4 < 3)
      r2 = C;
    else
      r2 = T;
    r3 = rest % 4;

    if (RandUniform(rd) > 0.5) // plus strand
    {
      n1 = r1;
      n2 = r2;
      n3 = r3;
      *newState = targetTriChange[rest / 4];
      // fprintf (stderr, "\nplus strand");
    } else // minus strand
    {
      n1 = complementBase[r3];
      n2 = complementBase[r2];
      n3 = complementBase[r1];
      *newState = complementBase[targetTriChange[rest / 4]];
      // fprintf (stderr, "\ncomplement reverse");
    }

    trindex = trinuc(n1, n2, n3);

#ifdef PRINT_TRIMUTATIONS
    fprintf(stderr, "\nnew chosen triMutation (id %d) = %c%c%c",
            chosenTriMutation, WhichNuc(r1), WhichNuc(r2), WhichNuc(r3));
    fprintf(stderr, "\n%c%c%c [%d] => ", WhichNuc(n1), WhichNuc(n2),
            WhichNuc(n3), triNucleotide[trindex].numAvailablePositions);
    fprintf(stderr, "%c%c%c", WhichNuc(n1), WhichNuc(*newState), WhichNuc(n3));
    fprintf(stderr, " (%c>%c)", WhichNuc(n2), WhichNuc(*newState));
#endif

    if (trials++ > 100 * numSites) {
      fprintf(stderr,
              "\n\nERROR: after %d trials cannot find an non-mutated trinuc "
              "site, increase the number of sites simulated\n\n",
              100 * numSites);
      exit(-1);
    }
  }

#ifdef PRINT_TRIMUTATIONS
  fprintf(stderr, "\nfinally chosen trimutation index = %d", chosenTriMutation);
#endif

  /* 3: Select the genomic position where this mutation will take place */
  /* note we assume for convenience that the context is constant along the
   * phylogeny */
  if (triNucleotide[trindex].numAvailablePositions == 1)
    index = 0;
  else if (triNucleotide[trindex].numAvailablePositions > 1)
    index = RandUniformInt(0, triNucleotide[trindex].numAvailablePositions, rd);
  else {
    fprintf(stderr, "\n\nWe have a problem, no available trinucleotides\n");
    exit(-1);
  }

  chosenPosition = triNucleotide[trindex].position[index] + 1;
#ifdef PRINT_TRIMUTATIONS
  fprintf(stderr, "\nposition selected = %d", chosenPosition + 1);
#endif

#ifdef PRINT_TRIMUTCOUNTER
  if (allSites[chosenPosition].numMutations == 0)
    triMutationsCounter[chosenTriMutation]++;
#endif

  /* 4: remove the chosen position from the available list */
  RemoveTriNucObservation(&triNucleotide[trindex], index);

  return chosenPosition;
}

/********************************** SimulateSignatureISMforSite
 * ***********************************/
/*	Simulates genetic signatures under an infinite sites model (ISM) for a
   given site. The branch where this mutation is placed is chosen according to
   its length.
 */

void SimulateSignatureISMforSite(TreeNode *p, int genome, int site,
                                 int newState, RandomDevice *rd) {
  static double cumBranchLength, uniform;
  int cell, anccell, ancState;

  if (p != NULL) {
    cell = p->label;

    if (p->isHealthyRoot == YES) {
      cumBranchLength = 0;
      uniform = RandUniform(rd) * totalTreeLength;
      // fprintf (stderr, "\nsimulating at root of position %d", site+1);
    } else {
      anccell = p->anc->label;
      ancState = data[genome][anccell][site];
      cumBranchLength += p->branchLength;

      if (cumBranchLength < uniform ||
          allSites[site].numMutations >
              0) /* => there will be no change; forcing ISM diploid */
      {
        data[genome][cell][site] = ancState;
        // fprintf (stderr, "\ncopying ancstate %c", WhichNuc(ancstate));
      } else /* => there will be change */
      {
        data[genome][cell][site] = newState;
        if (genome == MATERNAL)
          allSites[site].numMutationsMaternal++;
        else if (genome == PATERNAL)
          allSites[site].numMutationsPaternal++;
        allSites[site].numMutations++;
        numMU++;
        // fprintf (stderr, "\nchange at site %d from %c to %c  (numMU = %d)",
        // site+1, WhichNuc(ancState), WhichNuc(newState), numMU);
      }
    }
    SimulateSignatureISMforSite(p->left, genome, site, newState, rd);
    SimulateSignatureISMforSite(p->right, genome, site, newState, rd);
  }
}

/********************************** RemoveTriNucObservation
 * ***********************************/
/* Remove the occurrence of a trinucleotide along the genome */
void RemoveTriNucObservation(TriNucStr *trin, int site) {
  int i;

  for (i = site; i < trin->tempLength - 1; i++)
    trin->position[i] = trin->position[i + 1];

  trin->tempLength = --trin->numAvailablePositions;

  if (trin->numAvailablePositions > 0 &&
      trin->numAvailablePositions % 100 == 0) {
    // fprintf (stderr, "\nreallocating now (available positions =
    // %d)",trin->numAvailablePositions);
    trin->position =
        (int *)realloc(trin->position, trin->tempLength * sizeof(int));
    if (!trin->position) {
      fprintf(
          stderr,
          "Could not reallocate trinucleotide position array when removing\n");
      exit(-1);
    }
  }
  /*	else if (trin->numAvailablePositions == 0)
                 {
                 fprintf (stderr, "\nthis class has zero available positions
     now");
                 }
 */

  /* NOTE:  we could put this nucleotide class signature to zero, as it cannot
     mutate again, always considering the strand
                    But it does not seem necessary 240518 */
}

/************************************* SimulateMk2
 * **********************************************/
/* Simulates the nucleotide substitution process under the Mk2 model (see Lewis
   2001), also called Cavender-Farris-Neyman CFN  model or Jukes-Cantor (1969)
   model for two alleles */

void SimulateMk2(TreeNode *p, int genome, RandomDevice *rd) {
  int i;

  for (i = 0; i < numAltModelSites; i++)
    SimulateMk2forSite(p, genome, AltModelSites[i], rd);
}

/************************************* SimulateMk2ForSite
 * ***************************************/
/* Simulates the nucleotide substitution process for a given site under Mk2
 model (see Lewis 2001) with equal rates. 0 is the reference (healthy) allele */

void SimulateMk2forSite(TreeNode *p, int genome, int site, RandomDevice *rd) {
  double probOfChange, uniform, branchLength;
  int cell, anccell;

  if (p != NULL) {
    if (p->isHealthyRoot == NO) {
      cell = p->label;
      anccell = p->anc->label;

      if (doUserTree == YES)
        branchLength = p->branchLength;
      else {
        if (rateVarAmongSites == YES)
          branchLength =
              altModelMutationRate * p->length * allSites[site].rateMultiplier;
        else
          branchLength = altModelMutationRate * p->length;
      }

      probOfChange = 0.5 - 0.5 * exp(-2.0 * branchLength);

      uniform = RandUniform(rd);
      if (uniform >= probOfChange) /* => no change */
        data[genome][cell][site] = data[genome][anccell][site];
      else /* => there will be change */
      {
        if (data[genome][anccell][site] == 0)
          data[genome][cell][site] = 1;
        else
          data[genome][cell][site] = 0;

        if (genome == MATERNAL)
          allSites[site].numMutationsMaternal++;
        else if (genome == PATERNAL)
          allSites[site].numMutationsPaternal++;
        allSites[site].numMutations++;
        numMU++;
      }
    }
    SimulateMk2forSite(p->left, genome, site, rd);
    SimulateMk2forSite(p->right, genome, site, rd);
  }
}

/************************************* SimulateFiniteDNA
 * **********************************************/
/* Simulates the nucleotide substitution process under a 4-state Markov model
 * including JC, HKY, GTR and GTRnr */
/* Note that beta is set such that mean substitution rate will be 1.
        E.g., for J-C model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */

void SimulateFiniteDNA(TreeNode *p, int genome, RandomDevice *rd) {
  int i, j;

  if (doJC == YES) {
    beta = 4. / 3;
  } else if (doHKY == YES) {
    freqR = freq[A] + freq[G];
    freqY = freq[C] + freq[T];
    freqAG = freq[A] * freq[G];
    freqCT = freq[C] * freq[T];
    kappa = (titv * freqR * freqY) / (freqAG + freqCT);
    beta = 0.5 / (freqR * freqY + kappa * (freqAG + freqCT));
  } else if (doGTR == YES || doGTnR == YES) {
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
        Qij[i * 4 + j] = Mij[i][j] / Mij[2][3] * freq[j];
    mr = 0;
    for (i = 0; i < 4; i++) {
      Qij[i * 4 + i] = 0;
      Qij[i * 4 + i] =
          -(Qij[i * 4] + Qij[i * 4 + 1] + Qij[i * 4 + 2] + Qij[i * 4 + 3]);
      mr -= freq[i] * Qij[i * 4 + i];
    }
    EigenREV(Root, Cijk);
  }

  for (i = 0; i < numAltModelSites; i++)
    SimulateFiniteDNAforSite(p, genome, AltModelSites[i], rd);
}

/************************************* SimulateFiniteDNAforSite
 * **********************************************/
/* Simulates JC, HKY, GTR or GTRnr for a given site */
void SimulateFiniteDNAforSite(TreeNode *p, int genome, int site,
                              RandomDevice *rd) {
  double branchLength, Pij[4][4];
  int cell, anccell, ancstate, newstate;

  if (p != NULL) {
    if (p->isHealthyRoot == NO) {
      cell = p->label;
      anccell = p->anc->label;
      ancstate = data[genome][anccell][site];

      if (doUserTree == YES)
        branchLength = p->branchLength;
      else {
        if (rateVarAmongSites == YES)
          branchLength =
              altModelMutationRate * p->length * allSites[site].rateMultiplier;
        else
          branchLength = altModelMutationRate * p->length;
      }

      FillSubstitutionMatrix(Pij, branchLength);
      newstate = data[genome][cell][site] =
          ChooseUniformState(Pij[ancstate], rd);

      if (newstate != ancstate) {
        if (genome == MATERNAL)
          allSites[site].numMutationsMaternal++;
        else if (genome == PATERNAL)
          allSites[site].numMutationsPaternal++;
        allSites[site].numMutations++;
        numMU++;
      }
    }
    SimulateFiniteDNAforSite(p->left, genome, site, rd);
    SimulateFiniteDNAforSite(p->right, genome, site, rd);
  }
}

/********************* FillSubstitutionMatrix **********************/
/* Sets the apropriate model of nucleotide substitution   */
void FillSubstitutionMatrix(double ch_prob[4][4], double branchLength) {
  int i, j;

  if (branchLength < 1e-6) {
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        if (i == j)
          ch_prob[i][j] = 1.0;
        else
          ch_prob[i][j] = 0.0;
      }
    }
  } else if (doJC == YES)
    JCmodel(ch_prob, branchLength);
  else if (doHKY == YES)
    HKYmodel(ch_prob, branchLength);
  else if (doGTR == YES)
    GTRmodel(ch_prob, branchLength);
  else if (doGTnR == YES)
    GTRmodel(ch_prob, branchLength);
}

/*********************************** JC **************************************/
/*	JC performs Jukes-Cantor 69 correction */
/* Note that beta was set such that mean substitution rate will be 1.
  for the JC model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */
void JCmodel(double Pij[4][4], double branchLength) {
  int i, j;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      if (i == j)
        Pij[i][j] = 0.25 + 0.75 * exp(beta * -branchLength);
      else
        Pij[i][j] = 0.25 - 0.25 * exp(beta * -branchLength);
    }
  }
}

/*********************************** HKY **************************************/
/*	HKY performs Hasegawa-Kishino-Yano 85 correction */

void HKYmodel(double Pij[4][4], double branchLength) {
  int i, j;
  double AA, t, PIj;

  t = branchLength;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      if (j == A || j == G) /* purine */
        PIj = freqR;
      else
        PIj = freqY; /* pyrimidine */

      AA = 1 + PIj * (kappa - 1.0);

      if (i == j)
        Pij[i][j] = freq[j] + freq[j] * (1 / PIj - 1) * exp(-beta * t) +
                    ((PIj - freq[j]) / PIj) * exp(-beta * t * AA);
      else if ((i == A && j == G) || (i == C && j == T) || (i == G && j == A) ||
               (i == T && j == C)) /* transition */
        Pij[i][j] = freq[j] + freq[j] * (1 / PIj - 1) * exp(-beta * t) -
                    (freq[j] / PIj) * exp(-beta * t * AA);
      else /* transversion */
        Pij[i][j] = freq[j] * (1 - exp(-beta * t));
    }
  }
}

/*************** GTR **********************/
void GTRmodel(double Pij[4][4], double branchLength) {
  int i, j, k;
  double t, expt[4];

  t = branchLength;

  /* P(t)ij = SUM Cijk * exp{Root*t} */
  for (k = 1; k < 4; k++)
    expt[k] = exp(t * Root[k]);
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) {
      Pij[i][j] = Cijk[i * 4 * 4 + j * 4 + 0];
      for (k = 1; k < 4; k++)
        Pij[i][j] += Cijk[i * 4 * 4 + j * 4 + k] * expt[k];
    }
}

/*
 * added code
 */

/* end added code */

/********************************** EvolveCNLOHonTree
 * ***********************************/
/*	Simulates point copy-number neutral loss of heterozygosity events
  (CN-LOH) under an infinite haploid sites model (ISM). Haploid here means that
  site 33 in the maternal genome and site 33 in the paternal genome will be
  considered different sites, so both can mutate.

  This is a very simple mode, where CN_LOH evebts are simulated after the
  sequences have been evolved (i.e., SNVs are already in place)

   The number of point CNLOH will be distributed as a Poisson with parameter the
  sum of branch lengths (node length * deletion rate) over all sites.

       We will force it to get at most one CNLOH per site
*/
void EvolveCNLOHonTree(TreeNode *p, int genome, RandomDevice *rd) {
  int i, trials, numCNLOH, CNLOHsoFar;
  double totalCNLOHbranchSum;

  totalCNLOHbranchSum = numSites * totalTreeLength / mutationRate * CNLOHrate;
  numCNLOH =
      RandPoisson(totalCNLOHbranchSum,
                  rd); /* the number of CNLOH events will be distributed as
                            a Poisson with parameter totalCNLOHbranchSum */

  /* if the number of CNLOH events is bigger than the number of sites quit the
   * program and warn the user about violation of ISM */
  if (numISMCNLOH + numCNLOH > (2 * numSites)) {
    fprintf(stderr, "\n\nERROR: The haploid infinite sites model (ISM) for "
                    "CNLOH events has been violated. There will be");
    fprintf(stderr,
            "\nmore CNLOH events (%d existing + %d proposed]) than available "
            "sites under this model (%d).",
            numISMCNLOH, numCNLOH, 2 * numSites);
    fprintf(stderr, "\nTry using a smaller deletion rate.");
    fprintf(
        stderr,
        "\nIf using a user tree, try introducing smaller branch lengths.\n");
    exit(-1);
  }

  numISMCNLOH += numCNLOH;
  trials = 0;
  CNLOHsoFar = 0;
  while (CNLOHsoFar < numCNLOH) {
    i = RandUniformInt(0, numSites, rd); /* choose a site at random */
    while ((genome == MATERNAL && allSites[i].numCNLOHmaternal != 0) ||
           (genome == PATERNAL && allSites[i].numCNLOHpaternal != 0)) {
      i = RandUniformInt(0, numSites, rd);
      if (trials++ > 100 * numSites) {
        fprintf(stderr,
                "\n\n ERROR: after %d trials cannot find a site without CNLOH",
                100 * numSites);
        fprintf(stderr, "\nCNLOH = %d  CNLOH so far = %d\n\n", numCNLOH,
                CNLOHsoFar);
        for (i = 0; i < numSites; i++)
          fprintf(stderr, "\nsite %d has %d CNLOHs", i + 1,
                  allSites[i].numCNLOH);
        fprintf(stderr, "\n");
        exit(-1);
      }
    }

    SimulateCNLOHforSite(p, genome, i, rd);

    CNLOHsoFar++;

#ifdef MYDEBUG
    fprintf(stderr, "\nCNLOH = %d   nCNLOH so far = %d\n", numCNLOH,
            CNLOHsoFar);
            if (allSites[i]].numDeletions > 1)
                {
              fprintf(stderr, "\n\n ERROR: %d CNLOH in site %d",
                      allSites[i].numCNLOH, i + 1);
              exit(-1);
            }
            for (i = 0; i < numSites; i++)
              fprintf(stderr, "%2d[%d] ", i + 1, allSites[i].numCNLOH);
#endif
  }
}

/********************************** SimulateCNLOHforSite
 * ***********************************/
/*	Simulates a CNLOH event infinite sites model (ISM) for a given site. The
   branch where this mutation is placed is chosen according to its length.
  */
void SimulateCNLOHforSite(TreeNode *p, int genome, int site, RandomDevice *rd) {
  static double cumCNLOHbranchLength, uniform;
  int cell, anccell;

  if (p != NULL) {
    cell = p->label;

    if (p->isHealthyRoot == YES) {
      cumCNLOHbranchLength = 0;
      uniform = RandUniform(rd) * totalTreeLength / mutationRate * CNLOHrate;
      ;
    } else {
      anccell = p->anc->label;
      cumCNLOHbranchLength += p->branchLength / mutationRate * CNLOHrate;

      if ((cumCNLOHbranchLength <
           uniform) || /* => there will be no change at this branch */
          ((genome == MATERNAL) && (allSites[site].numCNLOHmaternal > 0)) ||
          ((genome == PATERNAL) && (allSites[site].numCNLOHpaternal > 0))) {
        data[genome][cell][site] = data[genome][anccell][site];
      } else /* => there will be a CN_LOH event */
      {
        if (genome == MATERNAL) {
          data[genome][cell][site] = data[PATERNAL][cell][site];
          allSites[site].numCNLOHmaternal++;
        } else {
          data[genome][cell][site] = data[MATERNAL][cell][site];
          allSites[site].numCNLOHpaternal++;
        }

        allSites[site].numCNLOH++;
        numCNLOH++;
      }
    }
    SimulateCNLOHforSite(p->left, genome, site, rd);
    SimulateCNLOHforSite(p->right, genome, site, rd);
  }
}

/********************************** EvolveDeletionsOnTree
 * ***********************************/
/*	Simulates point deletions under an infinite haploid sites model (ISM).
   Haploid here means that site 33 in the maternal genome and site 33 in the
  paternal genome will be considered different sites, so both can mutate.

  This is a very simple mode, where deletions are simulated after the sequences
  have been evolved (i.e., SNVs are already in place)

   The number of point deletions will be distributed as a Poisson with parameter
  the sum of branch lengths (node length * deletion rate) over all sites.

       We will force it to get at most one deletion per site
*/
void EvolveDeletionsOnTree(TreeNode *p, int genome, RandomDevice *rd) {
  int i, trials, numDeletions, deletionsSoFar;
  double totalDeletionBranchSum;

  /*
  totalDeletionBranchSum = 0;
  for (i=0; i<numSites; i++)
      {
      allSites[i].deletionBranchSum = SumBranches (healthyRoot) / mutationRate *
  deletionRate;; totalDeletionBranchSum += allSites[i].deletionBranchSum;
      }
      */

  totalDeletionBranchSum =
      numSites * totalTreeLength / mutationRate * deletionRate;
  numDeletions =
      RandPoisson(totalDeletionBranchSum,
                  rd); /* the number of deletions will be distributed as a
                            Poisson with parameter totalDeletionBranchSum */

  /* if the number of deletions is bigger than the number of sites quit the
   * program and warn the user about violation of ISM */
  if (numISMdeletions + numDeletions > (2 * numSites)) {
    fprintf(stderr, "\n\nERROR: The haploid infinite sites model (ISM) for "
                    "deletions has been violated. There will be");
    fprintf(stderr,
            "\nmore deletions (%d existing + %d proposed]) than available "
            "sites under this model (%d).",
            numISMdeletions, numDeletions, 2 * numSites);
    fprintf(stderr, "\nTry using a smaller deletion rate.");
    fprintf(
        stderr,
        "\nIf using a user tree, try introducing smaller branch lengths.\n");
    exit(-1);
  }

  numISMdeletions += numDeletions;
  trials = 0;
  deletionsSoFar = 0;
  while (deletionsSoFar < numDeletions) {
    i = RandUniformInt(0, numSites, rd); /* choose a site at random */
    while ((genome == MATERNAL && allSites[i].numDeletionsMaternal != 0) ||
           (genome == PATERNAL && allSites[i].numDeletionsPaternal != 0)) {
      i = RandUniformInt(0, numSites, rd);
      if (trials++ > 100 * numSites) {
        fprintf(stderr,
                "\n\n ERROR: after %d trials cannot find an undeleted site",
                100 * numSites);
        fprintf(stderr, "\ndeletions = %d   deletions so far = %d\n\n",
                numDeletions, deletionsSoFar);
        for (i = 0; i < numSites; i++)
          fprintf(stderr, "\nsite %d has %d deletions", i + 1,
                  allSites[i].numDeletions);
        fprintf(stderr, "\n");
        exit(-1);
      }
    }

    SimulateDeletionforSite(p, genome, i, rd);

    deletionsSoFar++;

#ifdef MYDEBUG
    fprintf(stderr, "\deletions = %d   deletions so far = %d\n", numDeletions,
            deletionsSoFar);
            if (allSites[i]].numDeletions > 1)
                {
              fprintf(stderr, "\n\n ERROR: %d deletions in site %d",
                      allSites[i].numDeletions, i + 1);
              exit(-1);
            }
            for (i = 0; i < numSites; i++)
              fprintf(stderr, "%2d[%d] ", i + 1, allSites[i].numDeletions);
#endif
  }
}

/********************************** SimulateDeletionforSite
 * ***********************************/
/*	Simulates a point deletion infinite sites model (ISM) for a given site.
   The branch where this mutation is placed is chosen according to its length.
  */
void SimulateDeletionforSite(TreeNode *p, int genome, int site,
                             RandomDevice *rd) {
  static double cumDeletionBranchLength, uniform;
  int cell, anccell;

  if (p != NULL) {
    cell = p->label;

    if (p->isHealthyRoot == YES) {
      cumDeletionBranchLength = 0;
      //            uniform = RandomUniform(seed) *
      //            allSites[site].deletionBranchSum;
      uniform = RandUniform(rd) * totalTreeLength / mutationRate * deletionRate;
      ;
    } else {
      anccell = p->anc->label;
      cumDeletionBranchLength += p->branchLength / mutationRate * deletionRate;

      if ((cumDeletionBranchLength <
           uniform) || /* => there will be no change at this branch */
          ((genome == MATERNAL) && (allSites[site].numDeletionsMaternal > 0)) ||
          ((genome == PATERNAL) && (allSites[site].numDeletionsPaternal > 0))) {
        data[genome][cell][site] = data[genome][anccell][site];
      } else /* => there will be a deletion */
      {
        if (data[genome][anccell][site] !=
            DELETION) /* checking all this might be excessive */
          data[genome][cell][site] = DELETION;
        else if (data[genome][anccell][site] == DELETION) {
          fprintf(stderr,
                  "\n\n ERROR: site %d in genome %d of cell %d cannot be "
                  "deleted twice under the ISM model\n",
                  site, genome, cell);
          exit(-1);
        } else {
          fprintf(stderr,
                  "\n\n ERROR: site %d in genome %d of cell %d has an unknow "
                  "state %d under the ISM model\n",
                  site, genome, anccell, data[genome][anccell][site]);
          exit(-1);
        }

        if (genome == MATERNAL)
          allSites[site].numDeletionsMaternal++;
        else if (genome == PATERNAL)
          allSites[site].numDeletionsPaternal++;
        allSites[site].numDeletions++;
        numDEL++;
      }
    }
    SimulateDeletionforSite(p->left, genome, site, rd);
    SimulateDeletionforSite(p->right, genome, site, rd);
  }
}

/********************* AllelicDropout  ************************/
/*
 Remove alleles from single chromosomes at a given ADO rate per genotype

 We assume that ADOrate is the product of the allele dropout rate as follows:
 A = genotype ADO; a = allele error
 A = 2a * (1 - a) + a^2 = 2a - a^2
 a = 1 - sqrt(1-A)
*/

void AllelicDropout(RandomDevice *rd) {
  int i, j;
  double alleleADOrateMean, alleleADOrateMod;
  double **alleleADOrate;

  // calculate ADO rate for a given allele given the mode of ADO
  if (ADOtype == LOCUS_DROPOUT)
    alleleADOrateMean = 1.0 - sqrt(1.0 - ADOrate);
  else if (ADOtype == SINGLE_ADO)
    alleleADOrateMean = ADOrate;

  /* allocate space for the ADO rates */
  alleleADOrate = (double **)calloc(numCells + 2, sizeof(double **));
  if (!alleleADOrate) {
    fprintf(stderr, "Could not allocate the alleleADOrate structure\n");
    exit(-1);
  }
  for (i = 0; i < numCells + 1; i++) {
    alleleADOrate[i] = (double *)calloc(numSites + 1, sizeof(double *));
    if (!alleleADOrate[i]) {
      fprintf(stderr, "Could not allocate the alleleADOrate[i] structure\n");
      exit(-1);
    }
  }

  // sample rate modifiers
  for (i = 0; i < numCells + 1; i++)
    for (j = 0; j < numSites; j++)
      alleleADOrate[i][j] = alleleADOrateMean;

  if (ADOvarAmongCells == YES)
    for (i = 0; i < numCells + 1; i++) {
      alleleADOrateMod = RandGamma(alphaADOcells, rd) / alphaADOcells;
      for (j = 0; j < numSites; j++)
        alleleADOrate[i][j] *= alleleADOrateMod;
    }
  if (ADOvarAmongSites == YES)
    for (j = 0; j < numSites; j++) {
      alleleADOrateMod = RandGamma(alphaADOsites, rd) / alphaADOsites;
      for (i = 0; i < numCells + 1; i++)
        alleleADOrate[i][j] *= alleleADOrateMod;
    }
  /*
          for (i=0; i<numCells+1; i++)
                  for (j=0; j<numSites; j++)
                          fprintf (stderr, "\ncell %d  site %d
     alleleADOrate[i][j] = %6.4f", i+1, j+1, alleleADOrate[i][j]);
  */

  /* added code */
  if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
        altModel == FiniteMu)) {
    for (i = 0; i < numCells + 1; i++) {
      for (j = 0; j < numSites; j++) {
        if (data[MATERNAL][i][j] != DELETION) {
          if (RandUniform(rd) < alleleADOrate[i][j]) {
            data[MATERNAL][i][j] = ADO;
            allSites[j].hasADO = YES;
          }
        }
        if (data[PATERNAL][i][j] != DELETION) {
          if (RandUniform(rd) < alleleADOrate[i][j]) {
            data[PATERNAL][i][j] = ADO;
            allSites[j].hasADO = YES;
          }
        }
      }
    }
  } else {
    int numADOFlag = 0;
    int alleleFlag;
    const int maxAlleleFlag = altModel == FiniteMuInDel ? 3 : 2;
    int droppableIdx[maxAlleleFlag];

    for (j = 0; j < numSites; ++j) {
      int thIndex = GetThreadIndex(j);
      int relaSite = j - threadStr[thIndex].start;
      const int readBlockIndex = GetReadBlockIndex(j);

      for (i = 0; i < numCells + 1; ++i) {
        /* determine the maximum number of alleles to drop out according to
         * ADOtype */
        if (ADOtype == LOCUS_DROPOUT) {
          numADOFlag = threadStr[thIndex].data[i][relaSite].numExistingAlleles;

          alleleFlag = 0;
          for (int k = 0; k < numADOFlag; ++k) {
            while (alleleFlag <= maxAlleleFlag) {
              if ((threadStr[thIndex].data[i][relaSite].alleles[alleleFlag] >=
                   0) &&
                  (cellExtd[i].site[j].observedNumExistingAlleles > 0)) {
                if (RandUniform(rd) < alleleADOrate[i][j]) {
                  if (noisy > 3)
                    fprintf(stderr,
                            "\ncell %d, site %d, allele %d is dropped out.\n",
                            i + 1, j + 1, alleleFlag);

                  threadStr[thIndex].data[i][relaSite].alleles[alleleFlag] =
                          ADOExtd;

                  allSites[j].hasADO = cellExtd[i].site[j].hasADO = YES;

                  cellExtd[i].site[j].numADO++;

                  threadStr[thIndex].data[i][relaSite].numExistingAlleles--;
                  cellExtd[i].site[j].observedNumExistingAlleles--;

                  int tmp = GetValIntArrToInt(
                          genotypeHashMapIntArrToInt,
                          threadStr[thIndex].data[i][relaSite].alleles, 0);
                  if (tmp >= 0) {
                    threadStr[thIndex].data[i][relaSite].code =
                    cellExtd[i].site[j].observedCode = tmp;
                  } else
                    exit(1);
                }
                break;
              }
              alleleFlag++;
            }
          }
        } else if (ADOtype == SINGLE_ADO && cellExtd[i].site[j].observedNumExistingAlleles > 0) {
          alleleFlag = 0;
          for (int k = 0; k < maxAlleleFlag; k++) {
            if (threadStr[thIndex].data[i][relaSite].alleles[k] >= 0) {
              droppableIdx[alleleFlag] = k;
              alleleFlag++;
            }
          }

          if (alleleFlag == 1)
            alleleFlag = droppableIdx[0];
          else
            alleleFlag = droppableIdx[GetAllele(RandUniform(rd), alleleFlag)];

          if (RandUniform(rd) < alleleADOrate[i][j]) {
            if (noisy > 3)
              fprintf(stderr,
                      "\ncell %d, site %d, allele %d is dropped out.\n",
                      i + 1, j + 1, alleleFlag);

            threadStr[thIndex].data[i][relaSite].alleles[alleleFlag] =
                    ADOExtd;

            allSites[j].hasADO = cellExtd[i].site[j].hasADO = YES;

            cellExtd[i].site[j].numADO++;

            threadStr[thIndex].data[i][relaSite].numExistingAlleles--;
            cellExtd[i].site[j].observedNumExistingAlleles--;

            int tmp = GetValIntArrToInt(
                    genotypeHashMapIntArrToInt,
                    threadStr[thIndex].data[i][relaSite].alleles, 0);
            if (tmp >= 0) {
              threadStr[thIndex].data[i][relaSite].code =
              cellExtd[i].site[j].observedCode = tmp;
            } else
              exit(1);
          }
        }

        if (cellExtd[i].site[j].observedNumExistingAlleles == 0)
          cellExtd[i].site[j].maxNumExistingAllelesPerStrand = 0;

        if (cellExtd[i].site[j].hasADO == YES &&
            cellExtd[i].site[j].trueNumExistingAlleles > 2 &&
            (cellExtd[i].site[j].numADO > 1 || RandUniform(rd) > 1.0 / 3))
          cellExtd[i].site[j].maxNumExistingAllelesPerStrand = 1;
      }
    }
  }
  /* end added code */

  /* original code */
  /*
  for (i = 0; i < numCells + 1; i++) {
    for (j = 0; j < numSites; j++) {
      if (data[MATERNAL][i][j] != DELETION) {
        if (RandomUniform(seed) < alleleADOrate[i][j]) {
          data[MATERNAL][i][j] = ADO;
          allSites[j].hasADO = YES;
        }
      }
      if (data[PATERNAL][i][j] != DELETION) {
        if (RandomUniform(seed) < alleleADOrate[i][j]) {
          data[PATERNAL][i][j] = ADO;
          allSites[j].hasADO = YES;
        }
      }
    }
  }
   */

  // free allele ADO rates
  for (i = 0; i < numCells + 1; i++) {
    free(alleleADOrate[i]);
    alleleADOrate[i] = NULL;
  }
  free(alleleADOrate);
  alleleADOrate = NULL;
}

/********************* GenotypeError  ************************/
/*
 Introduce false positive errors directly in the genotypes
 (as opposed to errors just in the reads) (note Eij=0 when i=j)

 We assume that the genotype error is the product of the allele error as
 follows: G = genotype error; a = allele error G = 2a - a^2 a = 1 - sqrt(1-E)
*/

void GenotypeError(RandomDevice *rd) {
  int i, j;
  double alleleError;

  alleleError = 1.0 - sqrt(1.0 - genotypingError);

  if (alphabet == DNA) {
    for (i = 0; i < numCells + 1; i++) {
      for (j = 0; j < numSites; j++) {
        if (data[MATERNAL][i][j] != ADO && data[MATERNAL][i][j] != DELETION &&
            RandUniform(rd) < alleleError) {
          data[MATERNAL][i][j] =
              ChooseUniformState(Eij[data[MATERNAL][i][j]], rd);
          allSites[j].hasGenotypeError = YES;
        }
        if (data[PATERNAL][i][j] != ADO && data[PATERNAL][i][j] != DELETION &&
            RandUniform(rd) < alleleError) {
          data[PATERNAL][i][j] =
              ChooseUniformState(Eij[data[PATERNAL][i][j]], rd);
          allSites[j].hasGenotypeError = YES;
        }
      }
    }
  } else /* for binary data */
  {
    for (i = 0; i < numCells + 1; i++) {
      for (j = 0; j < numSites; j++) {
        if (data[MATERNAL][i][j] != ADO && data[MATERNAL][i][j] != DELETION &&
            RandUniform(rd) < alleleError) {
          if (data[MATERNAL][i][j] == 0)
            data[MATERNAL][i][j] = 1;
          else if (data[MATERNAL][i][j] == 1)
            data[MATERNAL][i][j] = 0;
          allSites[j].hasGenotypeError = YES;
        }
        if (data[PATERNAL][i][j] != ADO && data[PATERNAL][i][j] != DELETION &&
            RandUniform(rd) < alleleError) {
          if (data[PATERNAL][i][j] == 0)
            data[PATERNAL][i][j] = 1;
          else if (data[PATERNAL][i][j] == 1)
            data[PATERNAL][i][j] = 0;
          allSites[j].hasGenotypeError = YES;
        }
      }
    }
  }
}

/************************* CountTrueVariants  ************************/
/* Count number of variants in a given data set, assuming there is no ADO (yet).
   A variant is defined as no-deletion genotype different from the healthy root
   genotype */

int CountTrueVariants() {
  int cell, site;
  int nVariants = 0;

  for (site = 0; site < numSites; site++) {
    for (cell = 0; cell < numCells + 1; cell++) {
      if ((data[MATERNAL][cell][site] != DELETION &&
           data[MATERNAL][cell][site] != data[MATERNAL][OUTGROUP_ROOT][site]) ||
          (data[PATERNAL][cell][site] != DELETION &&
           data[PATERNAL][cell][site] != data[PATERNAL][OUTGROUP_ROOT][site])) {
        allSites[site].isVariant = YES;
        variantSites[nVariants++] = site;
        break;
      }
    }
  }
  return nVariants;
}

/************************* CountAllelesInObservedGenotypes
 * ************************/
/* Identify reference and alternate alleles plus SNVs, in ingroup plus outgroup
   cell genotypes after introducing mutations, cnLOH, deletions, ADO and
   genotype errors A SNV is defined as no-deletion genotype different from the
   reference genotype */

int CountAllelesInObservedGenotypes() {
  int numAltAlleles, cell, site;
  int countA, countC, countG, countT, countADO, countDEL, nSNVs;

  nSNVs = 0;

  for (site = 0; site < numSites; site++) {
    countA = countC = countG = countT = countADO = countDEL = 0;
    for (cell = 0; cell < numCells + 1; cell++) {
      if (data[MATERNAL][cell][site] == A)
        countA++;
      else if (data[MATERNAL][cell][site] == C)
        countC++;
      else if (data[MATERNAL][cell][site] == G)
        countG++;
      else if (data[MATERNAL][cell][site] == T)
        countT++;
      else if (data[MATERNAL][cell][site] == ADO)
        countADO++;
      else if (data[MATERNAL][cell][site] == DELETION)
        countDEL++;

      if (data[PATERNAL][cell][site] == A)
        countA++;
      else if (data[PATERNAL][cell][site] == C)
        countC++;
      else if (data[PATERNAL][cell][site] == G)
        countG++;
      else if (data[PATERNAL][cell][site] == T)
        countT++;
      else if (data[PATERNAL][cell][site] == ADO)
        countADO++;
      else if (data[PATERNAL][cell][site] == DELETION)
        countDEL++;
    }

    allSites[site].countA = countA;
    allSites[site].countC = countC;
    allSites[site].countG = countG;
    allSites[site].countT = countT;
    allSites[site].countACGT = countA + countC + countG + countT;
    allSites[site].countDropped = countADO + countDEL;

    /* count number of alternate alleles, ignoring ADO or DELETION */
    numAltAlleles = 0;
    if (countA > 0 && allSites[site].referenceAllele != A)
      allSites[site].alternateAlleles[numAltAlleles++] = A;
    if (countC > 0 && allSites[site].referenceAllele != C)
      allSites[site].alternateAlleles[numAltAlleles++] = C;
    if (countG > 0 && allSites[site].referenceAllele != G)
      allSites[site].alternateAlleles[numAltAlleles++] = G;
    if (countT > 0 && allSites[site].referenceAllele != T)
      allSites[site].alternateAlleles[numAltAlleles++] = T;
    allSites[site].numAltAlleles = numAltAlleles;

    /* find out whether this site is a SNV */
    if (numAltAlleles > 0) {
      allSites[site].isSNV = YES;
      SNVsites[nSNVs++] = site;
    }
  }
  return nSNVs;
}

/************************* CountAllelesInMLGenotypes  ************************/
/*
 Used when NGS reads are generated and ML genotypes called.
 A SNV is defined as an ML genotype different from the reference genotype
 Note that ML genotypes cannot be ADO or DELETION
 */
int CountAllelesInMLGenotypes() {
  int numAltAllelesInCalls, i, j;
  int countA, countC, countG, countT, countN, nSNVs;

  nSNVs = 0;

  for (j = 0; j < numSites; j++) {
    countA = countC = countG = countT = countN = 0;
    for (i = 0; i < numCells + 1; i++) {
      if (cell[i].site[j].MLmatAllele == A)
        countA++;
      else if (cell[i].site[j].MLmatAllele == C)
        countC++;
      else if (cell[i].site[j].MLmatAllele == G)
        countG++;
      else if (cell[i].site[j].MLmatAllele == T)
        countT++;
      else
        countN++;

      if (cell[i].site[j].MLpatAllele == A)
        countA++;
      else if (cell[i].site[j].MLpatAllele == C)
        countC++;
      else if (cell[i].site[j].MLpatAllele == G)
        countG++;
      else if (cell[i].site[j].MLpatAllele == T)
        countT++;
      else
        countN++;
    }

    /* keep counts on ML genotypes */
    allSites[j].countA = countA;
    allSites[j].countC = countC;
    allSites[j].countG = countG;
    allSites[j].countT = countT;
    allSites[j].countN = countN;
    allSites[j].countACGT = countA + countC + countG + countT;
    allSites[j].countCellswithData = allSites[j].countACGT / 2.0;

    /* count number of alternate alleles in reads that appear in the ML
     * genotypes */
    numAltAllelesInCalls = 0;
    if (countA > 0 && allSites[j].referenceAllele != A)
      numAltAllelesInCalls++;
    if (countC > 0 && allSites[j].referenceAllele != C)
      numAltAllelesInCalls++;
    if (countG > 0 && allSites[j].referenceAllele != G)
      numAltAllelesInCalls++;
    if (countT > 0 && allSites[j].referenceAllele != T)
      numAltAllelesInCalls++;

    /* find out whether this site is a SNV */
    if (numAltAllelesInCalls > 0) {
      allSites[j].isSNV = YES;
      SNVsites[nSNVs++] = j;
    }
  }
  return nSNVs;
}

/* added code */
int CountTrueEvents() {
  int count[4]; // ACGT
  int countN, nSNVs, cellFlag, isDeleted, isInserted;
  int genotypeStructure[4];

  nSNVs = 0;

  for (int j = 0; j < numSites; j++) {
    count[0] = count[1] = count[2] = count[3] = countN = 0;
    allSites[j].countCellswithData = 0;
    isDeleted = isInserted = NO;

    for (int i = 0; i < numCells + 1; i++) {
      cellFlag = 0;
      GetGenotypeStructure(cellExtd[i].site[j].trueAlleles, genotypeStructure);

      if (cellExtd[i].site[j].trueNumExistingAlleles < 2)
        isDeleted = YES;
      else if (cellExtd[i].site[j].trueNumExistingAlleles > 2)
        isInserted = YES;

      if (genotypeStructure[A] > 0) {
        count[0] += genotypeStructure[A];
        cellFlag = 1;
      }

      if (genotypeStructure[C] > 0) {
        count[1] += genotypeStructure[C];
        cellFlag = 1;
      }

      if (genotypeStructure[G] > 0) {
        count[2] += genotypeStructure[G];
        cellFlag = 1;
      }

      if (genotypeStructure[T] > 0) {
        count[3] += genotypeStructure[T];
        cellFlag = 1;
      }

      if ((genotypeStructure[A] == 0) && (genotypeStructure[C] == 0) &&
          (genotypeStructure[G] == 0) && (genotypeStructure[T] == 0))
        countN++;

      if (cellFlag)
        allSites[j].countCellswithData++;
    }

    allSites[j].countA = count[0];
    allSites[j].countC = count[1];
    allSites[j].countG = count[2];
    allSites[j].countT = count[3];
    allSites[j].countN = countN;
    allSites[j].countACGT = SumInt(count, 4);

    int thIndex = GetThreadIndex(j);
    int relaSite = j - threadStr[thIndex].start;

    if (allSites[j].hasEventsOnLeaves == YES) {
      // record the true events
      SNVsites[nSNVs] = j;

      allSites[j].isSNV = threadStr[thIndex].isSNV[relaSite] =
          mutatedSites[nSNVs] =
              (allSites[j].countACGT - count[allSites[j].referenceAllele]) > 0;
      allSites[j].isDeleted = threadStr[thIndex].isDeleted[relaSite] =
          deletedSites[nSNVs] = isDeleted;
      allSites[j].isInserted = threadStr[thIndex].isInserted[relaSite] =
          insertedSites[nSNVs] = isInserted;

      nSNVs++; // "nSNVs" does not only include SNVs, but also sites only
               // contain insertions or deletions
    }
  }

  return nSNVs;
}

/********************* GenerateAnAlternativeNucleotide ********************/
// int GenerateAnAlternativeNucleotide(int referenceAllele, RandomDevice *rd) {
//   double randomValue;
//   int leftMargin;
//
//   while (1) {
//     randomValue = RandUniform(rd) * cumNucleotides[3];
//     for (int i = 0; i < 4; ++i) {
//       if (i == 0) {
//         leftMargin = 0;
//       } else {
//         leftMargin = cumNucleotides[i - 1];
//       }
//       if (randomValue > leftMargin && randomValue <= cumNucleotides[i] &&
//           i != referenceAllele) {
//         return i;
//       }
//     }
//   }
// }

/********************** ProcessCoverageProbAndDispersion *********************/
void ProcessCoverageProb(int replicate) {
  double mean, variance, probability, dispersion;
  int flags[4];

  flags[0] = flags[1] = flags[2] = flags[3] = 0;

  /* (re)allocate the third dimensional space for
   * (probSucTrialsDetailsAll)probSucTrialsDetails */
  for (int i = 0; i < 4; i++) {
    probNBDetails[i] = (double *)calloc(allelicZeroCoverageProportion[i][0] + 1,
                                        sizeof(double));
    if (!probNBDetails[i]) {
      fprintf(stderr, "Could not allocate probNBDetails[%d]\n", i);
      exit(-1);
    }

    if (fixedTree == YES && numDataSets > 1) {
      if (replicate == 0) {
        probNBDetailsAll[i] = (double *)calloc(
            allelicZeroCoverageProportion[i][0] + 1, sizeof(double));
        if (!probNBDetailsAll[i]) {
          fprintf(stderr, "Could not allocate probNBDetailsAll[%d]\n", i);
          exit(-1);
        }
      } else {
        double *ptr = (double *)realloc(
            probNBDetailsAll[i], (allelicZeroCoverageProportion[i][0] +
                                  allelicZeroCoverageProportionAll[i][0] + 1) *
                                     sizeof(double));
        if (!ptr) {
          fprintf(stderr,
                  "Could not allocate probNBDetailsAll[%d] for data set %d\n",
                  i, replicate);
          exit(-1);
        } else {
          probNBDetailsAll[i] = ptr;
        }
      }
    }
  }

  for (int i = 0; i < numCells + 1; ++i) {
    for (int j = 0; j < numSNVs; ++j) {
      mean = (cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles +
              numericalStabilizer) *
             allSites[SNVsites[j]].allelicCov * cellExtd[i].sizeFactor;
      variance =
          mean +
          pow(cellExtd[i].sizeFactor, 2) *
              pow(cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles +
                      numericalStabilizer,
                  2) *
              allSites[SNVsites[j]].allelicRawVar;
      probability = mean / variance;

      probNBDetails
          [cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles]
          [flags[cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles]] =
              probability;

      if (flags[cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles] ==
          0) {
        probNBInfo[cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles]
                  [0] = probNBInfo[cellExtd[i]
                                       .site[SNVsites[j]]
                                       .observedNumExistingAlleles][1] =
                      probability;

        if (fixedTree == YES && numDataSets > 1 && replicate == 0) {
          probNBInfoAll[cellExtd[i]
                            .site[SNVsites[j]]
                            .observedNumExistingAlleles][0] = probNBInfoAll
              [cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles][1] =
                  probability;
        }
      } else {
        if (probability <
            probNBInfo[cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles]
                      [0]) {
          probNBInfo[cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles]
                    [0] = probability;
        }

        if (probability >
            probNBInfo[cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles]
                      [1]) {
          probNBInfo[cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles]
                    [1] = probability;
        }
      }

      flags[cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles]++;
    }
  }

  for (int i = 0; i < 4; ++i) {
    ComputeMeanVarianceDoubleArr(probNBDetails[i],
                                 allelicZeroCoverageProportion[i][0],
                                 &probNBInfo[i][2], &probNBInfo[i][3], 0);

    if (fixedTree == YES && numDataSets > 1) {
      for (int j = 0; j < allelicZeroCoverageProportion[i][0]; ++j) {
        probNBDetailsAll[i][j + allelicZeroCoverageProportionAll[i][0]] =
            probNBDetails[i][j];
      }

      allelicZeroCoverageProportionAll[i][0] +=
          allelicZeroCoverageProportion[i][0];
      allelicZeroCoverageProportionAll[i][1] +=
          allelicZeroCoverageProportion[i][1];

      if (probNBInfo[i][0] < probNBInfoAll[i][0]) {
        probNBInfoAll[i][0] = probNBInfo[i][0];
      }

      if (probNBInfo[i][1] > probNBInfoAll[i][1]) {
        probNBInfoAll[i][1] = probNBInfo[i][1];
      }

      if (replicate == numDataSets - 1) {
        ComputeMeanVarianceDoubleArr(
            probNBDetailsAll[i], allelicZeroCoverageProportionAll[i][0],
            &probNBInfoAll[i][2], &probNBInfoAll[i][3], 0);
      }
    }
  }
}

/********************** PrintSummaryToLog *********************/
void PrintSummaryToLog(FILE *fp, int replicate, int toStderrFlag) {
  if (replicate == 0) {
    fprintf(fp, "\n\n\nEvolutionary events summary\n");
    fprintf(fp, "---------------------------\n");
  }
  fprintf(fp, "\n[Data set %d]", replicate + 1);
  if (!(numDataSets > 1 && fixedTree == YES)) {
    fprintf(fp, "\nNumber of coalescence events   =   %d", numCA);
  }
  fprintf(fp, "\nNumber of SNVs                 =   %d", numSNVs);
  fprintf(fp, "\nNumber of mutational events    =   %d", numMU);
  fprintf(fp, "\nNumber of deletion events      =   %d", numDEL);
  fprintf(fp, "\nNumber of insertion events     =   %d", numINS);

  fprintf(
      fp,
      "\n\nNumber and proportion of zero VAF (among SNV sites in all cells)  "
      "                      =   %ld, %.4f",
      zeroProportionVAF,
      (double)zeroProportionVAF / (numSNVs * (1 + numCells)));
  fprintf(fp,
          "\n Number and proportion of zero VAF with zero coverage           "
          "                        =   %ld, %.4f",
          zeroCoverage, (double)zeroCoverage / (numSNVs * (1 + numCells)));
  fprintf(fp,
          "\n  Number and proportion of complete deletion cases              "
          "                        =   %ld, %.4f",
          completeDeletions,
          (double)completeDeletions / (numSNVs * (1 + numCells)));
  fprintf(fp,
          "\n  Number and proportion of complete ADO cases                   "
          "                        =   %ld, %.4f",
          completeADOs, (double)completeADOs / (numSNVs * (1 + numCells)));
  fprintf(fp,
          "\n  Number and proportion of non-coverage cases                     "
          "                      =   %ld, %.4f",
          noCoverage, (double)noCoverage / (numSNVs * (1 + numCells)));
  fprintf(fp,
          "\n   Required number of successes                                   "
          "                      =   %.4f",
          pow(allelicCovMean, 2) / (4 * allelicRawVarMean));
  fprintf(fp, "\n   With 0 existing alleles:");
  fprintf(fp,
          "\n    Total number of cases for 0 existing alleles                  "
          "                      =   %d",
          allelicZeroCoverageProportion[0][0]);
  fprintf(fp,
          "\n     Probability: min, max, mean, var                             "
          "                      =   %.4f, %.4f, %.4f, %.4f",
          probNBInfo[0][0], probNBInfo[0][1], probNBInfo[0][2],
          probNBInfo[0][3]);
  fprintf(fp,
          "\n   With 1 existing alleles:                                       "
          "                      =   %d, %.4f",
          allelicZeroCoverageProportion[1][1],
          (double)allelicZeroCoverageProportion[1][1] /
              (numSNVs * (1 + numCells)));
  fprintf(fp,
          "\n    Total number of cases for 1 existing alleles; non-coverage "
          "proportion               =   %d; %.4f",
          allelicZeroCoverageProportion[1][0],
          (double)allelicZeroCoverageProportion[1][1] /
              allelicZeroCoverageProportion[1][0]);
  fprintf(fp,
          "\n     Probability: min, max, mean, var                             "
          "                      =   %.4f, %.4f, %.4f, %.4f",
          probNBInfo[1][0], probNBInfo[1][1], probNBInfo[1][2],
          probNBInfo[1][3]);
  fprintf(fp,
          "\n   With 2 existing alleles:                                       "
          "                      =   %d, %.4f",
          allelicZeroCoverageProportion[2][1],
          (double)allelicZeroCoverageProportion[2][1] /
              (numSNVs * (1 + numCells)));
  fprintf(fp,
          "\n    Total number of cases for 2 existing alleles; non-coverage "
          "proportion               =   %d; %.4f",
          allelicZeroCoverageProportion[2][0],
          (double)allelicZeroCoverageProportion[2][1] /
              allelicZeroCoverageProportion[2][0]);
  fprintf(fp,
          "\n     Probability: min, max, mean, var                             "
          "                      =   %.4f, %.4f, %.4f, %.4f",
          probNBInfo[2][0], probNBInfo[2][1], probNBInfo[2][2],
          probNBInfo[2][3]);
  fprintf(fp,
          "\n   With 3 existing alleles:                                       "
          "                      =   %d, %.4f",
          allelicZeroCoverageProportion[3][1],
          (double)allelicZeroCoverageProportion[3][1] /
              (numSNVs * (1 + numCells)));
  fprintf(fp,
          "\n    Total number of cases for 3 existing alleles; non-coverage "
          "proportion               =   %d; %.4f",
          allelicZeroCoverageProportion[3][0],
          (double)allelicZeroCoverageProportion[3][1] /
              allelicZeroCoverageProportion[3][0]);
  fprintf(fp,
          "\n     Probability: min, max, mean, var                             "
          "                      =   %.4f, %.4f, %.4f, %.4f",
          probNBInfo[3][0], probNBInfo[3][1], probNBInfo[3][2],
          probNBInfo[3][3]);

  fprintf(fp, "\n");

  if (numDataSets > 1 && fixedTree == YES && replicate == numDataSets - 1) {
    fprintf(fp, "\nSummary for all data sets with a fixed coalescence tree\n");
    fprintf(fp, "--------------------------------------------------------");
    fprintf(fp, "\nNumber of coalescence events   =   %d", (int)cumNumCA);
    fprintf(fp, "\nNumber of SNVs                 =   %d", (int)cumNumSNVs);
    fprintf(fp, "\nNumber of mutational events    =   %d", (int)cumNumMU);
    fprintf(fp, "\nNumber of deletion events      =   %d", (int)cumNumDEL);
    fprintf(fp, "\nNumber of insertion events     =   %d", (int)cumNumINS);

    fprintf(
        fp,
        "\n\nNumber and proportion of zero VAF (among SNV sites in all cells)  "
        "                      =   %ld, %.4f",
        zeroProportionVAFAll,
        (double)zeroProportionVAFAll / (cumNumSNVs * (1 + numCells)));
    fprintf(fp,
            "\n Number and proportion of zero VAF with zero coverage           "
            "                        =   %ld, %.4f",
            zeroCoverageAll,
            (double)zeroCoverageAll / (cumNumSNVs * (1 + numCells)));
    fprintf(fp,
            "\n  Number and proportion of complete deletion cases              "
            "                        =   %ld, %.4f",
            completeDeletionsAll,
            (double)completeDeletionsAll / (cumNumSNVs * (1 + numCells)));
    fprintf(fp,
            "\n  Number and proportion of complete ADO cases                   "
            "                        =   %ld, %.4f",
            completeADOsAll,
            (double)completeADOsAll / (cumNumSNVs * (1 + numCells)));
    fprintf(fp,
            "\n  Number and proportion of non-coverage cases                   "
            "                        =   %ld, %.4f",
            noCoverageAll,
            (double)noCoverageAll / (cumNumSNVs * (1 + numCells)));
    fprintf(
        fp,
        "\n   Required number of successes                                   "
        "                      =   %.4f",
        pow(allelicCovMean, 2) / (4 * allelicRawVarMean));
    fprintf(fp, "\n   With 0 existing alleles:");
    fprintf(
        fp,
        "\n    Total number of cases for 0 existing alleles                  "
        "                      =   %d",
        allelicZeroCoverageProportionAll[0][0]);
    fprintf(
        fp,
        "\n     Probability: min, max, mean, var                             "
        "                      =   %.4f, %.4f, %.4f, %.4f",
        probNBInfoAll[0][0], probNBInfoAll[0][1], probNBInfoAll[0][2],
        probNBInfoAll[0][3]);
    fprintf(fp,
            "\n   With 1 existing alleles:                                     "
            "                        =   %d, %.4f",
            allelicZeroCoverageProportionAll[1][1],
            (double)allelicZeroCoverageProportionAll[1][1] /
                (cumNumSNVs * (1 + numCells)));
    fprintf(fp,
            "\n    Total number of cases for 1 existing alleles; non-coverage "
            "proportion               =   %d; %.4f",
            allelicZeroCoverageProportionAll[1][0],
            (double)allelicZeroCoverageProportionAll[1][1] /
                allelicZeroCoverageProportionAll[1][0]);
    fprintf(
        fp,
        "\n     Probability: min, max, mean, var                             "
        "                      =   %.4f, %.4f, %.4f, %.4f",
        probNBInfoAll[1][0], probNBInfoAll[1][1], probNBInfoAll[1][2],
        probNBInfoAll[1][3]);
    fprintf(fp,
            "\n   With 2 existing alleles:                                     "
            "                        =   %d, %.4f",
            allelicZeroCoverageProportionAll[2][1],
            (double)allelicZeroCoverageProportionAll[2][1] /
                (cumNumSNVs * (1 + numCells)));
    fprintf(fp,
            "\n    Total number of cases for 2 existing alleles; non-coverage "
            "proportion               =   %d; %.4f",
            allelicZeroCoverageProportionAll[2][0],
            (double)allelicZeroCoverageProportionAll[2][1] /
                allelicZeroCoverageProportionAll[2][0]);
    fprintf(
        fp,
        "\n     Probability: min, max, mean, var                             "
        "                      =   %.4f, %.4f, %.4f, %.4f",
        probNBInfoAll[2][0], probNBInfoAll[2][1], probNBInfoAll[2][2],
        probNBInfoAll[2][3]);
    fprintf(fp,
            "\n   With 3 existing alleles:                                     "
            "                        =   %d, %.4f",
            allelicZeroCoverageProportionAll[3][1],
            (double)allelicZeroCoverageProportionAll[3][1] /
                (cumNumSNVs * (1 + numCells)));
    fprintf(fp,
            "\n    Total number of cases for 3 existing alleles; non-coverage "
            "proportion               =   %d; %.4f",
            allelicZeroCoverageProportionAll[3][0],
            (double)allelicZeroCoverageProportionAll[3][1] /
                allelicZeroCoverageProportionAll[3][0]);
    fprintf(
        fp,
        "\n     Probability: min, max, mean, var                             "
        "                      =   %.4f, %.4f, %.4f, %.4f",
        probNBInfoAll[3][0], probNBInfoAll[3][1], probNBInfoAll[3][2],
        probNBInfoAll[3][3]);

    fprintf(fp, "\n");
  }

  if (toStderrFlag) {
    if (replicate == 0) {
      fprintf(stderr, "\nEvolutionary events summary\n");
      fprintf(stderr, "---------------------------\n");
    }
    fprintf(stderr, "\n[Data set %d]", replicate + 1);
    if (!(numDataSets > 1 && fixedTree == YES)) {
      fprintf(stderr, "\nNumber of coalescence events   =   %d", numCA);
    }
    fprintf(stderr, "\nNumber of SNVs                 =   %d", numSNVs);
    fprintf(stderr, "\nNumber of mutational events    =   %d", numMU);
    fprintf(stderr, "\nNumber of deletion events      =   %d", numDEL);
    fprintf(stderr, "\nNumber of insertion events     =   %d", numINS);

    fprintf(
        stderr,
        "\n\nNumber and proportion of zero VAF (among SNV sites in all cells)  "
        "                      =   %ld, %.4f",
        zeroProportionVAF,
        (double)zeroProportionVAF / (numSNVs * (1 + numCells)));
    fprintf(stderr,
            "\n Number and proportion of zero VAF with zero coverage           "
            "                        =   %ld, %.4f",
            zeroCoverage, (double)zeroCoverage / (numSNVs * (1 + numCells)));
    fprintf(stderr,
            "\n  Number and proportion of complete deletion cases              "
            "                        =   %ld, %.4f",
            completeDeletions,
            (double)completeDeletions / (numSNVs * (1 + numCells)));
    fprintf(stderr,
            "\n  Number and proportion of complete ADO cases                   "
            "                        =   %ld, %.4f",
            completeADOs, (double)completeADOs / (numSNVs * (1 + numCells)));
    fprintf(stderr,
            "\n  Number and proportion of non-coverage cases                   "
            "                        =   %ld, %.4f",
            noCoverage, (double)noCoverage / (numSNVs * (1 + numCells)));
    fprintf(
        stderr,
        "\n   Required number of successes                                   "
        "                      =   %.4f",
        pow(allelicCovMean, 2) / (4 * allelicRawVarMean));
    fprintf(stderr, "\n   With 0 existing alleles:");
    fprintf(
        stderr,
        "\n    Total number of cases for 0 existing alleles                  "
        "                      =   %d",
        allelicZeroCoverageProportion[0][0]);
    fprintf(
        stderr,
        "\n     Probability: min, max, mean, var                             "
        "                      =   %.4f, %.4f, %.4f, %.4f",
        probNBInfo[0][0], probNBInfo[0][1], probNBInfo[0][2], probNBInfo[0][3]);
    fprintf(stderr,
            "\n   With 1 existing alleles:                                     "
            "                        =   %d, %.4f",
            allelicZeroCoverageProportion[1][1],
            (double)allelicZeroCoverageProportion[1][1] /
                (numSNVs * (1 + numCells)));
    fprintf(stderr,
            "\n    Total number of cases for 1 existing alleles; non-coverage "
            "proportion               =   %d; %.4f",
            allelicZeroCoverageProportion[1][0],
            (double)allelicZeroCoverageProportion[1][1] /
                allelicZeroCoverageProportion[1][0]);
    fprintf(
        stderr,
        "\n     Probability: min, max, mean, var                             "
        "                      =   %.4f, %.4f, %.4f, %.4f",
        probNBInfo[1][0], probNBInfo[1][1], probNBInfo[1][2], probNBInfo[1][3]);
    fprintf(stderr,
            "\n   With 2 existing alleles:                                     "
            "                        =   %d, %.4f",
            allelicZeroCoverageProportion[2][1],
            (double)allelicZeroCoverageProportion[2][1] /
                (numSNVs * (1 + numCells)));
    fprintf(stderr,
            "\n    Total number of cases for 2 existing alleles; non-coverage "
            "proportion               =   %d; %.4f",
            allelicZeroCoverageProportion[2][0],
            (double)allelicZeroCoverageProportion[2][1] /
                allelicZeroCoverageProportion[2][0]);
    fprintf(
        stderr,
        "\n     Probability: min, max, mean, var                             "
        "                      =   %.4f, %.4f, %.4f, %.4f",
        probNBInfo[2][0], probNBInfo[2][1], probNBInfo[2][2], probNBInfo[2][3]);
    fprintf(stderr,
            "\n   With 3 existing alleles:                                     "
            "                        =   %d, %.4f",
            allelicZeroCoverageProportion[3][1],
            (double)allelicZeroCoverageProportion[3][1] /
                (numSNVs * (1 + numCells)));
    fprintf(stderr,
            "\n    Total number of cases for 3 existing alleles; non-coverage "
            "proportion               =   %d; %.4f",
            allelicZeroCoverageProportion[3][0],
            (double)allelicZeroCoverageProportion[3][1] /
                allelicZeroCoverageProportion[3][0]);
    fprintf(
        stderr,
        "\n     Probability: min, max, mean, var                             "
        "                      =   %.4f, %.4f, %.4f, %.4f",
        probNBInfo[3][0], probNBInfo[3][1], probNBInfo[3][2], probNBInfo[3][3]);

    fprintf(stderr, "\n");

    if (numDataSets > 1 && fixedTree == YES && replicate == numDataSets - 1) {
      fprintf(stderr,
              "\nSummary for all data sets with a fixed coalescence tree\n");
      fprintf(stderr,
              "--------------------------------------------------------");
      fprintf(stderr, "\nNumber of coalescence events   =   %d", (int)cumNumCA);
      fprintf(stderr, "\nNumber of SNVs                 =   %d",
              (int)cumNumSNVs);
      fprintf(stderr, "\nNumber of mutational events    =   %d", (int)cumNumMU);
      fprintf(stderr, "\nNumber of deletion events      =   %d",
              (int)cumNumDEL);
      fprintf(stderr, "\nNumber of insertion events     =   %d",
              (int)cumNumINS);

      fprintf(stderr,
              "\n\nNumber and proportion of zero VAF (among SNV sites in all "
              "cells)  "
              "                      =   %ld, %.4f",
              zeroProportionVAFAll,
              (double)zeroProportionVAFAll / (cumNumSNVs * (1 + numCells)));
      fprintf(
          stderr,
          "\n Number and proportion of zero VAF with zero coverage           "
          "                        =   %ld, %.4f",
          zeroCoverageAll,
          (double)zeroCoverageAll / (cumNumSNVs * (1 + numCells)));
      fprintf(
          stderr,
          "\n  Number and proportion of complete deletion cases              "
          "                        =   %ld, %.4f",
          completeDeletionsAll,
          (double)completeDeletionsAll / (cumNumSNVs * (1 + numCells)));
      fprintf(
          stderr,
          "\n  Number and proportion of complete ADO cases                   "
          "                        =   %ld, %.4f",
          completeADOsAll,
          (double)completeADOsAll / (cumNumSNVs * (1 + numCells)));
      fprintf(
          stderr,
          "\n  Number and proportion of non-coverage cases                   "
          "                        =   %ld, %.4f",
          noCoverageAll, (double)noCoverageAll / (cumNumSNVs * (1 + numCells)));
      fprintf(
          stderr,
          "\n   Required number of successes                                   "
          "                      =   %.4f",
          pow(allelicCovMean, 2) / (4 * allelicRawVarMean));
      fprintf(stderr, "\n   With 0 existing alleles:");
      fprintf(
          stderr,
          "\n    Total number of cases for 0 existing alleles                  "
          "                      =   %d",
          allelicZeroCoverageProportionAll[0][0]);
      fprintf(
          stderr,
          "\n     Probability: min, max, mean, var                             "
          "                      =   %.4f, %.4f, %.4f, %.4f",
          probNBInfoAll[0][0], probNBInfoAll[0][1], probNBInfoAll[0][2],
          probNBInfoAll[0][3]);
      fprintf(stderr,
              "\n   With 1 existing alleles:                                   "
              "                          =   %d, %.4f",
              allelicZeroCoverageProportionAll[1][1],
              (double)allelicZeroCoverageProportionAll[1][1] /
                  (cumNumSNVs * (1 + numCells)));
      fprintf(stderr,
              "\n    Total number of cases for 1 existing alleles; "
              "non-coverage proportion               =   %d; %.4f",
              allelicZeroCoverageProportionAll[1][0],
              (double)allelicZeroCoverageProportionAll[1][1] /
                  allelicZeroCoverageProportionAll[1][0]);
      fprintf(
          stderr,
          "\n     Probability: min, max, mean, var                             "
          "                      =   %.4f, %.4f, %.4f, %.4f",
          probNBInfoAll[1][0], probNBInfoAll[1][1], probNBInfoAll[1][2],
          probNBInfoAll[1][3]);
      fprintf(stderr,
              "\n   With 2 existing alleles:                                   "
              "                          =   %d, %.4f",
              allelicZeroCoverageProportionAll[2][1],
              (double)allelicZeroCoverageProportionAll[2][1] /
                  (cumNumSNVs * (1 + numCells)));
      fprintf(stderr,
              "\n    Total number of cases for 2 existing alleles; "
              "non-coverage proportion               =   %d; %.4f",
              allelicZeroCoverageProportionAll[2][0],
              (double)allelicZeroCoverageProportionAll[2][1] /
                  allelicZeroCoverageProportionAll[2][0]);
      fprintf(
          stderr,
          "\n     Probability: min, max, mean, var                             "
          "                      =   %.4f, %.4f, %.4f, %.4f",
          probNBInfoAll[2][0], probNBInfoAll[2][1], probNBInfoAll[2][2],
          probNBInfoAll[2][3]);
      fprintf(stderr,
              "\n   With 3 existing alleles:                                   "
              "                          =   %d, %.4f",
              allelicZeroCoverageProportionAll[3][1],
              (double)allelicZeroCoverageProportionAll[3][1] /
                  (cumNumSNVs * (1 + numCells)));
      fprintf(stderr,
              "\n    Total number of cases for 3 existing alleles; "
              "non-coverage proportion               =   %d; %.4f",
              allelicZeroCoverageProportionAll[3][0],
              (double)allelicZeroCoverageProportionAll[3][1] /
                  allelicZeroCoverageProportionAll[3][0]);
      fprintf(
          stderr,
          "\n     Probability: min, max, mean, var                             "
          "                      =   %.4f, %.4f, %.4f, %.4f",
          probNBInfoAll[3][0], probNBInfoAll[3][1], probNBInfoAll[3][2],
          probNBInfoAll[3][3]);

      fprintf(stderr, "\n");
    }
  }
}
/* end added code */

/*************************** SumBranchlengths
 * **********************************/
/* Returns the sum of the branch lengths for a given tree/subtree */
double SumBranches(TreeNode *p) {
  static double sum;

  if (p != NULL) {
    if (p->anc == NULL)
      sum = 0;
    else
      sum += p->branchLength;
    SumBranches(p->left);
    SumBranches(p->right);
  }

  return sum;
}

/*************************** SumTreeHeights **********************************/
/* Returns the total height (sum of all paths from root to tips)
        and the number of tips for for a given tree/subtree
*/
/*
 double SumTreeHeights (TreeNode *p, double total, int *numTips)
        {
        if (p == NULL)
                return 0;

        if (p->left == NULL)
                {
                (*numTips)++;
                return total;
                }
        else return SumTreeHeights(p->left, total + p->left->branchLength,
 numTips) + SumTreeHeights(p->right, total + p->right->branchLength, numTips);
        }
*/

/********************* CompareGenotypes  ************************/
/*
 Compares two unphased genotypes (i.e, A/T = T/A)
*/
int CompareGenotypes(int a1, int a2, int b1, int b2) {
  int temp, equal;

  temp = 0;
  equal = YES;

  /* order the alleles to facilitate the comparisons */
  if (a1 > a2) {
    temp = a1;
    a1 = a2;
    a2 = temp;
  }

  if (b1 > b2) {
    temp = b1;
    b1 = b2;
    b2 = temp;
  }

  if (a1 != b1 || a2 != b2)
    equal = NO;

  return equal;
}

/********************* AllocateCellStructure  ************************/
/*
 This cell structure keeps information basically about read count functions. It
 was build initially for making doublets.
 */
void AllocateCellStructure() {
  int i, j, k;

  cell = (CellStr *)calloc(numCells + 1, sizeof(CellStr));
  if (!cell) {
    fprintf(stderr, "Could not allocate cell (%ld)\n",
            (numCells + 1) * sizeof(CellStr));
    exit(-1);
  }

  for (i = 0; i < numCells + 1; i++) {
    cell[i].site = (CellSiteStr *)calloc(numSites, sizeof(CellSiteStr));
    if (!cell[i].site) {
      fprintf(stderr, "Could not allocate the cell[%d].site structure\n", i);
      exit(-1);
    }

    if (doNGS == YES) {
      for (j = 0; j < numSites; j++) {
        cell[i].site[j].readCount = (int *)calloc(4, sizeof(int));
        if (!cell[i].site[j].readCount) {
          fprintf(
              stderr,
              "Could not allocate the cell[%d].site[%d].readCount structure\n",
              i, j);
          exit(-1);
        }

        cell[i].site[j].readCountDoublet = (int *)calloc(4, sizeof(int));
        if (!cell[i].site[j].readCountDoublet) {
          fprintf(stderr,
                  "Could not allocate the cell[%d].site[%d].readCountDoublet "
                  "structure\n",
                  i, j);
          exit(-1);
        }

        cell[i].site[j].genLike = (double **)calloc(4, sizeof(double *));
        if (!cell[i].site[j].genLike) {
          fprintf(
              stderr,
              "Could not allocate the cell[%d].site[%d].genLike structure\n", i,
              j);
          exit(-1);
        }
        for (k = 0; k < 4; k++) {
          cell[i].site[j].genLike[k] = (double *)calloc(4, sizeof(double));
          if (!cell[i].site[j].genLike[k]) {
            fprintf(stderr,
                    "Could not allocate the cell[%d].site[%d].genLike[%d] "
                    "structure\n",
                    i, j, k);
            exit(-1);
          }
        }

        cell[i].site[j].scaledGenLike = (double **)calloc(4, sizeof(double *));
        if (!cell[i].site[j].scaledGenLike) {
          fprintf(stderr,
                  "Could not allocate the cell[%d].site[%d].scaledGenLike "
                  "structure\n",
                  i, j);
          exit(-1);
        }
        for (k = 0; k < 4; k++) {
          cell[i].site[j].scaledGenLike[k] =
              (double *)calloc(4, sizeof(double));
          if (!cell[i].site[j].scaledGenLike[k]) {
            fprintf(stderr,
                    "Could not allocate the "
                    "cell[%d].site[%d].scaledGenLike[%d] structure\n",
                    i, j, k);
            exit(-1);
          }
        }

        cell[i].site[j].genLikeDoublet = (double **)calloc(4, sizeof(double *));
        if (!cell[i].site[j].genLikeDoublet) {
          fprintf(stderr,
                  "Could not allocate the cell[%d].site[%d].genLikeDoublet "
                  "structure\n",
                  i, j);
          exit(-1);
        }
        for (k = 0; k < 4; k++) {
          cell[i].site[j].genLikeDoublet[k] =
              (double *)calloc(4, sizeof(double));
          if (!cell[i].site[j].genLikeDoublet[k]) {
            fprintf(stderr,
                    "Could not allocate the "
                    "cell[%d].site[%d].genLikeDoublet[%d] structure\n",
                    i, j, k);
            exit(-1);
          }
        }

        cell[i].site[j].scaledGenLikeDoublet =
            (double **)calloc(4, sizeof(double *));
        if (!cell[i].site[j].scaledGenLikeDoublet) {
          fprintf(stderr,
                  "Could not allocate the "
                  "cell[%d].site[%d].scaledGenLikeDoublet structure\n",
                  i, j);
          exit(-1);
        }
        for (k = 0; k < 4; k++) {
          cell[i].site[j].scaledGenLikeDoublet[k] =
              (double *)calloc(4, sizeof(double));
          if (!cell[i].site[j].scaledGenLikeDoublet[k]) {
            fprintf(stderr,
                    "Could not allocate the "
                    "cell[%d].site[%d].scaledGenLikeDoublet[%d] structure\n",
                    i, j, k);
            exit(-1);
          }
        }
      }
    }
  }
}

/********************* GenerateReadCounts  ************************/
/*	For each individual SNV genotype the program will generate read counts
        given some sequencing depth or coverage. The number of reads follow a
   Poisson or a Negative Binomial distribution around mean coverage. Reads are
   randomly assigned to maternal/paternal chromosome, according to an allelic
   haploid coverage reduction (by default 0.5). Reads can contain some errors
   according to sequencing and amplification error parameters
*/

void GenerateReadCounts(RandomDevice *rd) {
  int i, j, equalGenotypes, countCalledGenotypes;
  int countReadsA, countReadsC, countReadsG, countReadsT;
  int numAltAlleles;
  int trueMaternalAllele, truePaternalAllele, MLmatAllele, MLpatAllele;
  double *probs, **ngsEij, **ampEijmat, **ampEijpat;
  CellSiteStr *c;

  countReadsA = countReadsC = countReadsG = countReadsT = 0;

  probs = (double *)calloc(4, sizeof(double));
  if (!probs) {
    fprintf(stderr, "Could not allocate the probs structure\n");
    exit(-1);
  }

  /* allocate error probabilities for NGS */
  ngsEij = (double **)calloc(4, sizeof(double **));
  if (!ngsEij) {
    fprintf(stderr, "Could not allocate the ngsEij structure\n");
    exit(-1);
  }
  for (i = 0; i < 4; i++) {
    ngsEij[i] = (double *)calloc(4, sizeof(double *));
    if (!ngsEij[i]) {
      fprintf(stderr, "Could not allocate the ngsEij[i] structure\n");
      exit(-1);
    }
  }

  /* allocate error probabilities for amplification */
  ampEijmat = (double **)calloc(4, sizeof(double **));
  if (!ampEijmat) {
    fprintf(stderr, "Could not allocate the ampEijmat structure\n");
    exit(-1);
  }
  for (i = 0; i < 4; i++) {
    ampEijmat[i] = (double *)calloc(4, sizeof(double *));
    if (!ampEijmat[i]) {
      fprintf(stderr, "Could not allocate the ampEijmat[i] structure\n");
      exit(-1);
    }
  }

  ampEijpat = (double **)calloc(4, sizeof(double **));
  if (!ampEijpat) {
    fprintf(stderr, "Could not allocate the ampEijpat structure\n");
    exit(-1);
  }
  for (i = 0; i < 4; i++) {
    ampEijpat[i] = (double *)calloc(4, sizeof(double *));
    if (!ampEijpat[i]) {
      fprintf(stderr, "Could not allocate the ampEijpat[i] structure\n");
      exit(-1);
    }
  }

  /* initialize ngsEij now as this will be contant across sites */
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) {
      if (i == j)
        ngsEij[i][j] = 1.0 - sequencingError;
      else
        ngsEij[i][j] = Eij[i][j] / cumEij[i][3] * sequencingError;
    }

  /* produce read counts for all sites */
  for (j = 0; j < numSites; j++) {
    for (i = 0; i < numCells + 1; i++) {
      SiteReadCounts(&cell[i].site[j], i, j, probs, ngsEij, ampEijmat,
                     ampEijpat, rd);
      if (cell[i].site[j].numReads > 0)
        GenotypeLikelihoods(&cell[i].site[j], i, j, probs, ngsEij, ampEijmat,
                            ampEijpat);
      else
        cell[i].site[j].MLmatAllele = cell[i].site[j].MLpatAllele = MISSING;
    }
  }

  /* add doublets and reads from them */
  if (meanDoubletRate > 0)
    MakeDoublets(probs, ngsEij, ampEijmat, ampEijpat, mainThreadRandDev);

  /* count alternate alleles in total reads */
  for (j = 0; j < numSites; j++) {
    countReadsA = countReadsC = countReadsG = countReadsT = 0;
    for (i = 0; i < numCells + 1; i++) {
      if (cell[i].hasDoublet == YES) {
        countReadsA += cell[i].site[j].readCountDoublet[A];
        countReadsC += cell[i].site[j].readCountDoublet[C];
        countReadsG += cell[i].site[j].readCountDoublet[G];
        countReadsT += cell[i].site[j].readCountDoublet[T];
      } else {
        countReadsA += cell[i].site[j].readCount[A];
        countReadsC += cell[i].site[j].readCount[C];
        countReadsG += cell[i].site[j].readCount[G];
        countReadsT += cell[i].site[j].readCount[T];
      }
    }

    numAltAlleles = 0;
    if (countReadsA > 0 && allSites[j].referenceAllele != A)
      allSites[j].alternateAlleles[numAltAlleles++] = A;
    if (countReadsC > 0 && allSites[j].referenceAllele != C)
      allSites[j].alternateAlleles[numAltAlleles++] = C;
    if (countReadsG > 0 && allSites[j].referenceAllele != G)
      allSites[j].alternateAlleles[numAltAlleles++] = G;
    if (countReadsT > 0 && allSites[j].referenceAllele != T)
      allSites[j].alternateAlleles[numAltAlleles++] = T;
    allSites[j].numAltAlleles = numAltAlleles;

    allSites[j].readCountA = countReadsA;
    allSites[j].readCountC = countReadsC;
    allSites[j].readCountG = countReadsG;
    allSites[j].readCountT = countReadsT;
    allSites[j].readCountACGT =
        countReadsA + countReadsC + countReadsG + countReadsT;
  }

  /* count ML genotyping errors */
  countMLgenotypeErrors = 0;
  countCalledGenotypes = 0;
  for (i = 0; i < numCells + 1; i++) {
    for (j = 0; j < numSites; j++) {
      c = &cell[i].site[j];
      trueMaternalAllele = c->trueMaternalAllele;
      truePaternalAllele = c->truePaternalAllele;
      if (cell[i].hasDoublet == NO) {
        MLmatAllele = c->MLmatAllele;
        MLpatAllele = c->MLpatAllele;
      } else {
        MLmatAllele = c->MLmatAlleleDoublet;
        MLpatAllele = c->MLpatAlleleDoublet;
      }

      if (MLmatAllele == MISSING)
        equalGenotypes = YES;
      else {
        countCalledGenotypes++;
        equalGenotypes = CompareGenotypes(
            trueMaternalAllele, truePaternalAllele, MLmatAllele, MLpatAllele);
      }

      if (equalGenotypes == NO)
        countMLgenotypeErrors++;
    }
  }

  cumCountMLgenotypeErrors +=
      countMLgenotypeErrors / (1.0 * countCalledGenotypes);

  free(probs);
  probs = NULL;
  free(ngsEij);
  ngsEij = NULL;
  free(ampEijmat);
  ampEijmat = NULL;
  free(ampEijpat);
  ampEijpat = NULL;
}

/********************* SiteReadCounts  ************************/
/*
 Generate read counts for a given genotype (cell and site)
 according to three different models considering sequencing and amplification
 error
 */
void SiteReadCounts(CellSiteStr *c, int i, int j, double *probs,
                    double **ngsEij, double **ampEijmat, double **ampEijpat,
                    RandomDevice *rd) {
  int k, l, m;
  int read, template;
  int numMaternalReads, numPaternalReads;
  int numReads;
  int *readCount;
  int maternalAllele, paternalAllele;
  int thereIsMaternalAllele, thereIsPaternalAllele;
  double maternalSiteAmplificationError, paternalSiteAmplificationError;
  int ampErrorMaternalAllele, ampErrorPaternalAllele, goodTemplate, badTemplate;
  int readsLeft, badReads;
  double cumProb;
  double uniform, goodTemplateProb, badTemplateProb, noSeqError;

  readCount = c->readCount;

  thereIsMaternalAllele = YES;
  thereIsPaternalAllele = YES;

  maternalAllele = data[MATERNAL][i][j];
  if (maternalAllele == ADO || maternalAllele == DELETION)
    thereIsMaternalAllele = NO;

  paternalAllele = data[PATERNAL][i][j];
  if (paternalAllele == ADO || paternalAllele == DELETION)
    thereIsPaternalAllele = NO;

  /* initialize some variables */
  noSeqError = 0;
  badTemplateProb = 0;
  goodTemplateProb = 0;
  ampErrorMaternalAllele = -9;
  ampErrorPaternalAllele = -9;
  badTemplate = -9;
  goodTemplate = -9;
  maternalSiteAmplificationError = 0;
  paternalSiteAmplificationError = 0;

  /* number of reads will follow a Poisson or a Negative Binomial distribution
   * (overdispersed coverage) around mean coverage, and are randomly assigned to
   * maternal/paternal chromosome */
  numMaternalReads = 0;
  numPaternalReads = 0;

  if (rateVarCoverage == YES) /* Negative Binomial */
  {
    if (thereIsMaternalAllele == YES &&
        thereIsPaternalAllele == YES) /* both alleles are present */
    {
      numMaternalReads =
          RandNegativeBinomial(allelicImbalance * coverage, alphaCoverage, rd);
      numPaternalReads = RandNegativeBinomial(
          (1.0 - allelicImbalance) * coverage, alphaCoverage, rd);
    } else if (thereIsMaternalAllele == YES &&
               thereIsPaternalAllele == NO) /* only maternal allele */
    {
      numMaternalReads = RandNegativeBinomial(
          haploidCoverageReduction * coverage, alphaCoverage, rd);
      numPaternalReads = 0;
    } else if (thereIsMaternalAllele == NO &&
               thereIsPaternalAllele == YES) /* only paternal allele */
    {
      numMaternalReads = 0;
      numPaternalReads = RandNegativeBinomial(
          haploidCoverageReduction * coverage, alphaCoverage, rd);
    } else /* no allele is present (locus dropout/deletion) */
    {
      numMaternalReads = 0;
      numPaternalReads = 0;
    }
  } else /* Poisson distributed */
  {
    if (thereIsMaternalAllele == YES &&
        thereIsPaternalAllele == YES) /* both alleles are present */
    {
      numMaternalReads = RandPoisson(allelicImbalance * coverage, rd);
      numPaternalReads = RandPoisson((1.0 - allelicImbalance) * coverage, rd);
    } else if (thereIsMaternalAllele == YES &&
               thereIsPaternalAllele == NO) /* only maternal allele */
    {
      numMaternalReads = RandPoisson(haploidCoverageReduction * coverage, rd);
      numPaternalReads = 0;
    } else if (thereIsMaternalAllele == NO &&
               thereIsPaternalAllele == YES) /* only paternal allele */
    {
      numMaternalReads = 0;
      numPaternalReads = RandPoisson(haploidCoverageReduction * coverage, rd);
    } else /* no allele is present (locus dropout/deletion) */
    {
      numMaternalReads = 0;
      numPaternalReads = 0;
    }
  }

  numReads = numMaternalReads + numPaternalReads;

  /* read count simulation */
  readCount[A] = readCount[C] = readCount[G] = readCount[T] = 0;

  if (numReads > 0) {
    /* choose amplification error for this site. In this case we will have a
     * proportion of alternative templates depending on whether the error
     * occurred early or late in the amplification (we use a Beta distribution
     * to modelize this proportion) */
    if (meanAmplificationError > 0) {
      maternalSiteAmplificationError =
          RandBeta(meanAmplificationError, varAmplificationError, rd);
      paternalSiteAmplificationError =
          RandBeta(meanAmplificationError, varAmplificationError, rd);
    }

    /* initialize ampEijmat */
    for (k = 0; k < 4; k++)
      for (l = 0; l < 4; l++) {
        if (k == l)
          ampEijmat[k][l] = 1.0 - maternalSiteAmplificationError;
        else
          ampEijmat[k][l] =
              Eij[k][l] / cumEij[k][3] * maternalSiteAmplificationError;
      }

    /* initialize ampEijpat */
    for (k = 0; k < 4; k++)
      for (l = 0; l < 4; l++) {
        if (k == l)
          ampEijpat[k][l] = 1.0 - paternalSiteAmplificationError;
        else
          ampEijpat[k][l] =
              Eij[k][l] / cumEij[k][3] * paternalSiteAmplificationError;
      }

    if (simulateOnlyTwoTemplates == YES) {
      /* choose base for wrong maternal template */
      ampErrorMaternalAllele = maternalAllele;
      if (maternalSiteAmplificationError > 0 && thereIsMaternalAllele == YES) {
        uniform = RandUniform(rd) * cumEij[maternalAllele][3];
        for (m = 0; m < 4; m++)
          if (uniform <= cumEij[maternalAllele][m]) {
            ampErrorMaternalAllele = m;
            break;
          }
      }
      /* choose base for wrong paternal template */
      ampErrorPaternalAllele = paternalAllele;
      if (paternalSiteAmplificationError > 0 && thereIsPaternalAllele == YES) {
        uniform = RandUniform(rd) * cumEij[paternalAllele][3];
        for (m = 0; m < 4; m++)
          if (uniform <= cumEij[paternalAllele][m]) {
            ampErrorPaternalAllele = m;
            break;
          }
      }
    }

    /* probabilities for the different types of reads  */
    if (sequencingError > 0 || maternalSiteAmplificationError > 0 ||
        paternalSiteAmplificationError > 0) {
      if (thereIsMaternalAllele == YES) {
        if (simulateOnlyTwoTemplates == YES) {
          /* probabilities for maternal reads  (assuming a chosen template )*/
          noSeqError = 1.0 - sequencingError;
          badTemplateProb = maternalSiteAmplificationError;
          goodTemplateProb = 1.0 - maternalSiteAmplificationError;
          badTemplate = ampErrorMaternalAllele;
          goodTemplate = maternalAllele;
          for (read = 0; read < 4; read++) {
            if (read == maternalAllele)
              probs[read] = goodTemplateProb * noSeqError +
                            badTemplateProb * ngsEij[badTemplate][read];
            else if (read == ampErrorMaternalAllele)
              probs[read] = goodTemplateProb * ngsEij[goodTemplate][read] +
                            badTemplateProb * noSeqError;
            else
              probs[read] = goodTemplateProb * ngsEij[goodTemplate][read] +
                            badTemplateProb * ngsEij[badTemplate][read];
          }
        } else {
          /* probabilities for maternal reads integrating over all potential
           * templates) */
          for (read = 0; read < 4; read++) {
            probs[read] = 0;
            for (template = 0; template < 4; template++) {
              probs[read] +=
                  ampEijmat[maternalAllele][template] * ngsEij[template][read];
            }
          }
        }

        cumProb = 0;
        readsLeft = numMaternalReads;
        /* multinomial assignment of reads */
        for (read = 0; read < 4; read++) {
          if (read != maternalAllele) {
            badReads = RandBinomial(probs[read] / (1 - cumProb), readsLeft, rd);
            readsLeft -= badReads;
            readCount[read] += badReads;
            cumProb += probs[read];
          }
        }
        readCount[maternalAllele] += readsLeft; /* good reads */
      }

      if (thereIsPaternalAllele == YES) {
        if (simulateOnlyTwoTemplates == YES) {
          /* probabilities for paternal reads  (assuming a chosen template )*/
          badTemplateProb = paternalSiteAmplificationError;
          goodTemplateProb = 1.0 - paternalSiteAmplificationError;
          badTemplate = ampErrorPaternalAllele;
          goodTemplate = paternalAllele;

          for (read = 0; read < 4; read++) {
            if (read == paternalAllele)
              probs[read] = goodTemplateProb * noSeqError +
                            badTemplateProb * ngsEij[badTemplate][read];
            else if (read == ampErrorPaternalAllele)
              probs[read] = goodTemplateProb * ngsEij[goodTemplate][read] +
                            badTemplateProb * noSeqError;
            else
              probs[read] = goodTemplateProb * ngsEij[goodTemplate][read] +
                            badTemplateProb * ngsEij[badTemplate][read];
          }
        } else {
          /* probabilities for paternal reads  (integrating over all potential
           * templates) */
          for (read = 0; read < 4; read++) {
            probs[read] = 0;
            for (template = 0; template <4; template ++)
              probs[read] +=
                  ampEijpat[paternalAllele][template] * ngsEij[template][read];
          }
        }

        cumProb = 0;
        readsLeft = numPaternalReads;
        /* multinomial assignment of reads */
        for (read = 0; read < 4; read++) {
          if (read != paternalAllele) {
            badReads = RandBinomial(probs[read] / (1 - cumProb), readsLeft, rd);
            readsLeft -= badReads;
            readCount[read] += badReads;
            cumProb += probs[read];
          }
        }
        readCount[paternalAllele] += readsLeft; // good reads
      }
    } else /* no sequencing or amplification error */
    {
      readCount[maternalAllele] += numMaternalReads;
      readCount[paternalAllele] += numPaternalReads;
    }
  } // numReads > 0

  c->numReads = numReads;
  c->readCount = readCount;
  c->maternalAllele = maternalAllele;
  c->paternalAllele = paternalAllele;
  c->thereIsMaternalAllele = thereIsMaternalAllele;
  c->thereIsPaternalAllele = thereIsPaternalAllele;
  c->maternalSiteAmplificationError = maternalSiteAmplificationError;
  c->paternalSiteAmplificationError = paternalSiteAmplificationError;
}

/* added code */
/********************* GenerateReadCountsExtd ************************/
/* For each site in each cell the program will generate read counts
   given some sequencing depth or coverage. The number of reads follow a
   Negative Binomial distribution around mean coverage. Reads are
   randomly assigned to each nucleotide following a beta-binomial
   distribution. Reads can contain some errors according to sequencing
   and amplification error parameters.
*/
void GenerateReadCountsExtd(int replicate, RandomDevice *rd) {
  int i, j;
  int countReadsA, countReadsC, countReadsG, countReadsT;
  int numAltAlleles;
  double *probs, **ngsEij, **ampEij;

  countReadsA = countReadsC = countReadsG = countReadsT = 0;

  probs = (double *)calloc(5, sizeof(double));
  if (!probs) {
    fprintf(stderr, "Could not allocate the probs structure\n");
    exit(-1);
  }

  /* allocate error probabilities for NGS */
  ngsEij = (double **)calloc(5, sizeof(double **));
  if (!ngsEij) {
    fprintf(stderr, "Could not allocate the ngsEij structure\n");
    exit(-1);
  }
  for (i = 0; i < 4; i++) {
    ngsEij[i] = (double *)calloc(5, sizeof(double *));
    if (!ngsEij[i]) {
      fprintf(stderr, "Could not allocate the ngsEij[i] structure\n");
      exit(-1);
    }
  }

  /* allocate error probabilities for amplification */
  ampEij = (double **)calloc(5, sizeof(double **));
  if (!ampEij) {
    fprintf(stderr, "Could not allocate the ampEij structure\n");
    exit(-1);
  }
  for (i = 0; i < 4; i++) {
    ampEij[i] = (double *)calloc(5, sizeof(double *));
    if (!ampEij[i]) {
      fprintf(stderr, "Could not allocate the ampEij[i] structure\n");
      exit(-1);
    }
  }

  /* initialize ngsEij now as this will be constant across sites */
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) {
      if (i == j)
        ngsEij[i][j] = 1.0 - sequencingError;
      else
        ngsEij[i][j] = Eij[i][j] / cumEij[i][3] * sequencingError;
    }

  /* if use fixed amplification error, initialize a global ampEij */
  if (fixedAmplificationError) {
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++) {
        if (i == j)
          ampEij[i][j] = 1.0 - meanAmplificationError;
        else
          ampEij[i][j] = Eij[i][j] / cumEij[i][3] * meanAmplificationError;
      }
  }

  /* initialize size factors according to normal distribution */
  RandGaussianArr(numCells + 1, meanSizeFactor, sqrt(varSizeFactor), rd,
                  sizeFactors, YES);

  /* initialize allelicCov and allelicRawVar according to normal distribution */
  RandGaussianArr(numSites, allelicCovMean, sqrt(allelicCovVar), rd, allelicCov,
                  YES);
  RandGaussianArr(numSites, allelicRawVarMean, sqrt(allelicRawVarVar), rd,
                  allelicRawVar, YES);

  /* produce read counts for all sites */
  for (i = 0; i < numCells + 1; i++) {
    cellExtd[i].sizeFactor = sizeFactors[i];
    for (j = 0; j < numSites; j++) {
      if (i == 0) {
        allSites[j].allelicCov = allelicCov[j];
        allSites[j].allelicRawVar = allelicRawVar[j];
      }

      SiteReadCountsExtd(&cellExtd[i].site[j], sizeFactors[i], j, probs,
                         ngsEij, ampEij, rd);
    }
  }

  /* add doublets and reads from them */
  if (meanDoubletRate > 0)
    MakeDoubletsExtd(probs, ngsEij, ampEij, mainThreadRandDev);

  /* count alternate alleles in total reads */
  for (j = 0; j < numSites; j++) {
    countReadsA = countReadsC = countReadsG = countReadsT = 0;
    for (i = 0; i < numCells + 1; i++) {
      countReadsA += cellExtd[i].site[j].readCount[A];
      countReadsC += cellExtd[i].site[j].readCount[C];
      countReadsG += cellExtd[i].site[j].readCount[G];
      countReadsT += cellExtd[i].site[j].readCount[T];

      if (cellExtd[i].isDoublet == YES)
      {
        countReadsA += cellExtd[cellExtd[i].mixedTumorCellNameIdx].site[j].readCount[A];
        countReadsC += cellExtd[cellExtd[i].mixedTumorCellNameIdx].site[j].readCount[C];
        countReadsG += cellExtd[cellExtd[i].mixedTumorCellNameIdx].site[j].readCount[G];
        countReadsT += cellExtd[cellExtd[i].mixedTumorCellNameIdx].site[j].readCount[T];
      }
    }

    numAltAlleles = 0;
    if (countReadsA > 0 && allSites[j].referenceAllele != A)
      allSites[j].alternateAlleles[numAltAlleles++] = A;
    if (countReadsC > 0 && allSites[j].referenceAllele != C)
      allSites[j].alternateAlleles[numAltAlleles++] = C;
    if (countReadsG > 0 && allSites[j].referenceAllele != G)
      allSites[j].alternateAlleles[numAltAlleles++] = G;
    if (countReadsT > 0 && allSites[j].referenceAllele != T)
      allSites[j].alternateAlleles[numAltAlleles++] = T;
    allSites[j].numAltAlleles = numAltAlleles;

    allSites[j].readCountA = countReadsA;
    allSites[j].readCountC = countReadsC;
    allSites[j].readCountG = countReadsG;
    allSites[j].readCountT = countReadsT;
    allSites[j].readCountACGT =
        countReadsA + countReadsC + countReadsG + countReadsT;
  }

  free(probs);
  probs = NULL;
  for (i = 0; i < 4; i++) {
    free(ngsEij[i]);
    ngsEij[i] = NULL;
    free(ampEij[i]);
    ampEij[i] = NULL;
  }
  free(ngsEij);
  ngsEij = NULL;
  free(ampEij);
  ampEij = NULL;
}

/********************* SiteReadCountsExtd ************************/
/*
 Generate read counts for a given genotype (cell and site)
 according to three different models considering sequencing and amplification
 error
 */
void SiteReadCountsExtd(CellSiteStrExtd *c, double sizeFactor, int j,
                        double *probs, double **ngsEij, double **ampEij,
                        RandomDevice *rd) {
  int k, l, read, template;
  int numReads;
  int *readCount;
  int numExistingAlleles, numRemainingAlleles, referenceAllele;
  double siteAmplificationError;
  int genotypeStructure[4];
  double randomNr;

  numReads = 0;
  readCount = c->readCount;
  numExistingAlleles = numRemainingAlleles = c->observedNumExistingAlleles;
  referenceAllele = allSites[j].referenceAllele;
  siteAmplificationError = 0;

  /* number of reads will follow a Negative Binomial distribution
   * (overdispersed coverage) around mean coverage */
  if (seqCovCellCoalStyle == YES) {
    while (numRemainingAlleles > 0)
    {
      numReads += RandNegativeBinomial(allSites[j].allelicCov, allSites[j].allelicRawVar, rd);
      numRemainingAlleles--;
    }
  } else {
    if (numExistingAlleles > 0)
    {
      numReads = GenerateSequencingCoverage(
                     numExistingAlleles, allSites[j].allelicCov, sizeFactor,
                     allSites[j].allelicRawVar, rd);
    }
  }

  /* read count simulation */
  readCount[A] = readCount[C] = readCount[G] = readCount[T] = 0;

  if (numReads > 0) {
    /* choose amplification error for this site. In this case we will have a
     * proportion of alternative templates depending on whether the error
     * occurred early or late in the amplification (we use a Beta distribution
     * to modelize this proportion) */
    if (meanAmplificationError > 0 && !fixedAmplificationError) {
      siteAmplificationError =
          RandBeta(meanAmplificationError, varAmplificationError, rd);
    }

    /* if amplification error is not fixed, initialize a site-wise ampEij */
    if (!fixedAmplificationError) {
      for (k = 0; k < 4; k++)
        for (l = 0; l < 4; l++) {
          if (k == l)
            ampEij[k][l] = 1.0 - siteAmplificationError;
          else
            ampEij[k][l] = Eij[k][l] / cumEij[k][3] * siteAmplificationError;
        }
    }

    /* probabilities for the different types of reads  */
    if (sequencingError > 0 || meanAmplificationError > 0) {
      if (numExistingAlleles > 0) /* beta-binomial assignment */
      {
        GetGenotypeStructure(c->observedAlleles, genotypeStructure);

        /* initialize nucleotideProportion */
        double nucleotideProportion[4];
        for (k = 0; k < 4; ++k) {
          nucleotideProportion[k] =
              genotypeStructure[k] / (double)numExistingAlleles;
        }

        /* initialize nucleotideProportionAfterAmp */
        double nucleotideProportionAfterAmp[4];
        for (read = 0; read < 4; read++) {
          nucleotideProportionAfterAmp[read] = 0;
          for (template = 0; template <4; template ++)
            nucleotideProportionAfterAmp[read] +=
                nucleotideProportion[template] * ampEij[template][read];
        }

        /* initialize the probability of observing each nucleotide after
         * sequencing
         */
        for (read = 0; read < 4; read++) {
          probs[read] = 0;
          for (template = 0; template <4; template ++)
            probs[read] +=
                nucleotideProportionAfterAmp[template] * ngsEij[template][read];
        }

        c->shapeController = shapeController[GetShapeController(
            c->observedAlleles, numExistingAlleles)];

        RandDirichletMultinomial(probs, c->shapeController, numReads,
                                 referenceAllele, 4, readCount, rd);
      } else /* random assignment if the number of existing alleles is 0 */
      {
        for (k = 0; k < numReads; ++k) {
          randomNr = RandUniform(rd) * cumNucleotides[3];
          for (l = 0; l < 4; ++l) {
            if (randomNr <= cumNucleotides[l]) {
              readCount[l]++;
              break;
            }
          }
        }
      }
    } else /* no sequencing and amplification error */
    {
      GetGenotypeStructure(c->observedAlleles, genotypeStructure);

      c->shapeController = shapeController[GetShapeController(
          c->observedAlleles, numExistingAlleles)];

      /* initialize nucleotideProportion */
      double nucleotideProportion[4];
      for (k = 0; k < 4; ++k) {
        nucleotideProportion[k] =
            genotypeStructure[k] / (double)numExistingAlleles;
      }

      RandDirichletMultinomial(nucleotideProportion, c->shapeController,
                               numReads, referenceAllele, 4, readCount, rd);
    }
  } // numReads > 0

  c->numReads = numReads;
  c->fixedAmplificationError = fixedAmplificationError;

  if (fixedAmplificationError) {
    if (meanAmplificationError) {
      c->hasAmplificationError = YES;
    } else {
      c->hasAmplificationError = NO;
    }
    c->amplificationError = meanAmplificationError;
  } else {
    c->hasAmplificationError = YES;
    c->amplificationError = siteAmplificationError;
  }
}

/********************* GenerateSequencingCoverage ********************/
int GenerateSequencingCoverage(int numExistingAlleles, double alleleCoverage,
                               double sizeFactor, double alleleRawVariance,
                               RandomDevice *rd) {
  double mean, variance, dispersion;

  mean =
      (numExistingAlleles + numericalStabilizer) * alleleCoverage * sizeFactor;
  variance = mean + pow(sizeFactor, 2) *
                        pow(numExistingAlleles + numericalStabilizer, 2) *
                        alleleRawVariance;
  dispersion = pow(mean, 2) / (variance - mean);
  return RandNegativeBinomial(mean, dispersion, rd);
}

/********************* RandomDirichletMultinomial ********************/
// void RandomDirichletMultinomial(double *means, double overdispersion,
//                                 int numTrials, int pivot, size_t size,
//                                 long int *seed, int *readCounts) {
//   double gammas[size], probs[size], sum = 0.0;
//   int badReads, readsLeft;
//
//   Normalize(means, size);
//
//   for (size_t i = 0; i < size; i++) {
//     gammas[i] = RandomGamma(means[i] * overdispersion, seed);
//     sum += gammas[i];
//
//     readCounts[i] = 0;
//   }
//
//   for (size_t i = 0; i < size; i++)
//     probs[i] = gammas[i] / sum;
//
//   readsLeft = numTrials;
//   for (size_t i = 0; i < size; i++) {
//     if (i != pivot) {
//       badReads = RandomBinomial(probs[i], numTrials, seed);
//
//       if (readsLeft - badReads >= 0)
//         readCounts[i] = badReads;
//       else {
//         badReads = readsLeft;
//         readCounts[i] = badReads;
//       }
//
//       readsLeft -= badReads;
//
//       if (readsLeft == 0)
//         break;
//     }
//   }
//   readCounts[pivot] = readsLeft;
// }

/********************* RandomBetaBinomial ********************/
// int RandomBetaBinomial(double mean, double overdispersion, int numTrials,
//                        long int *seed) {
//   double alpha, beta, gamma1, gamma2, prob;
//   int sum;
//
//   alpha = mean * overdispersion;
//   beta = overdispersion * (1 - mean);
//
//   gamma1 = RandomGamma(alpha, seed);
//   gamma2 = RandomGamma(beta, seed);
//   prob = gamma1 / (gamma1 + gamma2);
//
//   sum = 0;
//   for (int i = 0; i < numTrials; i++) {
//     if (RandomUniform(seed) < prob)
//       sum++;
//   }
//   return sum;
// }

/********************* GenotypeLikelihoods  ************************/
/*
 Calculate genotype likelihoods for a given site and cell given the simulated
 reads according to three different models

 **	Genotype likelihoods can be calculated given the read counts and the
 sequencing error (see Korneliussen 2013 BMCBioinf)
 */
// void GenotypeLikelihoodsExtd(CellSiteStrExtd *c, int i, int j, double
// **ngsEij,
//                              double **ampEij) {
//   int k, l;
//   int read, numReads;
//   const int *readCount;
//   int numExistingAlleles;
//   double amplificationError;
//   int template;
//   double maxLikelihood;
//   int MLCode;
//   double noSequencingError;
//   double *log10GenotypeLikelihood, *scaledLog10GenotypeLikelihood;
//   int genotypeStructure[4];
//   double nucleotideProportion[4];
//   double log10Likelihood;
//   int uniqueVals[NUM_CAN_ML_GENOTYPES], uniqueValsNum, endPosition, inFlag;
//   const int *candidateMLGenotypes;
//
//   numReads = c->numReads;
//   readCount = c->readCount;
//   log10GenotypeLikelihood = c->log10GenotypeLikelihood;
//   scaledLog10GenotypeLikelihood = c->scaledLog10GenotypeLikelihood;
//   numExistingAlleles = c->observedNumExistingAlleles;
//   amplificationError = c->amplificationError;
//   noSequencingError = 1 - sequencingError;
//   uniqueValsNum = 0;
//   candidateMLGenotypes = candidateMLGenotypesMD;
//
//   for (k = 0; k < NUM_CAN_ML_GENOTYPES; k++) {
//     uniqueVals[k] = 0;
//   }
//
//   /************************** MODEL 0 (GATK)
//    * *************************************/
//   /* MODEL 0: genotype log10 likelihoods according to the observed read
//    * assuming a single sequencing error independent of the nucleotides
//    * involved, and no amplification error */
//   /* extended from Korneliussen 2013 BMCBioinf
//    p(D|G=A1,A2) = multiply over reads [ 1/2 p(read i|A1) + 1/2 p(read i|A2)]
//    where p(read|A) e/3 if b!=A  or 1-e if b=A
//    N.B: ANGSD seems to use natural logs!
//    */
//   if (amplificationError == 0) {
//     for (k = 0; k < NUM_CAN_ML_GENOTYPES; ++k) {
//       GetGenotypeStructure(GetValIntToIntArr(genotypeHashMapIntToIntArr,
//                                              candidateMLGenotypes[k]),
//                            genotypeStructure);
//
//       /* initialize nucleotideProportion */
//       for (l = 0; l < 4; ++l) {
//         nucleotideProportion[l] =
//             genotypeStructure[l] / (double)numExistingAlleles;
//       }
//
//       log10GenotypeLikelihood[k] = -0.0;
//
//       for (read = 0; read < 4; ++read) // each unique read in readCount
//       {
//         if (readCount[read] > 0) {
//           log10Likelihood = -0.0;
//
//           for (template = 0; template <4; ++template) // each given
//           nucleotide
//           {
//             if (nucleotideProportion[template] > 0) {
//               if (read == template) {
//                 log10Likelihood +=
//                     nucleotideProportion[template] * noSequencingError;
//               } else {
//                 log10Likelihood +=
//                     nucleotideProportion[template] * sequencingError / 3.0;
//               }
//             }
//           }
//
//           log10GenotypeLikelihood[k] +=
//               readCount[read] * log10(log10Likelihood);
//         }
//       }
//     }
//
//     /* print information */
//     if (noisy > 3) {
//       fprintf(stderr, "\ncell %d site %d: (code) %d, (genotype)", i + 1, j +
//       1,
//               c->observedCode);
//       PrintGenotypeInLettersToFile(NULL, numExistingAlleles,
//       c->observedAlleles,
//                                    1);
//       fprintf(stderr, "  |  reads: A:%d C:%d G:%d T:%d  |  total:%d ",
//               readCount[A], readCount[C], readCount[G], readCount[T],
//               numReads);
//       fprintf(stderr, "\namplificationError = %.9f", amplificationError);
//       for (k = 0; k < NUM_CAN_ML_GENOTYPES; ++k) {
//         fprintf(stderr, "\n(GATK)log10 GenotypeLikelihood");
//         PrintGenotypeInLettersToFile(
//             NULL,
//             altModel == FiniteMuInDel
//                 ? GetExistingAllelesNrMDI(candidateMLGenotypes[k])
//                 : GetExistingAllelesNrMD(candidateMLGenotypes[k]),
//             GetValIntToIntArr(genotypeHashMapIntToIntArr,
//                               candidateMLGenotypes[k]),
//             1);
//         fprintf(stderr, " = %.9f", log10GenotypeLikelihood[k]);
//       }
//       fprintf(stderr, "\n");
//     }
//   } // model 0
//
//   /************************** MODEL 1 (4-templates)
//    * *************************************/
//   /* MODEL 1 (4-templates): genotype log10 likelihoods according to the
//    * observed read, assuming a single sequencing error independent of the
//    * nucleotides involved, and amplification error with all 4 templates */
//   if (amplificationError > 0) {
//     for (k = 0; k < NUM_CAN_ML_GENOTYPES; ++k) {
//       GetGenotypeStructure(GetValIntToIntArr(genotypeHashMapIntToIntArr,
//                                              candidateMLGenotypes[k]),
//                            genotypeStructure);
//
//       /* allocate and initialize nucleotideProportion */
//       for (l = 0; l < 4; ++l) {
//         nucleotideProportion[l] =
//             genotypeStructure[l] / (double)numExistingAlleles;
//       }
//
//       log10GenotypeLikelihood[k] = -0.0;
//
//       for (read = 0; read < 4; ++read) // each unique read in readCount
//       {
//         if (readCount[read] > 0) {
//           log10Likelihood = -0.0;
//
//           for (template = 0; template <4; ++template) // each given
//           nucleotide
//           {
//             if (nucleotideProportion[template] > 0) {
//               for (l = 0; l < 4; ++l) // each intermediate nucleotide
//               {
//                 log10Likelihood += ampEij[template][l] * ngsEij[l][read];
//               }
//
//               log10Likelihood *= nucleotideProportion[template];
//             }
//           }
//
//           log10GenotypeLikelihood[k] +=
//               readCount[read] * log10(log10Likelihood);
//         }
//       }
//     }
//
//     /* print information */
//     if (noisy > 3) {
//       fprintf(stderr, "\ncell %d site %d: (code) %d, (genotype)", i + 1, j +
//       1,
//               c->observedCode);
//       PrintGenotypeInLettersToFile(NULL, numExistingAlleles,
//       c->observedAlleles,
//                                    1);
//       fprintf(stderr, "  |  reads: A:%d C:%d G:%d T:%d  |  total:%d ",
//               readCount[A], readCount[C], readCount[G], readCount[T],
//               numReads);
//       fprintf(stderr, "\namplificationError = %.9f", amplificationError);
//       for (k = 0; k < NUM_CAN_ML_GENOTYPES; ++k) {
//         fprintf(stderr, "\n(4T)log10 GenotypeLikelihood");
//         PrintGenotypeInLettersToFile(
//             NULL,
//             altModel == FiniteMuInDel
//                 ? GetExistingAllelesNrMDI(candidateMLGenotypes[k])
//                 : GetExistingAllelesNrMD(candidateMLGenotypes[k]),
//             GetValIntToIntArr(genotypeHashMapIntToIntArr,
//                               candidateMLGenotypes[k]),
//             1);
//         fprintf(stderr, " = %.9f", log10GenotypeLikelihood[k]);
//       }
//       fprintf(stderr, "\n");
//     }
//   } // model 1
//
//   /* find max log10 likelihood  */
//   maxLikelihood = log10GenotypeLikelihood[0];
//   MLCode = candidateMLGenotypes[0];
//   for (k = 0; k < NUM_CAN_ML_GENOTYPES; k++) {
//     if (log10GenotypeLikelihood[k] > maxLikelihood) {
//       maxLikelihood = log10GenotypeLikelihood[k];
//       MLCode = candidateMLGenotypes[k];
//     }
//   }
//
//   /* rescale to likelihood ratios */
//   for (k = 0; k < NUM_CAN_ML_GENOTYPES; k++) {
//     scaledLog10GenotypeLikelihood[k] =
//         log10GenotypeLikelihood[k] - maxLikelihood;
//   }
//
//   if (noisy > 3) {
//     for (k = 0; k < NUM_CAN_ML_GENOTYPES; ++k) {
//       fprintf(stderr, "\nscalog10 GenotypeLikelihood");
//       PrintGenotypeInLettersToFile(
//           NULL,
//           altModel == FiniteMuInDel
//               ? GetExistingAllelesNrMDI(candidateMLGenotypes[k])
//               : GetExistingAllelesNrMD(candidateMLGenotypes[k]),
//           GetValIntToIntArr(genotypeHashMapIntToIntArr,
//                             candidateMLGenotypes[k]),
//           1);
//       fprintf(stderr, " = %.9f", scaledLog10GenotypeLikelihood[k]);
//     }
//     fprintf(stderr, "\n\nMaximum likelihood genotype:");
//     PrintGenotypeInLettersToFile(
//         NULL,
//         altModel == FiniteMuInDel
//             ? GetExistingAllelesNrMDI(candidateMLGenotypes[k])
//             : GetExistingAllelesNrMD(candidateMLGenotypes[k]),
//         GetValIntToIntArr(genotypeHashMapIntToIntArr, MLCode), 1);
//     fprintf(stderr, "\nMaximum likelihood genotype code: %d", MLCode);
//     fprintf(stderr, "\n");
//   }
//
//   /* get unique log-likelihoods in log10GenotypeLikelihood */
//   for (k = 0; k < NUM_CAN_ML_GENOTYPES; ++k) {
//     endPosition = uniqueValsNum;
//     inFlag = 0;
//     for (l = 0; l < endPosition; ++l) {
//       if (fabs(log10GenotypeLikelihood[k] -
//                log10GenotypeLikelihood[uniqueVals[l]]) < 1e-8) {
//         inFlag = 1;
//         break;
//       }
//     }
//
//     if (!inFlag) {
//       uniqueVals[uniqueValsNum++] = k;
//     }
//   }
//
//   /* allocate uniqueLog10GenotypeLikelihood and
//    * uniqueScaledLog10GenotypeLikelihood */
//   c->uniqueLog10GenotypeLikelihood =
//       (double *)calloc(uniqueValsNum + 1, sizeof(double));
//   if (!c->uniqueLog10GenotypeLikelihood) {
//     fprintf(stderr, "Could not allocate uniqueLog10GenotypeLikelihood\n");
//     exit(-1);
//   }
//   c->uniqueScaledLog10GenotypeLikelihood =
//       (double *)calloc(uniqueValsNum + 1, sizeof(double));
//   if (!c->uniqueScaledLog10GenotypeLikelihood) {
//     fprintf(stderr, "Could not allocate
//     uniqueScaledLog10GenotypeLikelihood\n"); exit(-1);
//   }
//
//   /* store unique likelihoods */
//   for (l = 0; l < uniqueValsNum; ++l) {
//     c->uniqueLog10GenotypeLikelihood[l] =
//         log10GenotypeLikelihood[uniqueVals[l]];
//     c->uniqueScaledLog10GenotypeLikelihood[l] =
//         scaledLog10GenotypeLikelihood[uniqueVals[l]];
//   }
//
//   c->uniqueGenotypeLikelihoodNr = uniqueValsNum;
//   c->MLCode = MLCode;
//   c->MLNumExistingAlleles = altModel == FiniteMuInDel
//                                 ? GetExistingAllelesNrMDI(MLCode)
//                                 : GetExistingAllelesNrMD(MLCode);
//   c->MLAlleles = GetValIntToIntArr(genotypeHashMapIntToIntArr, MLCode);
// }
/* end added code */

/********************* GenotypeLikelihoods  ************************/
/*
 Calculate genotype likelihoods for a given site and cell given the simulated
 reads according to three different models

 **	Genotype likelihoods can be calculated given the read counts and the
 sequencing error (see Korneliussen 2013 BMCBioinf)

 00 01 02 03 11 12 13 22 23 33
 AA AC AG AT CC CG CT GG GT TT
 */

void GenotypeLikelihoods(CellSiteStr *c, int i, int j, double *probs,
                         double **ngsEij, double **ampEijmat,
                         double **ampEijpat) {
  int k, l, a1, a2;
  int read, numReads;
  int *readCount;
  int maternalAllele, paternalAllele;
  int thereIsMaternalAllele, thereIsPaternalAllele;
  double maternalSiteAmplificationError, paternalSiteAmplificationError;
  int template, template1, template2;
  int debug_GL, fix_parameters;
  double maxLike, pReadGivenA1, pReadGivenA2;
  double homozygote_for_read, heterozygote_for_read, no_allele_for_read;
  double **genLike, **scaledGenLike;
  double MLmatAllele, MLpatAllele;

  numReads = c->numReads;
  readCount = c->readCount;
  genLike = c->genLike;
  scaledGenLike = c->scaledGenLike;
  maternalAllele = c->maternalAllele;
  paternalAllele = c->paternalAllele;
  thereIsMaternalAllele = c->thereIsMaternalAllele;
  thereIsPaternalAllele = c->thereIsPaternalAllele;
  maternalSiteAmplificationError = c->maternalSiteAmplificationError;
  paternalSiteAmplificationError = c->paternalSiteAmplificationError;

  debug_GL = NO;
  fix_parameters = NO;
  if (fix_parameters == YES) {
    maternalAllele = c->maternalAllele = A;
    paternalAllele = c->paternalAllele = C;
    thereIsMaternalAllele = c->thereIsMaternalAllele = YES;
    thereIsPaternalAllele = c->thereIsPaternalAllele = YES;
    readCount[A] = 76;
    readCount[C] = 71;
    readCount[G] = 6;
    readCount[T] = 6;

    readCount[A] = 3;
    readCount[C] = 2;
    readCount[G] = 0;
    readCount[T] = 0;

    numReads = c->numReads =
        readCount[A] + readCount[C] + readCount[G] + readCount[T];
    sequencingError = 0.0;
    maternalSiteAmplificationError = c->maternalSiteAmplificationError = 0.125;
    paternalSiteAmplificationError = c->paternalSiteAmplificationError = 0.125;

    /* reinitialize ngsEij  */
    for (k = 0; k < 4; k++)
      for (l = 0; l < 4; l++) {
        if (k == l)
          ngsEij[k][l] = 1.0 - sequencingError;
        else
          ngsEij[k][l] = Eij[k][l] / cumEij[k][3] * sequencingError;
      }

    /* reinitialize ampEijmat */
    for (k = 0; k < 4; k++)
      for (l = 0; l < 4; l++) {
        if (k == l)
          ampEijmat[k][l] = 1.0 - maternalSiteAmplificationError;
        else
          ampEijmat[k][l] =
              Eij[k][l] / cumEij[k][3] * maternalSiteAmplificationError;
      }

    /* reinitialize ampEijpat */
    for (k = 0; k < 4; k++)
      for (l = 0; l < 4; l++) {
        if (k == l)
          ampEijpat[k][l] = 1.0 - paternalSiteAmplificationError;
        else
          ampEijpat[k][l] =
              Eij[k][l] / cumEij[k][3] * paternalSiteAmplificationError;
      }
  }

  /************************** MODEL 0 (GATK)
   * *************************************/
  /* MODEL 0: genotype log10 likelihoods according to the observed read
   * assuming a single sequencing error independent of the nucleotides
   * involved, and no amplification error */
  /* extended from Korneliussen 2013 BMCBioinf
   p(D|G=A1,A2) = multiply over reads [ 1/2 p(read i|A1) + 1/2 p(read i|A2)]
   where p(read|A) e/3 if b!=A  or 1-e if b=A
   N.B: ANGSD seems to use natural logs!
   */
  if (maternalSiteAmplificationError == 0 &&
      paternalSiteAmplificationError == 0) {
    homozygote_for_read = log10(1.0 - sequencingError);
    no_allele_for_read = log10(sequencingError / 3.0);
    heterozygote_for_read =
        log10((1.0 - sequencingError) / 2.0 + sequencingError / 6.0);

    for (a1 = 0; a1 < 4; a1++)    // a1 is maternal
      for (a2 = a1; a2 < 4; a2++) // a2 is paternal
      {
        genLike[a1][a2] = -0.0;

        for (k = 0; k < 4; k++) {
          if (readCount[k] > 0) {
            if (k == a1 && k == a2)
              genLike[a1][a2] += readCount[k] * homozygote_for_read;
            else if (k != a1 && k != a2)
              genLike[a1][a2] += readCount[k] * no_allele_for_read;
            else
              genLike[a1][a2] += readCount[k] * heterozygote_for_read;
            // fprintf (stderr, "\nGATK read=%c gl[%c][%c] = %lf",
            // WhichNuc(k), WhichNuc(a1), WhichNuc(a2), genLike[a1][a2]);
          }
        }
      }

    if (debug_GL == YES) {
      fprintf(stderr, "\ncell %d site %d: %c%c", i + 1, j + 1,
              WhichNuc(maternalAllele), WhichNuc(paternalAllele));
      fprintf(stderr, "|\treads: A:%d C:%d G:%d T:%d  | total:%d ",
              readCount[A], readCount[C], readCount[G], readCount[T], numReads);
      fprintf(stderr, "\nmatAmpError = %f", maternalSiteAmplificationError);
      fprintf(stderr, "\npatAmpError = %f", paternalSiteAmplificationError);
      for (a1 = 0; a1 < 4; a1++)
        for (a2 = a1; a2 < 4; a2++)
          fprintf(stderr, "\n(GATK)log10 GL[%c][%c] = %lf", WhichNuc(a1),
                  WhichNuc(a2), genLike[a1][a2]);
      fprintf(stderr, "\n");
    }
  } // model 0

  /************************** MODEL 1 (4-templates)
   * *************************************/
  /* MODEL 1 (4-templates): genotype log10 likelihoods according to the
   * observed read, assuming a single sequencing error independent of the
   * nucleotides involved, and amplification error with all 4 templates */
  if (simulateOnlyTwoTemplates == NO && (maternalSiteAmplificationError > 0 ||
                                         paternalSiteAmplificationError > 0)) {
    for (a1 = 0; a1 < 4; a1++)    // a1 is maternal
      for (a2 = a1; a2 < 4; a2++) // a2 is paternal
      {
        genLike[a1][a2] = -0.0;
        for (read = 0; read < 4; read++) {
          pReadGivenA1 = pReadGivenA2 = 0.0;
          if (readCount[read] > 0) {
            for (template = 0; template <4; template ++) {
              pReadGivenA1 += ampEijmat[a1][template] * ngsEij[template][read];
              pReadGivenA2 += ampEijpat[a2][template] * ngsEij[template][read];
              // fprintf (stderr, "\n ampEijmat[%c]][%c]=%lf
              // ngsEij[%c]][%c]=%lf
              // ", WhichNuc(a1), WhichNuc(template), ampEijmat[a1][template],
              // WhichNuc(a2), WhichNuc(template), ngsEij[template][read]);
            }
            genLike[a1][a2] += readCount[read] *
                               log10(allelicImbalance * pReadGivenA1 +
                                     (1.0 - allelicImbalance) * pReadGivenA2);
          }
          // fprintf (stderr, "\n***read=%c gl[%c][%c] = %lf  (p1=%lf
          // p2=%lf)", WhichNuc(read), WhichNuc(a1), WhichNuc(a2),
          // genLike[a1][a2],pReadGivenA1,pReadGivenA2);
        }
      }

    if (debug_GL == YES) {
      fprintf(stderr, "\ncell %d site %d: %c%c", i + 1, j + 1,
              WhichNuc(maternalAllele), WhichNuc(paternalAllele));
      fprintf(stderr, "|\treads: A:%d C:%d G:%d T:%d  | total:%d", readCount[A],
              readCount[C], readCount[G], readCount[T], numReads);
      fprintf(stderr, "\nmatAmpError = %f", maternalSiteAmplificationError);
      fprintf(stderr, "\npatAmpError = %f", paternalSiteAmplificationError);
      for (a1 = 0; a1 < 4; a1++)
        for (a2 = a1; a2 < 4; a2++)
          fprintf(stderr, "\n(4T)log10 GL[%c][%c] = %lf", WhichNuc(a1),
                  WhichNuc(a2), genLike[a1][a2]);
      fprintf(stderr, "\n");
    }
  } // model 1

  /************************** MODEL 2 (2-templates)
   * *************************************/
  /* MODEL 2 (2-templates): genotype log10 likelihoods according to the
   * observed read, assuming a single sequencing error independent of the
   * nucleotides involved, and amplification error with up to 2 templates */
  if (simulateOnlyTwoTemplates == YES && (maternalSiteAmplificationError > 0 ||
                                          paternalSiteAmplificationError > 0)) {
    for (a1 = 0; a1 < 4; a1++)    // a1 is maternal
      for (a2 = a1; a2 < 4; a2++) // a2 is paternal
      {
        genLike[a1][a2] = -0.0;
        // fprintf (stderr, "\n");
        for (read = 0; read < 4; read++)
          if (readCount[read] > 0) {
            pReadGivenA1 = 0;
            pReadGivenA2 = 0;
            for (template1 = 0; template1 < 4; template1++)
              if (template1 != a1) {
                pReadGivenA1 +=
                    Eij[a1][template1] / cumEij[a1][3] *
                    ((1.0 - maternalSiteAmplificationError) * ngsEij[a1][read] +
                     maternalSiteAmplificationError * ngsEij[template1][read]);
                // fprintf (stderr, "\nread=%c t=[%c][%c] p1=%lf
                // Eijmat[%c][%c]=%lf", WhichNuc(read), WhichNuc(a1),
                // WhichNuc(template1), pReadGivenA1, WhichNuc(a1),
                // WhichNuc(template1), Eij[a1][template1]/cumEij[a1][3]);
              }
            for (template2 = 0; template2 < 4; template2++)
              if (template2 != a2) {
                pReadGivenA2 +=
                    Eij[a2][template2] / cumEij[a2][3] *
                    ((1.0 - paternalSiteAmplificationError) * ngsEij[a2][read] +
                     paternalSiteAmplificationError * ngsEij[template2][read]);
                // fprintf (stderr, "\nread=%c t=[%c][%c] p2=%lf
                // Eijpat[%c][%c]=%lf", WhichNuc(read), WhichNuc(a2),
                // WhichNuc(template2), pReadGivenA2, WhichNuc(a2),
                // WhichNuc(template2), Eij[a2][template2]/cumEij[a2][3]);
              }

            if (pReadGivenA1 > 0 || pReadGivenA2 > 0)
              genLike[a1][a2] += readCount[read] *
                                 log10(allelicImbalance * pReadGivenA1 +
                                       (1.0 - allelicImbalance) * pReadGivenA2);
            // fprintf (stderr, "\ngl[%c][%c] = %lf (p1=%lf p2=%lf)",
            // WhichNuc(a1), WhichNuc(a2), genLike[a1][a2], pReadGivenA1,
            // pReadGivenA2);
          } // read
      }     // a1,a2

    if (debug_GL == YES) {
      fprintf(stderr, "\ncell %d site %d: %c%c", i + 1, j + 1,
              WhichNuc(maternalAllele), WhichNuc(paternalAllele));
      fprintf(stderr, "|\treads: A:%d C:%d G:%d T:%d  | total:%d", readCount[A],
              readCount[C], readCount[G], readCount[T], numReads);
      fprintf(stderr, "\nmatAmpError = %f", maternalSiteAmplificationError);
      fprintf(stderr, "\npatAmpError = %f", paternalSiteAmplificationError);
      for (a1 = 0; a1 < 4; a1++)
        for (a2 = a1; a2 < 4; a2++)
          fprintf(stderr, "\n(2Tb)log10 GL[%c][%c] = %lf", WhichNuc(a1),
                  WhichNuc(a2), genLike[a1][a2]);
      fprintf(stderr, "\n");
    }
  } // model 2b

  /* find max log10 likelihood  */
  maxLike = genLike[A][A];
  MLmatAllele = MLpatAllele = A;
  for (a1 = 0; a1 < 4; a1++)
    for (a2 = a1; a2 < 4; a2++)
      if (genLike[a1][a2] > maxLike) {
        maxLike = genLike[a1][a2];
        MLmatAllele = a1;
        MLpatAllele = a2;
      }

  /* rescale to likelihood ratios */
  for (a1 = 0; a1 < 4; a1++)
    for (a2 = a1; a2 < 4; a2++) {
      scaledGenLike[a1][a2] = genLike[a1][a2] - maxLike;
      genLike[a2][a1] = genLike[a1][a2];
      scaledGenLike[a2][a1] = scaledGenLike[a1][a2];
    }
  // Qphred = -10 log(10) Perror)
  // Perror = 10 ^(-Qphred/10)

  if (debug_GL == YES) {
    // PrintSiteInfo (stderr, SNVsites[snv]);
    for (a1 = 0; a1 < 4; a1++)
      for (a2 = a1; a2 < 4; a2++)
        fprintf(stderr, "\nscalog10 GL[%c][%c] = %lf", WhichNuc(a1),
                WhichNuc(a2), scaledGenLike[a1][a2]);
    fprintf(stderr, "\n");
  }

  c->MLmatAllele = MLmatAllele;
  c->MLpatAllele = MLpatAllele;
}

/********************* MakeDoublets  ************************/
/*
 We simulate the possibility of amplifying/sequencing two cells at once
 by adding read counts from another sampled cell.
 Only tumor cells are considered.
*/

void MakeDoublets(double *probs, double **ngsEij, double **ampEijmat,
                  double **ampEijpat, RandomDevice *rd) {
  int c1, c2;
  int j, k, l;
  int debug_DL;
  int MLmatAllele, MLpatAllele;
  double maxLike;
  double doubletRate;
  CellStr celldoublet;

  debug_DL = NO;
  if (debug_DL == YES)
    fprintf(stderr, "\n\n\n********Doublet read counts and likelihoods");

  /* choose whether a tumor cell will become a doublet */
  for (c1 = 0; c1 < numCells; c1++) {
    doubletRate = RandBeta(meanDoubletRate, varDoubletRate, rd);
    if (RandUniform(rd) < doubletRate) {
      cell[c1].hasDoublet = YES;

      /* select a cell genotype from the tumor sample */
      c2 = RandUniformInt(0, numCells, rd);

      /* allocate temporal memory for this cell genotype */
      celldoublet.site =
          (CellSiteStr *)calloc(numSites + 1, sizeof(CellSiteStr));
      if (!celldoublet.site) {
        fprintf(stderr, "Could not allocate the celldoublet.site structure\n");
        exit(-1);
      }

      for (j = 0; j < numSites; j++) {
        celldoublet.site[j].readCount = (int *)calloc(5, sizeof(int));
        if (!celldoublet.site[j].readCount) {
          fprintf(stderr,
                  "Could not allocate the celldoublet.site[%d].readCount "
                  "structure\n",
                  j);
          exit(-1);
        }

        celldoublet.site[j].genLike = (double **)calloc(5, sizeof(double *));
        if (!celldoublet.site[j].genLike) {
          fprintf(stderr,
                  "Could not allocate the celldoublet.site[%d].genLike "
                  "structure\n",
                  j);
          exit(-1);
        }
        for (k = 0; k < 4; k++) {
          celldoublet.site[j].genLike[k] = (double *)calloc(5, sizeof(double));
          if (!celldoublet.site[j].genLike[k]) {
            fprintf(stderr,
                    "Could not allocate the celldoublet.site[%d].genLike[%d] "
                    "structure\n",
                    j, k);
            exit(-1);
          }
        }

        celldoublet.site[j].scaledGenLike =
            (double **)calloc(5, sizeof(double *));
        if (!celldoublet.site[j].scaledGenLike) {
          fprintf(stderr,
                  "Could not allocate the celldoublet.site[%d].scaledGenLike "
                  "structure\n",
                  j);
          exit(-1);
        }
        for (k = 0; k < 4; k++) {
          celldoublet.site[j].scaledGenLike[k] =
              (double *)calloc(5, sizeof(double));
          if (!celldoublet.site[j].scaledGenLike[k]) {
            fprintf(stderr,
                    "Could not allocate the "
                    "celldoublet.site[%d].scaledGenLike[%d] structure\n",
                    j, k);
            exit(-1);
          }
        }
      }

      /* get read counts and likelihoods for this cell genotype for all sites
       */
      for (j = 0; j < numSites;
           j++) // produce doublet read counts for all sites
      {
        if (debug_DL == YES)
          fprintf(stderr, "\n\n-- New celldoublet reads and likelihoods");
        SiteReadCounts(&celldoublet.site[j], c2, j, probs, ngsEij, ampEijmat,
                       ampEijpat, rd);
        GenotypeLikelihoods(&celldoublet.site[j], c2, j, probs, ngsEij,
                            ampEijmat, ampEijpat);

        if (debug_DL == YES) {
          fprintf(stderr, "\ncelltype %d site %d: %c%c", c2 + 1, j + 1,
                  WhichNuc(celldoublet.site[j].maternalAllele),
                  WhichNuc(celldoublet.site[j].paternalAllele));
          fprintf(stderr, "|\treads: A:%d C:%d G:%d T:%d  | total:%d ",
                  celldoublet.site[j].readCount[A],
                  celldoublet.site[j].readCount[C],
                  celldoublet.site[j].readCount[G],
                  celldoublet.site[j].readCount[T],
                  celldoublet.site[j].numReads);
          fprintf(stderr, "\nmatAmpError = %f",
                  celldoublet.site[j].maternalSiteAmplificationError);
          fprintf(stderr, "\npatAmpError = %f",
                  celldoublet.site[j].paternalSiteAmplificationError);
          for (k = 0; k < 4; k++)
            for (l = k; l < 4; l++)
              fprintf(stderr, "\nlog10 GL[%c][%c] = %lf", WhichNuc(k),
                      WhichNuc(l), celldoublet.site[j].genLike[k][l]);
          fprintf(stderr, "\n");
        }

        /* add doublet information */
        cell[c1].site[j].numReadsDoublet =
            cell[c1].site[j].numReads + celldoublet.site[j].numReads;
        for (k = 0; k < 4; k++)
          cell[c1].site[j].readCountDoublet[k] =
              cell[c1].site[j].readCount[k] + celldoublet.site[j].readCount[k];

        cell[c1].site[j].maternalAlleleDoublet =
            celldoublet.site[j].maternalAllele;
        cell[c1].site[j].paternalAlleleDoublet =
            celldoublet.site[j].paternalAllele;

        // fprintf(stderr,"\nln cell\t\tcelltype\tdoublet");
        /* adding the genotype log likelihoods of two sets of reads (c1 +
         * celldoublet) to get the total doublet likelihoods */
        for (k = 0; k < 4; k++)
          for (l = k; l < 4; l++) {
            if (isinf(cell[c1].site[j].genLike[k][l]) ||
                isinf(celldoublet.site[j].genLike[k][l]))
              cell[c1].site[j].genLikeDoublet[k][l] = -1.0 / 0.0;
            else
              cell[c1].site[j].genLikeDoublet[k][l] =
                  cell[c1].site[j].genLike[k][l] +
                  celldoublet.site[j].genLike[k][l];
          }

        if (debug_DL == YES) {
          fprintf(stderr, "\nCell %d site %d: %c%c", c1 + 1, j + 1,
                  WhichNuc(cell[c1].site[j].maternalAllele),
                  WhichNuc(cell[c1].site[j].paternalAllele));
          fprintf(stderr, "|\treads: A:%d C:%d G:%d T:%d  | total:%d ",
                  cell[c1].site[j].readCount[A], cell[c1].site[j].readCount[C],
                  cell[c1].site[j].readCount[G], cell[c1].site[j].readCount[T],
                  cell[c1].site[j].numReads);
          fprintf(stderr, "\n|reads doublet: A:%d C:%d G:%d T:%d  | total:%d ",
                  cell[c1].site[j].readCountDoublet[A],
                  cell[c1].site[j].readCountDoublet[C],
                  cell[c1].site[j].readCountDoublet[G],
                  cell[c1].site[j].readCountDoublet[T],
                  cell[c1].site[j].numReadsDoublet);
          fprintf(stderr, "\nmatAmpError = %f",
                  cell[c1].site[j].maternalSiteAmplificationError);
          fprintf(stderr, "\npatAmpError = %f",
                  cell[c1].site[j].paternalSiteAmplificationError);
          for (k = 0; k < 4; k++)
            for (l = k; l < 4; l++)
              fprintf(stderr, "\ndoublet log10 GL[%c][%c] = %lf", WhichNuc(k),
                      WhichNuc(l), cell[c1].site[j].genLikeDoublet[k][l]);
          fprintf(stderr, "\n");
        }

        /* reescale doublet genotype likelihood to the maximum likelihood and
         * assign ML alleles */
        /* find max log10 doublet likelihood  */
        maxLike = cell[c1].site[j].genLikeDoublet[A][A];
        MLmatAllele = MLpatAllele = A;
        for (k = 0; k < 4; k++)
          for (l = k; l < 4; l++)
            if (cell[c1].site[j].genLikeDoublet[k][l] > maxLike) {
              maxLike = cell[c1].site[j].genLikeDoublet[k][l];
              MLmatAllele = k;
              MLpatAllele = l;
            }

        /* rescale to likelihood ratios */
        for (k = 0; k < 4; k++)
          for (l = k; l < 4; l++)
            cell[c1].site[j].scaledGenLikeDoublet[k][l] =
                cell[c1].site[j].genLikeDoublet[k][l] - maxLike;

        cell[c1].site[j].MLmatAlleleDoublet = MLmatAllele;
        cell[c1].site[j].MLpatAlleleDoublet = MLpatAllele;

        if (debug_DL == YES) {
          for (k = 0; k < 4; k++)
            for (l = k; l < 4; l++)
              fprintf(stderr, "\nscaled doublet log10 GL[%c][%c] = %lf",
                      WhichNuc(k), WhichNuc(l),
                      cell[c1].site[j].scaledGenLikeDoublet[k][l]);
          fprintf(stderr, "\n");
        }
      }

      // free celldoublet
      for (j = 0; j < numSites; j++) {
        free(celldoublet.site[j].readCount);
        celldoublet.site[j].readCount = NULL;
        for (k = 0; k < 4; k++) {
          free(celldoublet.site[j].genLike[k]);
          celldoublet.site[j].genLike[k] = NULL;
          free(celldoublet.site[j].scaledGenLike[k]);
          celldoublet.site[j].scaledGenLike[k] = NULL;
        }
        free(celldoublet.site[j].genLike);
        celldoublet.site[j].genLike = NULL;
        free(celldoublet.site[j].scaledGenLike);
        celldoublet.site[j].scaledGenLike = NULL;
      }
    } // if we do doublet for c1
  }   // c1
}

/********************* MakeDoubletsExtd  ************************/
/*
 We simulate the possibility of amplifying/sequencing two cells at once
 by adding read counts from another sampled cell.
 Only tumor cells are considered.
*/

void MakeDoubletsExtd (double *probs, double **ngsEij, double **ampEij, RandomDevice *rd)
	{
	int		c1, c2;
	int 	j, doubletIdx;
	double	doubletRate;

  doubletIdx = numCells + 1;

	/* choose whether a tumor cell will become a doublet */
	for (c1 = 0; c1 < numCells; c1++)
		{
		doubletRate = RandBeta(meanDoubletRate, varDoubletRate, rd);
		if (RandUniform(rd) < doubletRate)
			{
			cellExtd[c1].isDoublet = cellExtd[doubletIdx].isDoublet = YES;

			/* select a cell genotype from the tumor sample */
		    c2 = RandUniformInt(0, numCells, rd);
		    cellExtd[c1].mixeTumorCellIdx = c2;
		    cellExtd[c1].mixedTumorCellNameIdx = doubletIdx;

		    cellExtd[doubletIdx].sizeFactor = sizeFactors[c2];
		    cellExtd[doubletIdx].mixeTumorCellIdx = c1;
		    cellExtd[doubletIdx].mixedTumorCellNameIdx = c1;

			/* get read counts and likelihoods for this cell genotype for all sites */
			for (j = 0; j < numSites; j++) // produce doublet read counts for all sites
			{
			    cellExtd[doubletIdx].site[j].shapeController = cellExtd[c2].site[j].shapeController;
			    cellExtd[doubletIdx].site[j].trueCode = cellExtd[c2].site[j].trueCode;
			    cellExtd[doubletIdx].site[j].observedCode = cellExtd[c2].site[j].observedCode;
			    cellExtd[doubletIdx].site[j].trueNumExistingAlleles = cellExtd[c2].site[j].trueNumExistingAlleles;
			    cellExtd[doubletIdx].site[j].observedNumExistingAlleles = cellExtd[c2].site[j].observedNumExistingAlleles;
			    cellExtd[doubletIdx].site[j].maxNumExistingAllelesPerStrand = cellExtd[c2].site[j].maxNumExistingAllelesPerStrand;
			    cellExtd[doubletIdx].site[j].trueAlleles = cellExtd[c2].site[j].trueAlleles;
			    cellExtd[doubletIdx].site[j].observedAlleles = cellExtd[c2].site[j].observedAlleles;
			    cellExtd[doubletIdx].site[j].pastEvolutionaryEvents = cellExtd[c2].site[j].pastEvolutionaryEvents;
			    cellExtd[doubletIdx].site[j].hasADO = cellExtd[c2].site[j].hasADO;
			    cellExtd[doubletIdx].site[j].numADO = cellExtd[c2].site[j].numADO;
			    cellExtd[doubletIdx].site[j].hasAmplificationError = cellExtd[c2].site[j].hasAmplificationError;
			    cellExtd[doubletIdx].site[j].fixedAmplificationError = cellExtd[c2].site[j].fixedAmplificationError;
			    cellExtd[doubletIdx].site[j].amplificationError = cellExtd[c2].site[j].amplificationError;

			    SiteReadCountsExtd(&cellExtd[doubletIdx].site[j], sizeFactors[c2], j, probs, ngsEij, ampEij, rd);
				}

		  doubletIdx++;

			} // if we do doublet for c1
		} //c1
	}

/********************* PrintCATG  ************************/
/**	CATG output format with read counts
        http://nbviewer.jupyter.org/gist/dereneaton/d2fd4e70d29f5ee5d195/testing_cat.ipynb#View-the-.cat-results-files

        The first line has the number of samples and the number of sites.
        Following this is a transposed data matrix with the consensus base
   calls where each row is a site and each column a different sample. To the
   right of each site is a tab-separated list of counts of the four bases at
   that site in the order C,A,T,G.

        12 44500
        1A0	1B0	1C0	1D0	2E0	2F0	2G0	2H0	3I0
        YCCCCCCCCCCC	10,0,10,0	20,0,0,0	20,0,0,0 20,0,0, ...
        GGGGGGGGGGGG	0,0,0,20	0,0,0,20	0,0,0,20 0,0,0,20 ...
        AAAAAAAAAAAA	0,20,0,0	0,20,0,0	0,20,0,0 0,20,0,0 ...
        CCCCCCCCCCCC	20,0,0,0	20,0,0,0	20,0,0,0 20,0,0,0 ...
        AAAAAAAAAAAA	0,20,0,0	0,20,0,0	0, 20,0,0 0,20,0,0 ...
        ...
*/
void PrintCATG(FILE *fp) {
  int i, j, snv;

  fprintf(fp, "%d %d\n", numCells + 1, numSNVs);
  for (i = 0; i < numCells; i++)
    fprintf(fp, "%s%04d ", inCellName, i + 1);

  fprintf(fp, "%s ", outCellName);

  for (snv = 0; snv < numSNVs; snv++) {
    j = SNVsites[snv];
    fprintf(fp, "\n");
    for (i = 0; i < numCells + 1; i++)
      fprintf(
          fp, "%c",
          WhichIUPAC(cell[i].site[j].MLmatAllele, cell[i].site[j].MLpatAllele));

    if (cell[i].hasDoublet == NO) {
      for (i = 0; i < numCells + 1; i++)
        fprintf(fp, "\t%d,%d,%d,%d", cell[i].site[j].readCount[C],
                cell[i].site[j].readCount[A], cell[i].site[j].readCount[T],
                cell[i].site[j].readCount[G]);
    } else {
      for (i = 0; i < numCells + 1; i++)
        fprintf(fp, "\t%d,%d,%d,%d", cell[i].site[j].readCountDoublet[C],
                cell[i].site[j].readCountDoublet[A],
                cell[i].site[j].readCountDoublet[T],
                cell[i].site[j].readCountDoublet[G]);
    }
  }

  fprintf(fp, "\n");
}

/********************* PrintVCF  ************************/
/**	VCF output format with read counts
 */
void PrintVCF(FILE *fp) {
  int i, j, k, l, snv, a1, a2;
  int maternalAllele, paternalAllele, referenceAllele;
  int trueMaternalAllele, truePaternalAllele;
  int maternalAlleleDoublet, paternalAlleleDoublet;
  char *fmt = NULL;
  char *fmt_f3_1 = "%3.1f";
  char *fmt_comma_f3_1 = ",%3.1f";

  /* print file header */
  fprintf(fp, "##fileformat=VCFv4.3");
  fprintf(fp, "\n##filedate=");
  PrintDate(fp);
  fprintf(fp, "##source=CellCoal SNV simulation");
  fprintf(fp, "\n##contig=<ID=1>");
  fprintf(fp, "\n##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral "
              "allele\">");
  fprintf(fp, "\n##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of "
              "samples with data\">");
  fprintf(fp, "\n##INFO=<ID=AF,Number=R,Type=Float,Description=\"True "
              "alternate/s allele frequency (ignores doublets)\">");
  fprintf(fp, "\n##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic "
              "SNV flag\">");
  fprintf(fp, "\n##FILTER=<ID=s50,Description=\"Less than half of samples have "
              "data\">");
  fprintf(fp, "\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Called ML "
              "genotype\">");
  fprintf(
      fp,
      "\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
  fprintf(fp, "\n##FORMAT=<ID=RC,Number=4,Type=Integer,Description=\"Read "
              "counts for ACGT\">");
  fprintf(fp, "\n##FORMAT=<ID=G10,Number=G,Type=Float,Description=\"Scaled "
              "log10 genotype likelihoods for all 10 genotypes (order is AA AC "
              "AG AT CC CG CT GG GT TT)\">");
  fprintf(fp, "\n##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Scaled "
              "log10 genotype likelihoods\">");
  fprintf(fp, "\n##FORMAT=<ID=ML,Number=1,Type=String,Description=\"Maximum "
              "likelihood genotype\">");
  fprintf(fp, "\n##FORMAT=<ID=NG,Number=1,Type=String,Description=\"Genotype "
              "in the NGS library (considers ADO/DEL; ignores doublets)\">");
  fprintf(fp,
          "\n##FORMAT=<ID=DG,Number=1,Type=String,Description=\"Genotype in "
          "the NGS library for the second cell if there is a doublet\">");
  fprintf(fp, "\n##FORMAT=<ID=TG,Number=1,Type=String,Description=\"True "
              "genotype (considers DEL; ignores ADO/doublets\">");
  fprintf(fp, "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

  for (i = 0; i < numCells; i++)
    fprintf(fp, "\ttumcell%04d", i + 1);
  fprintf(fp, "\thealthycell");

  /* print only apparent SNVs */
  for (snv = 0; snv < numSNVs; snv++) {
    j = SNVsites[snv];

    /* VFC: CHROMOSOME */
    fprintf(fp, "\n%d", 1);

    /* VFC: POSITION */
    fprintf(fp, "\t%d", j + 1);

    /* VFC: ID */
    if (allSites[j].isSNV == NO)
      fprintf(fp, "\tinv%06d", j + 1);
    else
      fprintf(fp, "\tsnv%06d", j + 1);

    /* VFC: REFERENCE allele(s) => healthy root alleles */
    referenceAllele = allSites[j].referenceAllele;
    fprintf(fp, "\t%c", WhichNuc(referenceAllele));

    /* VFC: ALTERNATE allele(s) */
    fprintf(fp, "\t");
    if (allSites[j].numAltAlleles == 0)
      fprintf(fp, ".");
    else {
      fprintf(fp, "%c", WhichNuc(allSites[j].alternateAlleles[0]));
      for (k = 1; k < allSites[j].numAltAlleles; k++)
        fprintf(fp, ",%c", WhichNuc(allSites[j].alternateAlleles[k]));
    }

    /* VFC: QUALITY */ /* Qphred probability that a SNV exists at this site; I
                          will put missing info here */
    fprintf(fp, "\t.");

    /* VFC: FILTER */
    if (allSites[j].countACGT <
        allSites[j]
            .countN) /* if less than half of the samples have missing data */
      fprintf(fp, "\ts50");
    else
      fprintf(fp, "\tPASS");

    /* VFC: INFO: AA (ancestral allele) ; ancestral allele is the reference
     * alleles in this case */
    fprintf(fp, "\tAA=%c", WhichNuc(referenceAllele));

    /* VFC: INFO: NS (number of samples with data) */
    fprintf(fp, ";NS=%d", allSites[j].countCellswithData);

    /* VFC: INFO: AF (alternate allele frequencies) */
    fprintf(fp, ";AF=");
    if (allSites[j].countA > 0 && referenceAllele != A)
      fprintf(fp, "%4.3f,", (double)allSites[j].countA / allSites[j].countACGT);
    if (allSites[j].countC > 0 && referenceAllele != C)
      fprintf(fp, "%4.3f,", (double)allSites[j].countC / allSites[j].countACGT);
    if (allSites[j].countG > 0 && referenceAllele != G)
      fprintf(fp, "%4.3f,", (double)allSites[j].countG / allSites[j].countACGT);
    if (allSites[j].countT > 0 && referenceAllele != T)
      fprintf(fp, "%4.3f,", (double)allSites[j].countT / allSites[j].countACGT);
    fseek(fp, -1, SEEK_CUR); /* get rid of the last comma */

    /* VFC: INFO : SOMATIC */
    fprintf(fp, ";SOMATIC");

    /* VFC: FORMAT */
    fprintf(fp, "\tGT:DP:RC:G10:GL:ML:NG:DG:TG");

    for (i = 0; i < numCells + 1; i++) {
      maternalAllele = cell[i].site[j].MLmatAllele;
      paternalAllele = cell[i].site[j].MLpatAllele;
      trueMaternalAllele = cell[i].site[j].trueMaternalAllele;
      truePaternalAllele = cell[i].site[j].truePaternalAllele;
      maternalAlleleDoublet = cell[i].site[j].maternalAlleleDoublet;
      paternalAlleleDoublet = cell[i].site[j].paternalAlleleDoublet;

      /* VFC: FORMAT: GT (genotypes) */
      if (maternalAllele == MISSING)
        fprintf(fp, "\t.");
      else if (maternalAllele == referenceAllele)
        fprintf(fp, "\t0");
      else {
        fprintf(fp, "\t");
        for (l = 0; l < allSites[j].numAltAlleles; l++) {
          if (allSites[j].alternateAlleles[l] == maternalAllele) {
            fprintf(fp, "%d", l + 1);
            break;
          }
        }
      }
      fprintf(fp, "|");
      if (paternalAllele == MISSING)
        fprintf(fp, ".");
      else if (paternalAllele == referenceAllele)
        fprintf(fp, "0");
      else {
        for (l = 0; l < allSites[j].numAltAlleles; l++) {
          if (allSites[j].alternateAlleles[l] == paternalAllele) {
            fprintf(fp, "%d", l + 1);
            break;
          }
        }
      }

      if (cell[i].hasDoublet == NO) {
        /* VFC: FORMAT: DP (read depth) */
        fprintf(fp, ":%d", cell[i].site[j].numReads);

        /* VFC: FORMAT: RC (read counts ACGT) */
        fprintf(fp, ":%d,%d,%d,%d", cell[i].site[j].readCount[A],
                cell[i].site[j].readCount[C], cell[i].site[j].readCount[G],
                cell[i].site[j].readCount[T]);

        /* VFC: FORMAT: G10 (all 10 genotype likelihoods) */
        fprintf(fp, ":");
        fmt = fmt_f3_1;
        for (k = 0; k < 4; k++)
          for (l = k; l < 4; l++) {
            fprintf(fp, fmt, cell[i].site[j].scaledGenLike[k][l]);
            fmt = fmt_comma_f3_1;
          }
        // fseek(fp, -1, SEEK_CUR); 	/* rewind to get rid of the last comma
        // */

        /* VFC: FORMAT: GL (genotype likelihoods) */
        fprintf(fp, ":");
        fmt = fmt_f3_1;
        for (k = 0; k <= allSites[j].numAltAlleles; k++)
          for (l = 0; l <= k; l++) {
            if (k == 0)
              a1 = allSites[j].referenceAllele;
            else
              a1 = allSites[j].alternateAlleles[k - 1];

            if (l == 0)
              a2 = allSites[j].referenceAllele;
            else
              a2 = allSites[j].alternateAlleles[l - 1];

            //						fprintf (fp,
            //"%3.1f,",
            // cell[i].site[j].scaledGenLike[a1][a2]);
            fprintf(fp, fmt, cell[i].site[j].scaledGenLike[a1][a2]);
            fmt = fmt_comma_f3_1;

            // fprintf (stderr, "\n gl[%c][%c] = %3.1f,", WhichNuc(a1),
            // WhichNuc(a2), cell[i].site[j].scaledGenLike[a1][a2]);
          }
        //					fseek(fp, -1, SEEK_CUR);
        ///* rewind to get rid of the last comma */

        /* VFC: FORMAT: ML (maximum likelihood genotype) */
        if (cell[i].site[j].numReads > 0) {
          fprintf(fp, ":");
          fprintf(fp, "%c/%c", WhichNuc(cell[i].site[j].MLmatAllele),
                  WhichNuc(cell[i].site[j].MLpatAllele));
        } else
          fprintf(fp, ":./.");

        /* VFC: FORMAT: NG (NGS SNV genotype) */
        fprintf(fp, ":");
        fprintf(fp, "%c|%c", WhichNuc(maternalAllele),
                WhichNuc(paternalAllele));

        /* VFC: FORMAT: DG (true SNV second cell genotype, in doublet) */
        fprintf(fp, ":");
        fprintf(fp, "./.");

        /* VFC: FORMAT: TG (true SNV genotype) */
        fprintf(fp, ":");
        fprintf(fp, "%c|%c", WhichNuc(trueMaternalAllele),
                WhichNuc(truePaternalAllele));

      } else // if this is a doublet
      {
        /* VFC: FORMAT: DP (read depth) */
        fprintf(fp, ":%d", cell[i].site[j].numReadsDoublet);

        /* VFC: FORMAT: RC (read counts ACGT) */
        fprintf(fp, ":%d,%d,%d,%d", cell[i].site[j].readCountDoublet[A],
                cell[i].site[j].readCountDoublet[C],
                cell[i].site[j].readCountDoublet[G],
                cell[i].site[j].readCountDoublet[T]);

        /* VFC: FORMAT: G10 (genotype likelihoods) */
        fprintf(fp, ":");
        fmt = fmt_f3_1;
        for (k = 0; k < 4; k++)
          for (l = k; l < 4; l++) {
            fprintf(fp, fmt, cell[i].site[j].scaledGenLikeDoublet[k][l]);
            fmt = fmt_comma_f3_1;
          }
        //				fseek(fp, -1, SEEK_CUR); 	/*
        // rewind to get rid of the last comma */

        /* VFC: FORMAT: GL (genotype likelihoods) */
        fprintf(fp, ":");
        fmt = fmt_f3_1;
        for (k = 0; k <= allSites[j].numAltAlleles; k++)
          for (l = 0; l <= k; l++) {
            if (k == 0)
              a1 = allSites[j].referenceAllele;
            else
              a1 = allSites[j].alternateAlleles[k - 1];

            if (l == 0)
              a2 = allSites[j].referenceAllele;
            else
              a2 = allSites[j].alternateAlleles[l - 1];
            fprintf(fp, fmt, cell[i].site[j].scaledGenLikeDoublet[a1][a2]);
            fmt = fmt_comma_f3_1;
            // fprintf (stderr, "\n gl[%c][%c] = %3.1f,", WhichNuc(a1),
            // WhichNuc(a2), cell[i].site[j].scaledGenLikeDoublet[a1][a2]);
          }
        // fseek(fp, -1, SEEK_CUR); 	/* rewind to get rid of the last comma
        // */

        /* VFC: FORMAT: ML (maximum likelihood genotype) */
        if (cell[i].site[j].numReadsDoublet > 0) {
          fprintf(fp, ":");
          fprintf(fp, "%c/%c", WhichNuc(cell[i].site[j].MLmatAlleleDoublet),
                  WhichNuc(cell[i].site[j].MLpatAlleleDoublet));
        } else
          fprintf(fp, ":./.");

        /* VFC: FORMAT: NG (NGS SNV genotype) */
        fprintf(fp, ":");
        fprintf(fp, "%c|%c", WhichNuc(maternalAllele),
                WhichNuc(paternalAllele));

        /* VFC: FORMAT: DG (true NGS SNV second cell genotype, in doublet) */
        fprintf(fp, ":");
        fprintf(fp, "%c|%c", WhichNuc(maternalAlleleDoublet),
                WhichNuc(paternalAlleleDoublet));

        /* VFC: FORMAT: TG (true SNV genotype) */
        fprintf(fp, ":");
        fprintf(fp, "%c|%c", WhichNuc(trueMaternalAllele),
                WhichNuc(truePaternalAllele));
      }
    } // cell
  }   // snv

  fprintf(fp, "\n");
}

/********************* PrepareGlobalFiles **********************/
/* Open global files to output results */
void PrepareGlobalFiles(int argc, char **argv) {
  char File[MAX_NAME];

  /* contains the simulated tree in Newick format */
  if (doPrintTree == YES) {
    sprintf(File, "%s/%s", resultsDir, treeFile);
    if ((fpTrees = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  /* contains a list of times and branch lenghts for all nodes in the
   * simulated trees */
  if (doPrintTimes == YES) {
    sprintf(File, "%s/%s", resultsDir, timesFile);
    if ((fpTimes = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  if (doSimulateData == YES) {
    /* contains SNV genotypes for every cell */
    if (doPrintSNVgenotypes == YES) {
      sprintf(File, "%s/%s", resultsDir, SNVgenotypesFile);
      if ((fpSNVgenotypes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains haplotypes for variable sites for every cell */
    if (doPrintSNVhaplotypes == YES) {
      sprintf(File, "%s/%s", resultsDir, SNVhaplotypesFile);
      if ((fpSNVhaplotypes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains error-free haplotypes for variable sites for every cell */
    if (doPrintTrueHaplotypes == YES) {
      sprintf(File, "%s/%s", resultsDir, trueHaplotypesFile);
      if ((fpTrueHaplotypes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains all genotypes (variable or invariable) for every cell */
    if (doPrintFullGenotypes == YES) {
      sprintf(File, "%s/%s", resultsDir, fullGenotypesFile);
      if ((fpFullGenotypes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains haplotypes for all sites for every cell */
    if (doPrintFullHaplotypes == YES) {
      sprintf(File, "%s/%s", resultsDir, fullHaplotypesFile);
      if ((fpFullHaplotypes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains reads counts and genotype likelihoods for every SNV and cell
     */
    if (doNGS == YES) {
      sprintf(File, "%s/%s", resultsDir, VCFfile);
      if ((fpVCF = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains reads counts for every SNV and cell */
    if (doPrintCATG == YES) {
      sprintf(File, "%s/%s", resultsDir, CATGfile);
      if ((fpCATG = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }
  }

  if (doPrintSNVgenotypes == YES) {
    /*fprintf (fpSNVgenotypes, "%s - ",PROGRAM_NAME);
     PrintDate (fpSNVgenotypes);
     fprintf (fpSNVgenotypes, "SNV genotypes\n");
     PrintCommandLine (fpSNVgenotypes, argc, argv);*/
    fprintf(fpSNVgenotypes, "%d\n", numDataSets);
  }

  if (doPrintSNVhaplotypes == YES) {
    /*fprintf (fpSNVhaplotypes, "%s - ",PROGRAM_NAME);
    PrintDate (fpSNVhaplotypes);
    fprintf (fpSNVhaplotypes, "SNV haplotypes\n");
    PrintCommandLine (fpSNVhaplotypes, argc, argv);*/
    fprintf(fpSNVhaplotypes, "%d\n", numDataSets);
  }

  if (doPrintFullGenotypes == YES) {
    /*fprintf (fpFullGenotypes, "%s - ",PROGRAM_NAME);
    PrintDate (fpFullGenotypes);
    fprintf (fpFullGenotypes, "Full genotypes\n");
    PrintCommandLine (fpFullGenotypes, argc, argv);*/
    fprintf(fpFullGenotypes, "%d\n", numDataSets);
  }

  if (doPrintFullHaplotypes == YES) {
    /*fprintf (fpFullHaplotypes, "%s - ",PROGRAM_NAME);
    PrintDate (fpFullHaplotypes);
    fprintf (fpFullHaplotypes, "Full haplotypes\n");
    PrintCommandLine (fpFullHaplotypes, argc, argv);*/
    fprintf(fpFullHaplotypes, "%d\n", numDataSets);
  }

  if (doPrintTrueHaplotypes == YES) {
    /*fprintf (fpTrueHaplotypes, "%s - ",PROGRAM_NAME);
     PrintDate (fpTrueHaplotypes);
     fprintf (fpTrueHaplotypes, "True haplotypes\n");
     PrintCommandLine (fpTrueHaplotypes, argc, argv);*/
    fprintf(fpTrueHaplotypes, "%d\n", numDataSets);
  }

  if (doNGS == YES) {
    /*fprintf (fpVCF, "%s - ",PROGRAM_NAME);
    PrintDate (fpVCF);
    fprintf (fpVCF, "Read counts and genotype likelihoods\n");
    PrintCommandLine (fpVCF, argc, argv);*/
    fprintf(fpVCF, "%d\n", numDataSets);
  }

  if (doPrintCATG == YES) {
    /*fprintf (fpCATG, "%s - ",PROGRAM_NAME);
    PrintDate (fpCATG);
    fprintf (fpCATG, "Read counts\n");
    PrintCommandLine (fpCATG, argc, argv);*/
    fprintf(fpCATG, "%d\n", numDataSets);
  }
}

/* added code */
/********************* PrepareGlobalFilesExtd **********************/
/* Open global files to output results */
void PrepareGlobalFilesExtd(int argc, char **argv) {
  char File[MAX_NAME];

  /* contains a list of times and branch lengths for all nodes in the
   * simulated trees */
  if (doPrintTimes == YES) {
    sprintf(File, "%s/%s", resultsDir, timesFile);
    if ((fpTimes = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  /* contains all genotypes (variable or invariable) for every cell */
  if (doPrintFullGenotypes == YES) {
    sprintf(File, "%s/%s", resultsDir, fullGenotypesFile);
    if ((fpFullGenotypes = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  if (doPrintFullGenotypes == YES) {
    /*fprintf (fpFullGenotypes, "%s - ",PROGRAM_NAME);
    PrintDate (fpFullGenotypes);
    fprintf (fpFullGenotypes, "Full genotypes\n");
    PrintCommandLine (fpFullGenotypes, argc, argv);*/
    fprintf(fpFullGenotypes, "%d\n", numDataSets);
  }
}
/* end added code */

/********************* PrepareSeparateFiles **********************/
/* Open individual files to output results */

void PrepareSeparateFiles(int replicate) {
  /* contains the simulated tree in Newick format */
  if (doPrintTree == YES) {
    sprintf(File, "%s/%s", resultsDir, treeDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s.%04d", resultsDir, treeDir, treeFile,
            replicate + 1);
    if ((fpTrees = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  if (doPrintTimes == YES) {
    sprintf(File, "%s/%s", resultsDir, timesDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s.%04d", resultsDir, timesDir, timesFile,
            replicate + 1);
    if ((fpTimes = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  if (doSimulateData == YES) {
    /* contains haplotypes for variable sites for every cell */
    if (doPrintSNVhaplotypes == YES) {
      sprintf(File, "%s/%s", resultsDir, SNVhaplotypesDir);
      mkdir(File, S_IRWXU);
      sprintf(File, "%s/%s/%s.%04d", resultsDir, SNVhaplotypesDir,
              SNVhaplotypesFile, replicate + 1);
      if ((fpSNVhaplotypes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains reference haplotypes (before errors)  for every cell */
    if (doPrintTrueHaplotypes == YES) {
      sprintf(File, "%s/%s", resultsDir, trueHaplotypesDir);
      mkdir(File, S_IRWXU);
      sprintf(File, "%s/%s/%s.%04d", resultsDir, trueHaplotypesDir,
              trueHaplotypesFile, replicate + 1);
      if ((fpTrueHaplotypes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains all genotypes (variable or invariable) for every cell */
    if (doPrintFullGenotypes == YES) {
      sprintf(File, "%s/%s", resultsDir, fullGenotypesDir);
      mkdir(File, S_IRWXU);
      sprintf(File, "%s/%s/%s.%04d", resultsDir, fullGenotypesDir,
              fullGenotypesFile, replicate + 1);
      if ((fpFullGenotypes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains haplotypes for all sites for every cell */
    if (doPrintFullHaplotypes == YES) {
      sprintf(File, "%s/%s", resultsDir, fullHaplotypesDir);
      mkdir(File, S_IRWXU);
      sprintf(File, "%s/%s/%s.%04d", resultsDir, fullHaplotypesDir,
              fullHaplotypesFile, replicate + 1);
      if ((fpFullHaplotypes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains reads counts and log10 normalized genotype likelihoods for
     * every SNV and cell */
    if (doNGS == YES) {
      sprintf(File, "%s/%s", resultsDir, VCFdir);
      mkdir(File, S_IRWXU);
      sprintf(File, "%s/%s/%s.%04d", resultsDir, VCFdir, VCFfile,
              replicate + 1);
      if ((fpVCF = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }

    /* contains reads counts for every SNV and cell */
    if (doPrintCATG == YES) {
      sprintf(File, "%s/%s", resultsDir, CATGdir);
      mkdir(File, S_IRWXU);
      sprintf(File, "%s/%s/%s.%04d", resultsDir, CATGdir, CATGfile,
              replicate + 1);
      if ((fpCATG = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }
  }
}

/* added code */
/********************* PrepareSeparateFoldersAndFiles **********************/
void PrepareSeparateFoldersAndFiles(int replicate) {
  /* contains the simulated tree in Newick format */
  if (doPrintTree == YES && (altModel == FiniteMuInDel ||
                             altModel == FiniteMuDel || altModel == FiniteMu)) {
    sprintf(File, "%s/%s", resultsDir, treeDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s", resultsDir, treeDir, treeWithTrunkDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s/%s.%04d", resultsDir, treeDir, treeWithTrunkDir,
            treeFile, replicate + 1);
    if ((fpTrees = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }

    sprintf(File, "%s/%s", resultsDir, treeDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s", resultsDir, treeDir, treeWithoutTrunkDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s/%s.%04d", resultsDir, treeDir, treeWithoutTrunkDir,
            treeWithoutTrunkFile, replicate + 1);
    sprintf(treeWithoutTrunkFileCopy, "%s/%s/%s.%04d", treeDir,
            treeWithoutTrunkDir, treeWithoutTrunkFile, replicate + 1);
    if ((fpTreeWithoutTrunk = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }

    sprintf(File, "%s/%s", resultsDir, treeDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s", resultsDir, treeDir, treeWithHealthyTipDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s/%s.%04d", resultsDir, treeDir,
            treeWithHealthyTipDir, treeWithHealthyTipFile, replicate + 1);
    if ((fpTreeWithHealthyTip = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  sprintf(File, "%s/%s", resultsDir, cellNamesDir);
  mkdir(File, S_IRWXU);
  sprintf(File, "%s/%s/%s.%04d", resultsDir, cellNamesDir, cellNamesFile,
          replicate + 1);
  if ((fpCellNames = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  /* contains SNV genotypes */
  if (doPrintSNVgenotypes == YES) {
    sprintf(File, "%s/%s", resultsDir, SNVgenotypesDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s.%04d", resultsDir, SNVgenotypesDir,
            SNVgenotypesFile, replicate + 1);
    if ((fpSNVgenotypes = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  /* contains parallel mutations */
  sprintf(File, "%s/%s", resultsDir, parallelMutDir);
  mkdir(File, S_IRWXU);
  sprintf(File, "%s/%s/%s.%04d", resultsDir, parallelMutDir, parallelMutFile,
          replicate + 1);
  if ((fpParallelMut = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  /* contains true SNV genotypes (alphabetic, numeric, ternary matrix) */
  sprintf(File, "%s/%s", resultsDir, trueSNVDir);
  mkdir(File, S_IRWXU);

  sprintf(File, "%s/%s/%s.%04d", resultsDir, trueSNVDir, trueSNVGenotypesAlFile,
          replicate + 1);
  if ((fpTrueSNVGenotypesAl = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  sprintf(File, "%s/%s/%s.%04d", resultsDir, trueSNVDir, trueSNVGenotypesNrFile,
          replicate + 1);
  if ((fpTrueSNVGenotypesNr = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  sprintf(File, "%s/%s/%s.%04d", resultsDir, trueSNVDir,
          trueSNVTernaryMatrixFile, replicate + 1);
  if ((fpTrueSNVTernaryMatrix = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  sprintf(File, "%s/%s/%s.%04d", resultsDir, trueSNVDir, trueSNVCovFile,
          replicate + 1);
  if ((fpTrueSNVCov = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  sprintf(File, "%s/%s/%s.%04d", resultsDir, trueSNVDir, evolutionaryEventsFile,
          replicate + 1);
  if ((fpEvolutionaryEvents = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  sprintf(File, "%s/%s", resultsDir, SNVSitesNamesDir);
  mkdir(File, S_IRWXU);
  sprintf(File, "%s/%s/%s.%04d", resultsDir, SNVSitesNamesDir,
          SNVSitesNamesFile, replicate + 1);
  if ((fpSNVSitesNames = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  sprintf(File, "%s/%s", resultsDir, SNVAllelicSeqInfoDir);
  mkdir(File, S_IRWXU);
  sprintf(File, "%s/%s/%s.%04d", resultsDir, SNVAllelicSeqInfoDir,
          SNVAllelicSeqInfoFile, replicate + 1);
  if ((fpSNVAllelicSeqInfo = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  sprintf(File, "%s/%s", resultsDir, adoStatesDir);
  mkdir(File, S_IRWXU);
  sprintf(File, "%s/%s/%s.%04d", resultsDir, adoStatesDir, adoStatesFile,
          replicate + 1);
  if ((fpAdoStates = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  sprintf(File, "%s/%s", resultsDir, sizeFactorsDir);
  mkdir(File, S_IRWXU);
  sprintf(File, "%s/%s/%s.%04d", resultsDir, sizeFactorsDir, sizeFactorsFile,
          replicate + 1);
  if ((fpSizeFactors = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }
}

void PrepareRawReadSeqDataFiles(int replicate, int numTumCells)
{
  /* contains raw sequencing data */
  sprintf(File, "%s/%s", resultsDir, rawSeqDataDir);
  mkdir(File, S_IRWXU);
  sprintf(File, "%s/%s/%04d", resultsDir, rawSeqDataDir, replicate + 1);
  mkdir(File, S_IRWXU);

  sprintf(File, "%s/%s/%04d/%s.fasta", resultsDir, rawSeqDataDir, replicate + 1,
          rawSeqDataRefFile);
  sprintf(rawDataFiles[0].fileName, "%s.fasta", rawSeqDataRefFile);
  if ((rawDataFiles[0].fp = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  for (int i = 1; i < numTumCells + 2; i++) {
    if (healthyTip->label != i - 1) {
      sprintf(File, "%s/%s/%04d/%s.sam", resultsDir, rawSeqDataDir,
              replicate + 1, cellNamesExtd[i - 1]);
      sprintf(rawDataFiles[i].fileName, "%s.sam", cellNamesExtd[i - 1]);
    } else {
      sprintf(File, "%s/%s/%04d/healthycell.sam", resultsDir, rawSeqDataDir,
              replicate + 1);
      sprintf(rawDataFiles[i].fileName, "healthycell.sam");
    }

    if ((rawDataFiles[i].fp = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  // dummy bam file names
  if (replicate == 0) {
    sprintf(File, "%s/%s/%s", resultsDir, rawSeqDataDir, dummyBamDataDir);
    mkdir(File, S_IRWXU);

    FILE *fpTMP;

    for (int i = 0; i < numTumCells + 1; i++) {
      if (healthyTip->label == i) {
        sprintf(File, "%s/%s/%s/healthycell.bam", resultsDir, rawSeqDataDir,
                dummyBamDataDir);

        if ((fpTMP = fopen(File, "w")) == NULL) {
          fprintf(stderr, "Can't open \"%s\"\n", File);
          exit(-1);
        }
        fclose(fpTMP);
      } else {
        sprintf(File, "%s/%s/%s/%s.bam", resultsDir, rawSeqDataDir,
                dummyBamDataDir, cellNamesExtd[i]);
        if ((fpTMP = fopen(File, "w")) == NULL) {
          fprintf(stderr, "Can't open \"%s\"\n", File);
          exit(-1);
        }
        fclose(fpTMP);
      }
    }
  }

  // MPILEUP
  sprintf(File, "%s/%s", resultsDir, mpileupDataDir);
  mkdir(File, S_IRWXU);
  sprintf(File, "%s/%s/%04d", resultsDir, mpileupDataDir, replicate + 1);
  mkdir(File, S_IRWXU);

  sprintf(File, "%s/%s/%04d/%s.pileup", resultsDir, mpileupDataDir,
          replicate + 1, mpileupDataNMFile);
  if ((fpMpileupDataNM = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  sprintf(File, "%s/%s/%04d/%s.pileup", resultsDir, mpileupDataDir,
          replicate + 1, mpileupDataMFile);
  if ((fpMpileupDataM = fopen(File, "w")) == NULL) {
    fprintf(stderr, "Can't open \"%s\"\n", File);
    exit(-1);
  }

  if (meanDoubletRate > 0) {
    sprintf(File, "%s/%s/%04d/%s.pileup", resultsDir, mpileupDataDir,
            replicate + 1, mpileupDataNUFile);
    if ((fpMpileupDataNU = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }

    sprintf(File, "%s/%s/%04d/%s.pileup", resultsDir, mpileupDataDir,
            replicate + 1, mpileupDataUFile);
    if ((fpMpileupDataU = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }
}

/********************* PrepareSeparateFilesExtd **********************/
/* Open individual files to output results */
void PrepareSeparateFilesExtd(int replicate) {
  /* contains the simulated tree in Newick format */
  if (doPrintTree == YES &&
      !(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
        altModel == FiniteMu)) {
    sprintf(File, "%s/%s", resultsDir, treeDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s.%04d", resultsDir, treeDir, treeFile,
            replicate + 1);
    if ((fpTrees = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  if (doPrintTimes == YES) {
    if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu) ||
        ((altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu) &&
         (fixedTree == NO || (fixedTree == YES && replicate == 0)))) {
      sprintf(File, "%s/%s", resultsDir, timesDir);
      mkdir(File, S_IRWXU);
      sprintf(File, "%s/%s/%s.%04d", resultsDir, timesDir, timesFile,
              replicate + 1);
      if ((fpTimes = fopen(File, "w")) == NULL) {
        fprintf(stderr, "Can't open \"%s\"\n", File);
        exit(-1);
      }
    }
  }

  /* contains all genotypes (variable or invariable) for every cell */
  if (doPrintFullGenotypes == YES) {
    sprintf(File, "%s/%s", resultsDir, fullGenotypesDir);
    mkdir(File, S_IRWXU);
    sprintf(File, "%s/%s/%s.%04d", resultsDir, fullGenotypesDir,
            fullGenotypesFile, replicate + 1);
    if ((fpFullGenotypes = fopen(File, "w")) == NULL) {
      fprintf(stderr, "Can't open \"%s\"\n", File);
      exit(-1);
    }
  }

  PrepareSeparateFoldersAndFiles(replicate);
}
/* end added code */

/**************** Print Trees ***************/
/*  Print trees to treefile in Newick format */
void PrintTree(TreeNode *treeRoot, FILE *fp) {
  WriteTree(treeRoot, fp);
  fprintf(fp, ";\n");
}

/******************* WriteTree ****************/
/* Writes a given tree starting on a particular node (usually the root, but
 * not necessarily) */
void WriteTree(TreeNode *p, FILE *fp) {
  if (p != NULL) {
    if (p->isHealthyTip == YES) {
      int tmpLabel = Label(p);
      sprintf(cellNamesExtd[tmpLabel - 1], "healthycell");
      if (p->branchLength != 0.0)
        fprintf(fp, "healthycell:%8.18f", p->branchLength);
      else
        fprintf(fp, "healthycell:0.0");
    } else if (p->left == NULL && p->right == NULL) {
      if (doUserTree == NO) {
        int tmpLabel = Label(p);
        sprintf(cellNamesExtd[tmpLabel - 1], "tumcell%04d", tmpLabel);
        fprintf(fp, "%s:", cellNamesExtd[tmpLabel - 1]);
        if (p->branchLength != 0.0)
          fprintf(fp, "%8.18f", p->branchLength);
        else
          fprintf(fp, "0.0");
      } else {
        fprintf(fp, "%s:", p->name);
        if (p->branchLength != 0.0)
          fprintf(fp, "%8.18f", p->branchLength);
        else
          fprintf(fp, "0.0");
      }
    } else {
      fprintf(fp, "(");
      WriteTree(p->left, fp);
      if (!(p->isHealthyRoot == YES && doPrintHealthyTip == NO)) {
        fprintf(fp, ",");
        WriteTree(p->right, fp);
      }
      if (p->branchLength != 0.0)
        fprintf(fp, "):%8.18f", p->branchLength);
      else
        fprintf(fp, "):0.0");
    }
  }
}

/********************* PrintTimes **********************/
/*	Prints to a given file a detailed description of
        the tree as lost of nodes, times, branch lengths.
 */

void PrintTimes(int listPosition, FILE *fp) {
  fprintf(fp, "\n-------------------- Nodes -------------------");
  fprintf(fp, "\n       class  label  index  (left right anc) |         time   "
              "  time length    branch length");
  fprintf(fp, "\n--------------------------------------------------------------"
              "------------------------------\n");
  ListTimes(listPosition, fp);
}

/********************** ListTimes ************************/
/*	Writes a given tree description from ListTimes
        It does not list superfluous nodes
 */

void ListTimes(int position, FILE *fp) {
  TreeNode *p;
  do {
    p = treeNodes + position;

    if (p->isHealthyTip == YES)
      fprintf(fp,
              "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       "
              "%10.4lf\n",
              "outgroupTip", Label(p), Index(p), Index(p->left),
              Index(p->right), Index(p->anc), p->time, p->length,
              p->branchLength);
    else if (p->anc != NULL && p->left != NULL && p->right != NULL &&
             p->anc->anc == NULL)
      fprintf(fp,
              "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       "
              "%10.4lf\n",
              "ingroupMRCA", Label(p), Index(p), Index(p->left),
              Index(p->right), Index(p->anc), p->time, p->length,
              p->branchLength);
    else if (p->anc != NULL && p->left != NULL && p->right != NULL)
      fprintf(fp,
              "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       "
              "%10.4lf\n",
              "internal", Label(p), Index(p), Index(p->left), Index(p->right),
              Index(p->anc), p->time, p->length, p->branchLength);
    else if (p->anc != NULL && p->left == NULL && p->right == NULL)
      fprintf(fp,
              "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       "
              "%10.4lf\n",
              "tip", Label(p), Index(p), Index(p->left), Index(p->right),
              Index(p->anc), p->time, p->length, p->branchLength);
    else if (p->anc == NULL && p->left != NULL && p->right != NULL)
      fprintf(fp,
              "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       "
              "%10.4lf\n",
              "outgroupRoot", Label(p), Index(p), Index(p->left),
              Index(p->right), Index(p->anc), p->time, p->length,
              p->branchLength);

    position++;

    if (position > 500000) /* just in case the tree is weird */
    {
      fprintf(stderr, "\n\nERROR: Problems when listing times in the tree\n\n");
      exit(-1);
    }
  } while (p->isHealthyTip == NO);
}

int GetMaxNumExistingAllelesPerStrand() {
  if (altModel == FiniteMuInDel)
    return 2;
  else if (altModel == FiniteMuDel || altModel == FiniteMu)
    return 1;
  return 1;
}

/************************ PrintSiteInfo ***********************/
/* Prints info for the site struct */

static void PrintSiteInfo(FILE *fp, int i) {
  int j;
  fprintf(fp, "\n\nSite %d", i + 1);
  fprintf(fp, "\n isSNV = %d", allSites[i].isSNV);
  fprintf(fp, "\n isSNP = %d", allSites[i].isSNP);
  fprintf(fp, "\n numMutations = %d", allSites[i].numMutations);
  fprintf(fp, "\n  numMutationsMaternal = %d",
          allSites[i].numMutationsMaternal);
  fprintf(fp, "\n  numMutationsPaternal = %d",
          allSites[i].numMutationsPaternal);
  fprintf(fp, "\n numDeletions = %d", allSites[i].numDeletions);
  fprintf(fp, "\n  numDeletionsMaternal = %d",
          allSites[i].numDeletionsMaternal);
  fprintf(fp, "\n  numDeletionsPaternal = %d",
          allSites[i].numDeletionsPaternal);
  fprintf(fp, "\n numCNLOH = %d", allSites[i].numCNLOH);
  fprintf(fp, "\n  numCNLOHmaternal = %d", allSites[i].numCNLOHmaternal);
  fprintf(fp, "\n  numCNLOHpaternal = %d", allSites[i].numCNLOHpaternal);
  fprintf(fp, "\n hasADO = %d", allSites[i].hasADO);
  fprintf(fp, "\n hasGenotypingError = %d", allSites[i].hasGenotypeError);
  fprintf(fp, "\n referenceAllele = %c", WhichNuc(allSites[i].referenceAllele));
  fprintf(fp, "\n numAltAlleles (in reads) = %d   (",
          allSites[i].numAltAlleles);
  for (j = 0; j < allSites[i].numAltAlleles; j++)
    fprintf(fp, " %c", WhichNuc(allSites[i].alternateAlleles[j]));
  fprintf(fp, " )");
  fprintf(fp, "\n alleleA = %d", allSites[i].countA);
  fprintf(fp, "\n alleleC = %d", allSites[i].countC);
  fprintf(fp, "\n alleleG = %d", allSites[i].countG);
  fprintf(fp, "\n alleleT = %d", allSites[i].countT);
  fprintf(fp, "\n alleleN = %d", allSites[i].countN);
  fprintf(fp, "\n countDropped = %d", allSites[i].countDropped);
  fprintf(fp, "\n readCountA = %d", allSites[i].readCountA);
  fprintf(fp, "\n readCountC = %d", allSites[i].readCountC);
  fprintf(fp, "\n readCountG = %d", allSites[i].readCountG);
  fprintf(fp, "\n readCountT = %d", allSites[i].readCountT);
  fprintf(fp, "\n rateMultiplier = %f", allSites[i].rateMultiplier);
}

/* added code */
/************************ PrintSiteInfoExtd ***********************/
/* Prints info for the site struct */
void PrintSiteInfoExtd(FILE *fp, int i) {
  int j;
  fprintf(fp, "\n\nSite %d", i + 1);
  fprintf(fp, "\n referenceAllele = %c",
          WhichNucExtd(allSites[i].referenceAllele, MISSING_CHAR));
  fprintf(fp, "\n isSNV = %d", allSites[i].isSNV);
  fprintf(fp, "\n isVariant = %d", allSites[i].isVariant);
  fprintf(fp, "\n MLVariantNucleotide = %d", allSites[i].MLVariantNucleotide);

  fprintf(fp, "\n numMutations = %d", allSites[i].numMutations);
  fprintf(fp, "\n isDeleted = %d", allSites[i].isDeleted);
  fprintf(fp, "\n numDeletions = %d", allSites[i].numDeletions);
  fprintf(fp, "\n isInserted = %d", allSites[i].isInserted);
  fprintf(fp, "\n numInsertions = %d", allSites[i].numInsertions);

  fprintf(fp, "\n hasADO = %d", allSites[i].hasADO);

  fprintf(fp, "\n numAltAlleles (in reads) = %d   (",
          allSites[i].numAltAlleles);
  for (j = 0; j < allSites[i].numAltAlleles; j++)
    fprintf(fp, " %c",
            WhichNucExtd(allSites[i].alternateAlleles[j], MISSING_CHAR));
  fprintf(fp, " )");

  fprintf(fp, "\n alleleA = %d", allSites[i].countA);
  fprintf(fp, "\n alleleC = %d", allSites[i].countC);
  fprintf(fp, "\n alleleG = %d", allSites[i].countG);
  fprintf(fp, "\n alleleT = %d", allSites[i].countT);
  fprintf(fp, "\n alleleN = %d", allSites[i].countN);

  fprintf(fp, "\n readCountA = %d", allSites[i].readCountA);
  fprintf(fp, "\n readCountC = %d", allSites[i].readCountC);
  fprintf(fp, "\n readCountG = %d", allSites[i].readCountG);
  fprintf(fp, "\n readCountT = %d", allSites[i].readCountT);

  fprintf(fp, "\n rateMultiplier = %f", allSites[i].rateMultiplier);
}
/* end added code */

/************************ PrintSNVGenotypes ***********************/
/* Prints observed/ML genotypes at variable sites (SNVs) to a file */
static void PrintSNVGenotypes(FILE *fp) {
  int i, j, k;

  if (doPrintAncestors == YES && doNGS == NO)
    fprintf(fp, "%d %d\n", numCells + 3, numSNVs);
  else
    fprintf(fp, "%d %d\n", numCells + 1, numSNVs);

  /* site information */
  fprintf(fp, "%d", SNVsites[0] + 1);
  for (i = 1; i < numSNVs; i++)
    fprintf(fp, " %d", SNVsites[i] + 1);
  fprintf(fp, "\n");
  k = 0;

  if (alphabet == DNA) {
    for (i = 0; i < numCells + 1; i++) {
      if (i == numCells)
        fprintf(fp, "%s ", outCellName);
      else {
        if (doUserTree == NO)
          fprintf(fp, "%s%04d ", inCellName, i + 1);
        else
          fprintf(fp, "%-*s", stringPrecision, cellNames[i]);
      }
      for (j = 0; j < numSNVs; j++) {
        if (doNGS == NO)
          fprintf(fp, " %c%c", WhichNuc(data[MATERNAL][i][SNVsites[j]]),
                  WhichNuc(data[PATERNAL][i][SNVsites[j]]));
        else if (cell[i].hasDoublet == NO)
          fprintf(fp, " %c%c", WhichNuc(cell[i].site[SNVsites[j]].MLmatAllele),
                  WhichNuc(cell[i].site[SNVsites[j]].MLpatAllele));
        else
          fprintf(fp, " %c%c",
                  WhichNuc(cell[i].site[SNVsites[j]].MLmatAlleleDoublet),
                  WhichNuc(cell[i].site[SNVsites[j]].MLpatAlleleDoublet));
      }
      fprintf(fp, "\n");
    }

    if (doPrintAncestors == YES & doNGS == NO) {
      fprintf(fp, "%s ", outRootCellName);
      for (j = 0; j < numSNVs; j++)
        fprintf(fp, " %c%c",
                WhichNuc(data[MATERNAL][OUTGROUP_ROOT][SNVsites[j]]),
                WhichNuc(data[PATERNAL][OUTGROUP_ROOT][SNVsites[j]]));
      fprintf(fp, "\n%s ", inRootCellName);
      for (j = 0; j < numSNVs; j++)
        fprintf(fp, " %c%c",
                WhichNuc(data[MATERNAL][INGROUP_ROOT][SNVsites[j]]),
                WhichNuc(data[PATERNAL][INGROUP_ROOT][SNVsites[j]]));
      fprintf(fp, "\n");
    }
  } else // binary data
  {
    for (i = 0; i < numCells + 1; i++) {
      if (i == numCells)
        fprintf(fp, "%s ", outCellName);
      else {
        if (doUserTree == NO)
          fprintf(fp, "%s%04d ", inCellName, i + 1);
        else
          fprintf(fp, "%-*s", stringPrecision, cellNames[i]);
      }
      for (j = 0; j < numSNVs; j++)
        fprintf(fp, " %c%c", WhichMut(data[MATERNAL][i][SNVsites[j]]),
                WhichMut(data[PATERNAL][i][SNVsites[j]]));
      fprintf(fp, "\n");
    }

    if (doPrintAncestors == YES) {
      fprintf(fp, "%s ", outRootCellName);
      for (j = 0; j < numSNVs; j++)
        fprintf(fp, " %d%d", data[MATERNAL][OUTGROUP_ROOT][SNVsites[j]],
                data[PATERNAL][OUTGROUP_ROOT][SNVsites[j]]);
      fprintf(fp, "\n%s ", inRootCellName);
      for (j = 0; j < numSNVs; j++)
        fprintf(fp, " %d%d", data[MATERNAL][INGROUP_ROOT][SNVsites[j]],
                data[PATERNAL][INGROUP_ROOT][SNVsites[j]]);
      fprintf(fp, "\n");
    }
  }
}

/* added code */
/************************ PrintTumCellNames ***********************/
void PrintTumCellNames(FILE *fp, char col_sep, int includeNormalCell) {
  int doublet_with, doublet_as;

  fprintf(fp, "cell%cdoublet_with%cdoublet_as\n", col_sep, col_sep);

  for (int i = 0; i < numCells + 1; i++) {
    doublet_with = cellExtd[i].mixeTumorCellIdx;
    doublet_as = cellExtd[i].mixedTumorCellNameIdx;

    if (includeNormalCell > 0 || healthyTip->label != i)
      fprintf(fp, "%s%c%s%c%s\n", cellNamesExtd[i], col_sep,
        doublet_with == -1 ? "" : cellNamesExtd[doublet_with], col_sep,
        doublet_as == -1 ? "" : cellNamesExtd[doublet_as]);
  }
} // PrintCellNames

/************************ PrintSNVSitesNames ***********************/
void PrintSNVSitesNames(FILE *fp) {
  fprintf(fp, "#chromosome\tsite\tSNVs(0: NO; 1: YES)\tdeletions(0: NO; 1: YES)\tinsertions(0: NO; 1: YES)\n");
  for (int i = 0; i < numSNVs; i++) {
    fprintf(fp, "chr1\t%d\t%d\t%d\t%d",
            SNVsites[i] + 1, mutatedSites[i], deletedSites[i], insertedSites[i]);

    if (i < numSNVs - 1) {
      fprintf(fp, "\n");
    }
  }
}

/************************ PrintSNVAllelicInformation ***********************/
void PrintSNVAllelicInformation(FILE *fp) {
  fprintf(fp, "#coverage raw_variance\n");
  for (int i = 0; i < numSNVs; i++) {
    fprintf(fp, "%f %f", allSites[SNVsites[i]].allelicCov,
            allSites[SNVsites[i]].allelicRawVar);

    if (i < numSNVs - 1) {
      fprintf(fp, "\n");
    }
  }
}

/************************ PrintAdoStates ***********************/
void PrintAdoStates(FILE *fp) {
  fprintf(fp, "#The number of ADOs occurred.\n#CHROM\tPOS\n");

  for (int j = 0; j < numSites; j++) {

    fprintf(fp, "chr1\t%d", j + 1);

    for (int i = 0; i < numCells + 1; i++)
      fprintf(fp, "\t%d", cellExtd[i].site[j].numADO);

    if (j < numSites - 1) {
      fprintf(fp, "\n");
    }
  }
}

/************************ PrintParallelMut ***********************/
void PrintParallelMut(FILE *fp) {
  fprintf(fp, "#site_index\tcells (in the same order as they are listed in "
              "'cell_names')\n");
  fprintf(fp, "#0: no parallel mutation. 1: parallel mutation\n");

  for (int i = 0; i < numThreads; i++) {
    for (int j = 0; j < threadStr[i].numSubSites; j++) {
      if (threadStr[i].hasParallelMutations[j] == YES) {
        fprintf(fp, "chr1_%d", threadStr[i].start + j);

        for (int k = 0; k < numCells + 1; k++) {
          if (healthyTip->label != k) {
            fprintf(fp, " %d", threadStr[i].parallelMutations[j][k]);
          }
        }

        fprintf(fp, "\n");
      }
    }
  }
}

/************************ PrintSizeFactors ***********************/
void PrintSizeFactors(FILE *fp) {
  for (int i = 0; i < numCells; i++)
    fprintf(fp, "%f\n", sizeFactors[i]);
}

/************************ PrintMPileUp ***********************/
void PrintMPileUp(FILE *fpNM, FILE *fpNU, FILE *fpM, FILE *fpU) {
  int max_doublet_idx = numCells;

  for (int j = 0; j < numSites; j++) {

    // with healthy cell; merged doublets
    if (fpNM != NULL)
      fprintf(fpNM, "chr1\t%d\t%c", j + 1,
              WhichNucExtd(allSites[j].referenceAllele, MISSING_CHAR));

    // with healthy cell; unmerged doublets
    if (fpNU != NULL)
      fprintf(fpNU, "chr1\t%d\t%c", j + 1,
              WhichNucExtd(allSites[j].referenceAllele, MISSING_CHAR));

    // without healthy cell; merged doublets
    if (fpM != NULL)
      fprintf(fpM, "chr1\t%d\t%c", j + 1,
              WhichNucExtd(allSites[j].referenceAllele, MISSING_CHAR));

    // without healthy cell; unmerged doublets
    if (fpU != NULL)
      fprintf(fpU, "chr1\t%d\t%c", j + 1,
              WhichNucExtd(allSites[j].referenceAllele, MISSING_CHAR));

    for (int i = 0; i < 2 * numCells + 1; i++) {
      if (j == 0 && i < numCells + 1 && cellExtd[i].isDoublet == YES)
        max_doublet_idx++;

      if (i > max_doublet_idx) break;

      const int cov = cellExtd[i].site[j].numReads;
      char reads[cov == 0 ? 2 : cov + 1];
      char quals[cov == 0 ? 2 : cov + 1];
      FillMpileupEntry(allSites[j].referenceAllele, cellExtd[i].site + j, NULL, cov, reads, quals);

      if (fpNU != NULL)
        fprintf(fpNU, "\t%d\t%s\t%s", cov, reads, quals);
      if (fpU != NULL && healthyTip->label != i)
        fprintf(fpU, "\t%d\t%s\t%s", cov, reads, quals);

      if (i < numCells + 1)
      {
        int merged_cov = cov;
        char *merged_reads = NULL, *merged_quals = NULL;

        if (cellExtd[i].isDoublet == YES)
        {
          merged_cov += cellExtd[cellExtd[i].mixedTumorCellNameIdx].site[j].numReads;

          merged_reads = (char *)malloc((merged_cov == 0 ? 3 : merged_cov + 2) * sizeof(char));
          merged_quals = (char *)malloc((merged_cov == 0 ? 3 : merged_cov + 2) * sizeof(char));

          FillMpileupEntry(allSites[j].referenceAllele,
            cellExtd[i].site + j,
            cellExtd[cellExtd[i].mixedTumorCellNameIdx].site + j,
            cov, merged_reads, merged_quals);
        } else
        {
          merged_reads = reads;
          merged_quals = quals;
        }

        if (fpNM != NULL)
          fprintf(fpNM, "\t%d\t%s\t%s", merged_cov, merged_reads, merged_quals);
        if (fpM != NULL && healthyTip->label != i)
          fprintf(fpM, "\t%d\t%s\t%s", merged_cov, merged_reads, merged_quals);

        if (cellExtd[i].isDoublet == YES)
        {
          free(merged_reads);
          merged_reads = NULL;
          free(merged_quals);
          merged_quals = NULL;
        }
      }

    }

    if (j < numSites - 1) {
      if (fpNM != NULL)
        fprintf(fpNM, "\n");
      if (fpNU != NULL)
        fprintf(fpNU, "\n");
      if (fpM != NULL)
        fprintf(fpM, "\n");
      if (fpU != NULL)
        fprintf(fpU, "\n");
    }
  }
}

/************************ PrintRawData ***********************/
void PrintRawData(RawDataFile *files) {
  // print reference genome
  fprintf(stderr, "\n>>> Printing reference sequence to %s...\n",
          files[0].fileName);
  fprintf(files[0].fp, ">chr1\n");
  for (int i = 0; i < numSites; i++) {
    if (i % 80 == 0 && i != 0 && i != numSites - 1)
      fprintf(files[0].fp, "\n");

    fprintf(files[0].fp, "%c",
            WhichNucExtd(allSites[i].referenceAllele, MISSING_CHAR));
  }
  fprintf(stderr, ">>> Done.\n\n");

  // print SAM
  fprintf(stderr, ">>> Printing reads of all cells...\n");
  PrintRawDataSam();
}

/************************ PrintRawDataFASTQ ***********************/
void PrintRawDataSam() {
  pthread_attr_t attr;
  double t_start, t_end;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  t_start = omp_get_wtime();
  for (long i = 0; i < numThreads; i++) {
    pthread_create(&threads[i], &attr, ThreadedPrintRawDataSam,
                   (void *)(rawDataPrinters + i));
  }
  pthread_attr_destroy(&attr);

  for (long i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }
  t_end = omp_get_wtime();

  fprintf(stderr, "Multi-threading finished in %f seconds.\n", t_end - t_start);
}

/************************ ThreadedPrintRawDataSam ***********************/
void *ThreadedPrintRawDataSam(void *arg) {
  ThreadedRawDataPrinter *printer;
  double t_start, t_end;
  const int maxNumExistingAllelesPerStrand =
      GetMaxNumExistingAllelesPerStrand();

  printer = (ThreadedRawDataPrinter *)arg;
  fprintf(stderr, "Thread %ld is starting...\n", printer->id);

  char readEnds[maxNumExistingAllelesPerStrand * readLengths[0] + 1],
      readInternal[maxNumExistingAllelesPerStrand * readLengths[1] + 1],
      cigarStrEnds[maxNumExistingAllelesPerStrand * readLengths[0] + 1],
      cigarStrInternal[maxNumExistingAllelesPerStrand * readLengths[1] + 1],
      qualEnds[maxNumExistingAllelesPerStrand * readLengths[0] + 1],
      qualInternal[maxNumExistingAllelesPerStrand * readLengths[1] + 1];

  int cigarCntEnds[maxNumExistingAllelesPerStrand * readLengths[0]],
      cigarCntInternal[maxNumExistingAllelesPerStrand * readLengths[1]];

  // coverage, A, C, G, T
  int cellReadsBlock[readLengths[1]][5];

  t_start = omp_get_wtime();
  for (int i = 0; i < printer->numCells; i++) {
    const unsigned cellIndex = printer->cellIndices[i];
    char cellName[MAX_NAME];
    FILE *fp = printer->files[i].fp;

    GetCellName(printer->files[i].fileName, cellName, MAX_NAME);

    fprintf(stderr, ">>> [Thread %ld] Printing cell %d to %s...\n", printer->id,
            cellIndex + 1, printer->files[i].fileName);

    // print header of a SAM file
    fprintf(fp, "@HD\tVN:1.6\tSO:coordinate\n");
    fprintf(fp, "@SQ\tSN:chr1\tLN:%d\n", numSites);
    fprintf(fp, "@RG\tID:%s\tSM:%s\tPL:Illumina\n", cellName, cellName);

    for (int j = 0; j < numReadBlocks; j++) {
      char *reads, *cigarStr, *quals;
      int *cigarCnt, *counts;

      if (j == 0 || j == numReadBlocks - 1) {
        reads = readEnds;
        cigarStr = cigarStrEnds;
        cigarCnt = cigarCntEnds;
        quals = qualEnds;
      } else {
        reads = readInternal;
        cigarStr = cigarStrInternal;
        cigarCnt = cigarCntInternal;
        quals = qualInternal;
      }

      int initializeReads = YES, iteration, cigarStrIdx, startPos, readPos,
          nucNumToWrite;
      do {
        // reset the number of iterations
        iteration = INT_MIN;

        // reset cigar string counter
        cigarStrIdx = -1;
        memset(cigarCntEnds, 0,
               maxNumExistingAllelesPerStrand * readLengths[0] * sizeof(int));
        memset(cigarCntInternal, 0,
               maxNumExistingAllelesPerStrand * readLengths[1] * sizeof(int));

        startPos = -1;
        readPos = 0;
        for (int k = startSiteBlocks[j]; k < startSiteBlocks[j + 1]; k++) {
          if (initializeReads == YES) {
            cellReadsBlock[k - startSiteBlocks[j]][0] =
                cellExtd[cellIndex].site[k].numReads;
            cellReadsBlock[k - startSiteBlocks[j]][1] =
                cellExtd[cellIndex].site[k].readCount[0];
            cellReadsBlock[k - startSiteBlocks[j]][2] =
                cellExtd[cellIndex].site[k].readCount[1];
            cellReadsBlock[k - startSiteBlocks[j]][3] =
                cellExtd[cellIndex].site[k].readCount[2];
            cellReadsBlock[k - startSiteBlocks[j]][4] =
                cellExtd[cellIndex].site[k].readCount[3];
          }

          nucNumToWrite =
              cellExtd[cellIndex].site[k].maxNumExistingAllelesPerStrand;
          if (cellExtd[cellIndex].site[k].maxNumExistingAllelesPerStrand > 1 &&
              cellExtd[cellIndex].site[k].observedNumExistingAlleles -
                      cellExtd[cellIndex]
                          .site[k]
                          .maxNumExistingAllelesPerStrand >
                  0 &&
              RandUniform(printer->rd) < 0.5)
            nucNumToWrite =
                cellExtd[cellIndex].site[k].observedNumExistingAlleles -
                cellExtd[cellIndex].site[k].maxNumExistingAllelesPerStrand;

          int nucIdx = 0;
          do {
            char tmp = GetANucleotide(cellReadsBlock[k - startSiteBlocks[j]],
                                      printer->rd);

            if (tmp != 'N') {
              if (iteration <= cellReadsBlock[k - startSiteBlocks[j]][0] ||
                  nucIdx > 0)
                iteration = cellReadsBlock[k - startSiteBlocks[j]][0];

              reads[readPos] = tmp;
              quals[readPos] = '~';
              readPos++;

              if (startPos < 0)
                startPos = k;
            }

            if (tmp != 'N' && nucIdx > 0) {
              cigarStrIdx++;
              cigarStr[cigarStrIdx] = 'I';
            } else if (cigarStrIdx == -1 ||
                       (tmp == 'N' && cigarStr[cigarStrIdx] != 'D') ||
                       (tmp != 'N' && cigarStr[cigarStrIdx] != 'M')) {
              cigarStrIdx++;
              cigarStr[cigarStrIdx] = tmp == 'N' ? 'D' : 'M';
            }

            cigarCnt[cigarStrIdx]++;
            nucIdx++;
          } while (nucIdx < nucNumToWrite);
        }

        reads[readPos] = '\0';
        quals[readPos] = '\0';

        // write non-empty read to the file
        if (readPos > 0) {
          // 1: QNAME
          fprintf(fp, "*\t");

          // 2: FLAG
          fprintf(fp, "0\t");

          // 3: RNAME
          fprintf(fp, "chr1\t");

          // 4: POS (1-based)
          fprintf(fp, "%d\t", startPos + 1);

          // 5: MAPQ (assuming perfect mapping)
          fprintf(fp, "254\t");

          // 6: CIGAR
          for (int ii = 0; ii < cigarStrIdx + 1; ii++) {
            if ((ii == 0 || ii == cigarStrIdx) && cigarStr[ii] == 'D')
              continue;

            fprintf(fp, "%d%c", cigarCnt[ii], cigarStr[ii]);
          }
          fprintf(fp, "\t");

          // 7: RNEXT
          fprintf(fp, "*\t");

          // 8: PNEXT
          fprintf(fp, "0\t");

          // 9: TLEN
          fprintf(fp, "0\t");

          // 10: SEQ
          fprintf(fp, "%s\t", reads);

          // 11: QUAL (assuming highest quality of base read)
          fprintf(fp, "%s\t", quals);

          // 12 (optional tags): RG
          fprintf(fp, "RG:Z:%s", cellName);

          initializeReads = NO;

          if (j < numReadBlocks - 1 || iteration != 0)
            fprintf(fp, "\n");
        }
      } while (iteration > 0);
    }
  }
  t_end = omp_get_wtime();

  fprintf(stderr, "Thread %ld finished in %f seconds.\n", printer->id,
          t_end - t_start);

  pthread_exit(NULL);
}

char GetANucleotide(int *readCounts, RandomDevice *rd) {
  if (readCounts[0] < 0) {
    fprintf(stderr, "Error: read depth = %d < 0", readCounts[0]);
    exit(1);
  } else if (readCounts[0] == 0)
    return 'N';
  else {
    double probs[4];
    double cumProbs[5];

    cumProbs[0] = 0.0;
    cumProbs[4] = 1.0;
    for (int i = 0; i < 4; i++) {
      probs[i] = ((double)readCounts[i + 1]) / readCounts[0];
      cumProbs[i + 1] = cumProbs[i] + probs[i];
    }

    double value = RandUniform(rd);
    int category;
    for (int i = 0; i < 4; i++) {
      if (value >= cumProbs[i] && value < cumProbs[i + 1])
        category = i;
    }

    switch (category) {
    case 0:
      readCounts[0]--;
      readCounts[1]--;
      return 'A';
    case 1:
      readCounts[0]--;
      readCounts[2]--;
      return 'C';
    case 2:
      readCounts[0]--;
      readCounts[3]--;
      return 'G';
    case 3:
      readCounts[0]--;
      readCounts[4]--;
      return 'T';
    default:
      return 'N';
    }
  }
}

/************************ PrintTrueSNVCov ***********************/
void PrintTrueSNVCov(FILE *fp) {
  for (int j = 0; j < numSNVs; j++) {
    int k = 0;
    for (int i = 0; i < numCells + 1; i++) {
      if (i != healthyTip->label) {
        fprintf(fp, "%d", cellExtd[i].site[SNVsites[j]].numReads);
        k++;
      }

      if (k < numCells)
        fprintf(fp, "\t");
    }

    if (j < numSNVs - 1)
      fprintf(fp, "\n");
  }
}

/************************ PrintTrueSNV ***********************/
void PrintTrueSNV(FILE *fpAl, FILE *fpNr, FILE *fpMt) {
  for (int i = 0; i < numSNVs; i++) {
    for (int j = 0; j < numCells + 1; j++) {
      if (healthyTip->label != j) {
        const int *alleles = GetValIntToIntArr(
            genotypeHashMapIntToIntArr, cellExtd[j].site[SNVsites[i]].trueCode);

        PrintGenotypeInLettersToFile(
            fpAl, cellExtd[j].site[SNVsites[i]].trueNumExistingAlleles, alleles,
            NO, MISSING_CHAR_ALT, j > 0);

        PrintGenotypeInNumbersToFile(
            fpNr, fpMt, allSites[SNVsites[i]].referenceAllele,
            cellExtd[j].site[SNVsites[i]].trueNumExistingAlleles, alleles, NO,
            MISSING_CHAR_ALT, j > 0, 1);
      }
    }

    if (i < numSNVs - 1) {
      fprintf(fpAl, "\n");
      fprintf(fpNr, "\n");
      fprintf(fpMt, "\n");
    }
  }
}

void PrintEvolutionaryEventsSummary(FILE *fp) {
  /*
   * min [0] and max [1] number of:
   * mutations, back mutations, deletions, insertions
   */
  unsigned int info[4][2];

  info[0][0] = info[1][0] = info[2][0] = info[3][0] = UINT_MAX;
  info[0][1] = info[1][1] = info[2][1] = info[3][1] = 0;

  // Gather information
  for (int i = 0; i < numSNVs; i++) {
    for (int j = 0; j < numCells + 1; j++) {
      if (healthyTip->label != j) {
        for (int k = 0; k < 4; k++) {
          if (cellExtd[j].site[SNVsites[i]].pastEvolutionaryEvents[k] <
              info[k][0])
            info[k][0] =
                cellExtd[j].site[SNVsites[i]].pastEvolutionaryEvents[k];
          if (cellExtd[j].site[SNVsites[i]].pastEvolutionaryEvents[k] >
              info[k][1])
            info[k][1] =
                cellExtd[j].site[SNVsites[i]].pastEvolutionaryEvents[k];
        }
      }
    }
  }

  // Print information
  fprintf(fp, "#Cell-wise number of mutations = [%d, %d]\n", info[0][0],
          info[0][1]);
  fprintf(fp, "#Cell-wise number of back mutations = [%d, %d]\n", info[1][0],
          info[1][1]);
  fprintf(fp, "#Cell-wise number of deletions = [%d, %d]\n", info[2][0],
          info[2][1]);
  fprintf(fp, "#Cell-wise number of insertions = [%d, %d]\n", info[3][0],
          info[3][1]);
  fprintf(fp, "#Format of tab separated entries: mutations,back "
              "mutations,deletions,insertions\n");
  for (int i = 0; i < numSNVs; i++) {
    int cellsPrinted = 0;
    for (int j = 0; j < numCells + 1; j++) {
      if (healthyTip->label != j) {
        for (int k = 0; k < 4; k++) {
          fprintf(fp, "%d",
                  cellExtd[j].site[SNVsites[i]].pastEvolutionaryEvents[k]);
          if (k < 3)
            fprintf(fp, ",");
        }

        if (cellsPrinted++ < numCells - 1)
          fprintf(fp, "\t");
      }
    }

    if (i < numSNVs - 1)
      fprintf(fp, "\n");
  }
}

/************************ PrintSNVGenotypesExtd ***********************/
/* Prints observed genotypes at variable sites (SNVs) to a file */
void PrintSNVGenotypesExtd(FILE *fp, int toStderrFlag) {
  int i, j;

  if (doPrintAncestors == YES) {
    fprintf(fp, "%d %d\n", numCells + 1, numSNVs);
    if (toStderrFlag) {
      fprintf(stderr, "\n\n%d %d\n", numCells + 1, numSNVs);
    }
  }

  /* site information */
  fprintf(fp, "%d", SNVsites[0] + 1);
  if (toStderrFlag) {
    fprintf(stderr, "%d", SNVsites[0] + 1);
  }
  for (i = 1; i < numSNVs; i++) {
    fprintf(fp, " %d", SNVsites[i] + 1);
    if (toStderrFlag) {
      fprintf(stderr, " %d", SNVsites[i] + 1);
    }
  }
  fprintf(fp, "\n");
  if (toStderrFlag) {
    fprintf(stderr, "\n");
  }

  for (i = 0; i < numCells + 1; i++) {
    if (i == numCells) {
      fprintf(fp, "%s ", outCellName);
      if (toStderrFlag) {
        fprintf(stderr, "%s ", outCellName);
      }
    } else {
      if (doUserTree == NO) {
        fprintf(fp, "%s%04d ", inCellName, i + 1);
        if (toStderrFlag) {
          fprintf(stderr, "%s%04d ", inCellName, i + 1);
        }
      } else {
        fprintf(fp, "%-*s", stringPrecision, cellNames[i]);
        if (toStderrFlag) {
          fprintf(stderr, "%-*s", stringPrecision, cellNames[i]);
        }
      }
    }

    for (j = 0; j < numSNVs; j++) {
      PrintGenotypeInLettersToFile(
          fp, cellExtd[i].site[SNVsites[j]].observedNumExistingAlleles,
          cellExtd[i].site[SNVsites[j]].observedAlleles, NO, MISSING_CHAR_ALT,
          1);
    }
    fprintf(fp, "\n");
    if (toStderrFlag) {
      fprintf(stderr, "\n");
    }
  }

  if (doPrintAncestors == YES) {
    fprintf(fp, "%s ", outRootCellName);
    if (toStderrFlag) {
      fprintf(stderr, "%s ", outRootCellName);
    }
    for (j = 0; j < numSNVs; j++) {
      int thIndex = GetThreadIndex(SNVsites[j]);
      int relaSite = SNVsites[j] - threadStr[thIndex].start;
      PrintGenotypeInLettersToFile(
          fp,
          threadStr[thIndex].data[OUTGROUP_ROOT][relaSite].numExistingAlleles,
          threadStr[thIndex].data[OUTGROUP_ROOT][relaSite].alleles, NO,
          MISSING_CHAR_ALT, 1);
    }
    fprintf(fp, "\n%s ", inRootCellName);
    if (toStderrFlag) {
      fprintf(stderr, "\n%s ", inRootCellName);
    }

    for (j = 0; j < numSNVs; j++) {
      int thIndex = GetThreadIndex(SNVsites[j]);
      int relaSite = SNVsites[j] - threadStr[thIndex].start;
      PrintGenotypeInLettersToFile(
          fp,
          threadStr[thIndex].data[INGROUP_ROOT][relaSite].numExistingAlleles,
          threadStr[thIndex].data[INGROUP_ROOT][relaSite].alleles, NO,
          MISSING_CHAR_ALT, 1);
    }
    fprintf(fp, "\n");
    if (toStderrFlag) {
      fprintf(stderr, "\n");
    }
  }
}
/* end added code */

/***************************** PrintSNVHaplotypes
 * *******************************/
/* Prints obsserved/ML haplotypes at variable sites (SNVs), to a file */

static void PrintSNVHaplotypes(FILE *fp, int PrintTrueVariants) {
  int i, j;
  int *sites;

  if (PrintTrueVariants == YES)
    sites = variantSites;
  else
    sites = SNVsites;

  if (doPrintIUPAChaplotypes == YES) {
    if (doPrintAncestors == YES)
      fprintf(fp, "%d %d\n", numCells + 3, numSNVs);
    else
      fprintf(fp, "%d %d\n", numCells + 1, numSNVs);
  } else {
    if (doPrintAncestors == YES && doNGS == NO)
      fprintf(fp, "%d %d\n", 2 * (numCells + 3), numSNVs);
    else
      fprintf(fp, "%d %d\n", 2 * (numCells + 1), numSNVs);
  }

  /* site information */
  fprintf(fp, "%d", sites[0] + 1);
  for (i = 1; i < numSNVs; i++)
    fprintf(fp, " %d", sites[i] + 1);
  fprintf(fp, "\n");

  if (alphabet == DNA) {
    if (doPrintIUPAChaplotypes == YES) {
      for (i = 0; i < numCells + 1; i++) {
        /* print haplotype */
        if (i == numCells)
          fprintf(fp, "%s ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04d ", inCellName, i + 1);
          else
            fprintf(fp, "%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSNVs; j++) {
          if (doNGS == NO)
            fprintf(fp, "%c",
                    WhichIUPAC(data[MATERNAL][i][sites[j]],
                               data[PATERNAL][i][sites[j]]));
          else if (cell[i].hasDoublet == NO)
            fprintf(fp, "%c",
                    WhichIUPAC(cell[i].site[SNVsites[j]].MLmatAllele,
                               cell[i].site[SNVsites[j]].MLpatAllele));
          else
            fprintf(fp, "%c",
                    WhichIUPAC(cell[i].site[SNVsites[j]].MLmatAlleleDoublet,
                               cell[i].site[SNVsites[j]].MLpatAlleleDoublet));
        }
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES && doNGS == NO) {
        fprintf(fp, "%s ", outRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c",
                  WhichIUPAC(data[MATERNAL][OUTGROUP_ROOT][sites[j]],
                             data[PATERNAL][OUTGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n%s ", inRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c",
                  WhichIUPAC(data[MATERNAL][INGROUP_ROOT][sites[j]],
                             data[PATERNAL][INGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n");
      }
    } else // print maternal and paternal DNA haplotypes
    {
      for (i = 0; i < numCells + 1; i++) {
        /* print maternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sm  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dm  ", inCellName, i + 1);
          else
            fprintf(fp, "m%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSNVs; j++) {
          if (doNGS == NO)
            fprintf(fp, "%c", WhichNuc(data[MATERNAL][i][SNVsites[j]]));
          else if (cell[i].hasDoublet == NO)
            fprintf(fp, "%c", WhichNuc(cell[i].site[SNVsites[j]].MLmatAllele));
          else
            fprintf(fp, "%c",
                    WhichNuc(cell[i].site[SNVsites[j]].MLmatAlleleDoublet));
        }
        fprintf(fp, "\n");

        /* print paternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sp  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dp  ", inCellName, i + 1);
          else
            fprintf(fp, "p%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSNVs; j++) {
          if (doNGS == NO)
            fprintf(fp, "%c", WhichNuc(data[PATERNAL][i][SNVsites[j]]));
          else if (cell[i].hasDoublet == NO)
            fprintf(fp, "%c", WhichNuc(cell[i].site[SNVsites[j]].MLpatAllele));
          else
            fprintf(fp, "%c",
                    WhichNuc(cell[i].site[SNVsites[j]].MLpatAlleleDoublet));
        }
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES && doNGS == NO) {
        fprintf(fp, "%sm  ", outRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichNuc(data[MATERNAL][OUTGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n%sp  ", outRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichNuc(data[PATERNAL][OUTGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n%sm  ", inRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichNuc(data[MATERNAL][INGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n%sp  ", inRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichNuc(data[PATERNAL][INGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n");
      }
    }
  } else // print binary haplotypes
  {
    if (doPrintIUPAChaplotypes == YES) // print binary consensus SNV haplotypes
    {
      for (i = 0; i < numCells + 1; i++) {
        if (i == numCells)
          fprintf(fp, "%s ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04d ", inCellName, i + 1);
          else
            fprintf(fp, "%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c",
                  WhichConsensusBinary(data[MATERNAL][i][sites[j]],
                                       data[PATERNAL][i][sites[j]]));
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES) {
        fprintf(fp, "%s ", outRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(
              fp, "%c",
              WhichConsensusBinary(data[MATERNAL][OUTGROUP_ROOT][sites[j]],
                                   data[PATERNAL][OUTGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n%s ", inRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c",
                  WhichConsensusBinary(data[MATERNAL][INGROUP_ROOT][sites[j]],
                                       data[PATERNAL][INGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n");
      }
    } else // print maternal and paternal binary SNV haplotypes
    {
      for (i = 0; i < numCells + 1; i++) {
        /* print maternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sm  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dm  ", inCellName, i + 1);
          else
            fprintf(fp, "m%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichMut(data[MATERNAL][i][sites[j]]));
        fprintf(fp, "\n");

        /* print paternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sp  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dp  ", inCellName, i + 1);
          else
            fprintf(fp, "p%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichMut(data[PATERNAL][i][sites[j]]));
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES) {
        fprintf(fp, "%sm  ", outRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichMut(data[MATERNAL][OUTGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n%sp  ", outRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichMut(data[PATERNAL][OUTGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n%sm  ", inRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichMut(data[MATERNAL][INGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n%sp  ", inRootCellName);
        for (j = 0; j < numSNVs; j++)
          fprintf(fp, "%c", WhichMut(data[PATERNAL][INGROUP_ROOT][sites[j]]));
        fprintf(fp, "\n");
      }
    }
  }
}

/********************* PrintFullGenotypes ********************/
/* Prints observed/ML genotypes at all sites to a file */

static void PrintFullGenotypes(FILE *fp) {
  int i, j;

  if (doPrintAncestors == YES)
    fprintf(fp, "%d %d\n", numCells + 3, numSites);
  else
    fprintf(fp, "%d %d\n", numCells + 1, numSites);

  if (alphabet == DNA) {
    for (i = 0; i < numCells + 1; i++) {
      if (i == numCells)
        fprintf(fp, "%s ", outCellName);
      else {
        if (doUserTree == NO)
          fprintf(fp, "%s%04d ", inCellName, i + 1);
        else
          fprintf(fp, "%-*s", stringPrecision, cellNames[i]);
      }
      for (j = 0; j < numSites; j++) {
        if (doNGS == NO)
          fprintf(fp, " %c%c", WhichNuc(data[MATERNAL][i][j]),
                  WhichNuc(data[PATERNAL][i][j]));
        else if (cell[i].hasDoublet == NO)
          fprintf(fp, " %c%c", WhichNuc(cell[i].site[j].MLmatAllele),
                  WhichNuc(cell[i].site[j].MLpatAllele));
        else
          fprintf(fp, " %c%c", WhichNuc(cell[i].site[j].MLmatAlleleDoublet),
                  WhichNuc(cell[i].site[j].MLpatAlleleDoublet));
      }
      fprintf(fp, "\n");
    }

    /* print ancestral sequences, but not for the NGS case */
    if (doPrintAncestors == YES && doNGS == NO) {
      fprintf(fp, "%s ", outRootCellName);
      for (j = 0; j < numSites; j++)
        fprintf(fp, " %c%c", WhichNuc(data[MATERNAL][OUTGROUP_ROOT][j]),
                WhichNuc(data[PATERNAL][OUTGROUP_ROOT][j]));
      fprintf(fp, "\n%s ", inRootCellName);

      for (j = 0; j < numSites; j++)
        fprintf(fp, " %c%c", WhichNuc(data[MATERNAL][INGROUP_ROOT][j]),
                WhichNuc(data[PATERNAL][INGROUP_ROOT][j]));
      fprintf(fp, "\n");
    }
  } else // binary alphabet
  {
    for (i = 0; i < numCells + 1; i++) {
      if (i == numCells)
        fprintf(fp, "%s ", outCellName);
      else {
        if (doUserTree == NO)
          fprintf(fp, "%s%04d ", inCellName, i + 1);
        else
          fprintf(fp, "%-*s", stringPrecision, cellNames[i]);
      }
      for (j = 0; j < numSites; j++)
        fprintf(fp, " %c%c", WhichMut(data[0][i][j]), WhichMut(data[1][i][j]));
      fprintf(fp, "\n");
    }

    if (doPrintAncestors == YES) {
      fprintf(fp, "%s ", outRootCellName);
      for (j = 0; j < numSites; j++)
        fprintf(fp, " %d%d", data[MATERNAL][OUTGROUP_ROOT][j],
                data[PATERNAL][OUTGROUP_ROOT][j]);
      fprintf(fp, "\n%s ", inRootCellName);
      for (j = 0; j < numSites; j++)
        fprintf(fp, " %d%d", data[MATERNAL][INGROUP_ROOT][j],
                data[PATERNAL][INGROUP_ROOT][j]);
      fprintf(fp, "\n");
    }
  }
}

/* added code */
/********************* PrintFullGenotypesExtd ********************/
/* Prints observed/ML genotypes at all sites to a file */
void PrintFullGenotypesExtd(FILE *fp, int toStderrFlag) {
  int i, j;

  if (doPrintAncestors == YES) {
    fprintf(fp, "%d %d\n", numCells + 3, numSites);
    if (toStderrFlag) {
      fprintf(stderr, "\n\n%d %d\n", numCells + 3, numSites);
    }
  } else {
    fprintf(fp, "%d %d\n", numCells + 1, numSites);
    if (toStderrFlag) {
      fprintf(stderr, "\n\n%d %d\n", numCells + 1, numSites);
    }
  }

  for (i = 0; i < numCells + 1; i++) {
    if (i == numCells) {
      fprintf(fp, "%s ", outCellName);
      if (toStderrFlag) {
        fprintf(stderr, "%s ", outCellName);
      }
    } else {
      if (doUserTree == NO) {
        fprintf(fp, "%s%04d ", inCellName, i + 1);
        if (toStderrFlag) {
          fprintf(stderr, "%s%04d ", inCellName, i + 1);
        }
      } else {
        fprintf(fp, "%-*s", stringPrecision, cellNames[i]);
        if (toStderrFlag) {
          fprintf(stderr, "%-*s", stringPrecision, cellNames[i]);
        }
      }
    }

    for (j = 0; j < numSites; j++) {
      PrintGenotypeInLettersToFile(
          fp, cellExtd[i].site[j].observedNumExistingAlleles,
          cellExtd[i].site[j].observedAlleles, NO, MISSING_CHAR_ALT, 1);
    }
    fprintf(fp, "\n");
    if (toStderrFlag) {
      fprintf(stderr, "\n");
    }
  }

  /* print ancestral sequences */
  if (doPrintAncestors == YES) {
    fprintf(fp, "%s ", outRootCellName);
    if (toStderrFlag) {
      fprintf(stderr, "%s ", outRootCellName);
    }
    for (j = 0; j < numSites; j++) {
      int thIndex = GetThreadIndex(j);
      int relaSite = j - threadStr[thIndex].start;
      PrintGenotypeInLettersToFile(
          fp,
          threadStr[thIndex].data[OUTGROUP_ROOT][relaSite].numExistingAlleles,
          threadStr[thIndex].data[OUTGROUP_ROOT][relaSite].alleles, NO,
          MISSING_CHAR_ALT, 1);
    }
    fprintf(fp, "\n%s ", inRootCellName);
    if (toStderrFlag) {
      fprintf(stderr, "\n%s ", inRootCellName);
    }

    for (j = 0; j < numSites; j++) {
      int thIndex = GetThreadIndex(j);
      int relaSite = j - threadStr[thIndex].start;
      PrintGenotypeInLettersToFile(
          fp,
          threadStr[thIndex].data[INGROUP_ROOT][relaSite].numExistingAlleles,
          threadStr[thIndex].data[INGROUP_ROOT][relaSite].alleles, NO,
          MISSING_CHAR_ALT, 1);
    }
    fprintf(fp, "\n");
    if (toStderrFlag) {
      fprintf(stderr, "\n");
    }
  }
}
/* end added code */

/***************************** PrintFullHaplotypes
 * *******************************/
/* Prints observed/ML haplotypes for all sites (variable + invariable) to a
 * file
 */

static void PrintFullHaplotypes(FILE *fp) {
  int i, j;

  if (doPrintIUPAChaplotypes == NO) {
    if (doPrintAncestors == YES && doNGS == NO)
      fprintf(fp, "%d %d\n", 2 * (numCells + 3), numSites);
    else
      fprintf(fp, "%d %d\n", 2 * (numCells + 1), numSites);
  } else {
    if (doPrintAncestors == YES && doNGS == NO)
      fprintf(fp, "%d %d\n", numCells + 3, numSites);
    else
      fprintf(fp, "%d %d\n", numCells + 1, numSites);
  }

  if (alphabet == DNA) {
    if (doPrintIUPAChaplotypes == YES) {
      for (i = 0; i < numCells + 1; i++) {
        /* print IUPAC haplotype */
        if (i == numCells)
          fprintf(fp, "%s ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04d ", inCellName, i + 1);

          else
            fprintf(fp, "%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++) {
          if (doNGS == NO)
            fprintf(fp, "%c",
                    WhichIUPAC(data[MATERNAL][i][j], data[PATERNAL][i][j]));
          else if (cell[i].hasDoublet == NO)
            fprintf(fp, "%c",
                    WhichIUPAC(cell[i].site[j].MLmatAllele,
                               cell[i].site[j].MLpatAllele));
          else
            fprintf(fp, "%c",
                    WhichIUPAC(cell[i].site[j].MLmatAlleleDoublet,
                               cell[i].site[j].MLpatAlleleDoublet));
        }
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES && doNGS == NO) {
        fprintf(fp, "%s ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c",
                  WhichIUPAC(data[MATERNAL][OUTGROUP_ROOT][j],
                             data[PATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%s ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c",
                  WhichIUPAC(data[MATERNAL][INGROUP_ROOT][j],
                             data[PATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n");
      }
    } else // print maternal and paternal DNA haplotypes
    {
      for (i = 0; i < numCells + 1; i++) {
        /* print maternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sm  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dm  ", inCellName, i + 1);

          else
            fprintf(fp, "m%-*s ", stringPrecision, cellNames[i]);
        }

        for (j = 0; j < numSites; j++) {
          if (doNGS == NO)
            fprintf(fp, "%c", WhichNuc(data[MATERNAL][i][j]));
          else if (cell[i].hasDoublet == NO)
            fprintf(fp, "%c", WhichNuc(cell[i].site[j].MLmatAllele));
          else
            fprintf(fp, "%c", WhichNuc(cell[i].site[j].MLmatAlleleDoublet));
        }
        fprintf(fp, "\n");

        /* print paternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sp  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dp  ", inCellName, i + 1);
          else
            fprintf(fp, "p%-*s ", stringPrecision, cellNames[i]);
        }

        for (j = 0; j < numSites; j++) {
          if (doNGS == NO)
            fprintf(fp, "%c", WhichNuc(data[PATERNAL][i][j]));
          else if (cell[i].hasDoublet == NO)
            fprintf(fp, "%c", WhichNuc(cell[i].site[j].MLpatAllele));
          else
            fprintf(fp, "%c", WhichNuc(cell[i].site[j].MLpatAlleleDoublet));
        }
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES && doNGS == NO) {
        fprintf(fp, "%sm  ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichNuc(data[MATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%sp  ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichNuc(data[PATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%sm  ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichNuc(data[MATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n%sp  ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichNuc(data[PATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n");
      }
    }
  } else // print binary haplotypes
  {
    if (doPrintIUPAChaplotypes == YES) // print binary consensus haplotypes
    {
      for (i = 0; i < numCells + 1; i++) {
        if (i == numCells)
          fprintf(fp, "%s ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04d ", inCellName, i + 1);
          else
            fprintf(fp, "%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++)
          fprintf(
              fp, "%c",
              WhichConsensusBinary(data[MATERNAL][i][j], data[PATERNAL][i][j]));
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES) {
        fprintf(fp, "%s ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c",
                  WhichConsensusBinary(data[MATERNAL][OUTGROUP_ROOT][j],
                                       data[PATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%s ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c",
                  WhichConsensusBinary(data[MATERNAL][INGROUP_ROOT][j],
                                       data[PATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n");
      }
    } else // print maternal and paternal binary haplotypes
    {
      for (i = 0; i < numCells + 1; i++) {
        /* print maternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sm  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dm  ", inCellName, i + 1);
          else
            fprintf(fp, "m%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[MATERNAL][i][j]));
        fprintf(fp, "\n");

        /* print paternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sp  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dp  ", inCellName, i + 1);
          else
            fprintf(fp, "p%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[PATERNAL][i][j]));
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES) {
        fprintf(fp, "%sm  ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[MATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%sp  ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[PATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%sm  ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[MATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n%sp  ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[PATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n");
      }
    }
  }
}

/***************************** PrintTrueFullHaplotypes
 * *******************************/
/* Prints observed/ML haplotypes for all sites (variable + invariable) to a
 * file
 */

static void PrintTrueFullHaplotypes(FILE *fp) {
  int i, j;

  if (doPrintIUPAChaplotypes == NO) {
    if (doPrintAncestors == YES)
      fprintf(fp, "%d %d\n", 2 * (numCells + 3), numSites);
    else
      fprintf(fp, "%d %d\n", 2 * (numCells + 1), numSites);
  } else {
    if (doPrintAncestors == YES)
      fprintf(fp, "%d %d\n", numCells + 3, numSites);
    else
      fprintf(fp, "%d %d\n", numCells + 1, numSites);
  }

  if (alphabet == DNA) {
    if (doPrintIUPAChaplotypes == YES) {
      for (i = 0; i < numCells + 1; i++) {
        /* print IUPAC haplotype */
        if (i == numCells)
          fprintf(fp, "%s ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04d ", inCellName, i + 1);
          else
            fprintf(fp, "%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++) {
          fprintf(fp, "%c",
                  WhichIUPAC(data[MATERNAL][i][j], data[PATERNAL][i][j]));
        }
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES) {
        fprintf(fp, "%s ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c",
                  WhichIUPAC(data[MATERNAL][OUTGROUP_ROOT][j],
                             data[PATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%s ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c",
                  WhichIUPAC(data[MATERNAL][INGROUP_ROOT][j],
                             data[PATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n");
      }
    } else // print maternal and paternal DNA haplotypes
    {
      for (i = 0; i < numCells + 1; i++) {
        /* print maternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sm  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dm  ", inCellName, i + 1);
          else
            fprintf(fp, "m%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++) {
          fprintf(fp, "%c", WhichNuc(data[MATERNAL][i][j]));
        }
        fprintf(fp, "\n");

        /* print paternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sp  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dp  ", inCellName, i + 1);
          else
            fprintf(fp, "p%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++) {
          fprintf(fp, "%c", WhichNuc(data[PATERNAL][i][j]));
        }
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES) {
        fprintf(fp, "%sm  ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichNuc(data[MATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%sp  ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichNuc(data[PATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%sm  ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichNuc(data[MATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n%sp  ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichNuc(data[PATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n");
      }
    }
  } else // print binary haplotypes
  {
    if (doPrintIUPAChaplotypes == YES) // print binary consensus haplotypes
    {
      for (i = 0; i < numCells + 1; i++) {
        if (i == numCells)
          fprintf(fp, "%s ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04d ", inCellName, i + 1);
          else
            fprintf(fp, "%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++)
          fprintf(
              fp, "%c",
              WhichConsensusBinary(data[MATERNAL][i][j], data[PATERNAL][i][j]));
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES) {
        fprintf(fp, "%s ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c",
                  WhichConsensusBinary(data[MATERNAL][OUTGROUP_ROOT][j],
                                       data[PATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%s ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c",
                  WhichConsensusBinary(data[MATERNAL][INGROUP_ROOT][j],
                                       data[PATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n");
      }
    } else // print maternal and paternal binary haplotypes
    {
      for (i = 0; i < numCells + 1; i++) {
        /* print maternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sm  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dm  ", inCellName, i + 1);
          else
            fprintf(fp, "m%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[MATERNAL][i][j]));
        fprintf(fp, "\n");

        /* print paternal haplotype */
        if (i == numCells)
          fprintf(fp, "%sp  ", outCellName);
        else {
          if (doUserTree == NO)
            fprintf(fp, "%s%04dp  ", inCellName, i + 1);
          else
            fprintf(fp, "p%-*s ", stringPrecision, cellNames[i]);
        }
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[PATERNAL][i][j]));
        fprintf(fp, "\n");
      }

      if (doPrintAncestors == YES) {
        fprintf(fp, "%sm  ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[MATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%sp  ", outRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[PATERNAL][OUTGROUP_ROOT][j]));
        fprintf(fp, "\n%sm  ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[MATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n%sp  ", inRootCellName);
        for (j = 0; j < numSites; j++)
          fprintf(fp, "%c", WhichMut(data[PATERNAL][INGROUP_ROOT][j]));
        fprintf(fp, "\n");
      }
    }
  }
}

/********************* WhichIUPAC ************************/
/* Returns the IUPAC representation of the genotype */
/*
 UPAC nucleotide code	Base
 A	Adenine
 C	Cytosine
 G	Guanine
 T (or U)	Thymine (or Uracil)
 R	A or G
 Y	C or T
 S	G or C
 W	A or T
 K	G or T
 M	A or C
 B	C or G or T
 D	A or G or T
 H	A or C or T
 V	A or C or G
 N	unknown state
 . or -	gap

This is what we do:

 A/A => A
 A/C => M
 A/G => R
 A/T => W
 A/_ => a
 A/- => a

 C/A => M
 C/C => C
 C/G => S
 C/T => Y
 C/_ => c
 C/- => c

 G/A => R
 G/C => S
 G/G => G
 G/T => K
 G/_ => g
 G/- => g

 T/A => W
 T/C => Y
 T/G => K
 T/T => T
 T/_ => t
 T/- => t

 _/A => a
 _/C => c
 _/G => g
 _/T => t
 _/_ => -
 _/- => -


 -/A => a
 -/C => c
 -/G => g
 -/T => t
 -/_ => -
 -/- => -

 */

char WhichIUPAC(int allele1, int allele2) {
  if (allele1 == 0) {
    if (allele2 == 0) // AA
      return ('A');
    else if (allele2 == 1) // AC
      return ('M');
    else if (allele2 == 2) // AG
      return ('R');
    else if (allele2 == 3) // AT
      return ('W');
    else if (allele2 == ADO) // A_
      return ('a');
    else if (allele2 == DELETION) // A–
      return ('a');
    else
      return ('N');
  } else if (allele1 == 1) {
    if (allele2 == 0) // CA
      return ('M');
    else if (allele2 == 1) // CC
      return ('C');
    else if (allele2 == 2) // CG
      return ('S');
    else if (allele2 == 3) // CT
      return ('Y');
    else if (allele2 == ADO) // C?
      return ('c');
    else if (allele2 == DELETION) // C–
      return ('c');
    else
      return ('N');
  } else if (allele1 == 2) {
    if (allele2 == 0) // GA
      return ('R');
    else if (allele2 == 1) // GC
      return ('S');
    else if (allele2 == 2) // GG
      return ('G');
    else if (allele2 == 3) // GT
      return ('K');
    else if (allele2 == ADO) // G?
      return ('g');
    else if (allele2 == DELETION) // G–
      return ('g');
    else
      return ('N');
  } else if (allele1 == 3) {
    if (allele2 == 0) // TA
      return ('W');
    else if (allele2 == 1) // TC
      return ('Y');
    else if (allele2 == 2) // TG
      return ('K');
    else if (allele2 == 3) // TT
      return ('T');
    else if (allele2 == ADO) // T?
      return ('t');
    else if (allele2 == DELETION) // T–
      return ('t');
    else
      return ('N');
  } else if (allele1 == ADO) {
    if (allele2 == 0) //_A
      return ('a');
    else if (allele2 == 1) //_C
      return ('c');
    else if (allele2 == 2) //_G
      return ('g');
    else if (allele2 == 3) //_T
      return ('t');
    else if (allele2 == ADO) //__
      return ('-');
    else if (allele2 == DELETION) //_-
      return ('-');
    else
      return ('N');
  } else if (allele1 == DELETION) {
    if (allele2 == 0) //-A
      return ('a');
    else if (allele2 == 1) //-C
      return ('c');
    else if (allele2 == 2) //-G
      return ('g');
    else if (allele2 == 3) //-T
      return ('t');
    else if (allele2 == ADO) //-_
      return ('-');
    else if (allele2 == DELETION) //--
      return ('-');
    else
      return ('N');
  } else
    return ('N');
}

/********************* WhichConsensusBinary ************************/
/* Returns a consensus representation of the binary genotype */
/*
 0/0 => 0
 0/1 => 1
 1/0 => 1
 1/1 => 2

 0/_ => 0
 _/0 => 0

 1/_ => 2
 _/1 => 2

 _/_ => -
 */

char WhichConsensusBinary(int allele1, int allele2) {
  if (allele1 == 0) {
    if (allele2 == 0) // 00
      return ('0');
    else if (allele2 == 1) // 01
      return ('1');
    else if (allele2 == ADO) // 0?
      return ('0');
    else if (allele2 == DELETION) // 0-
      return ('0');
    else
      return ('N');
  } else if (allele1 == 1) {
    if (allele2 == 0) // 10
      return ('1');
    else if (allele2 == 1) // 11
      return ('2');
    else if (allele2 == ADO) // 1?
      return ('2');
    else if (allele2 == DELETION) // 1-
      return ('2');
    else
      return ('N');
  } else if (allele1 == ADO) {
    if (allele2 == 0) //?0
      return ('0');
    else if (allele2 == 1) //?1
      return ('2');
    else if (allele2 == ADO) //??
      return ('-');
    else if (allele2 == DELETION) //?-
      return ('-');
    else
      return ('N');
  } else if (allele1 == DELETION) {
    if (allele2 == 0) //-0
      return ('0');
    else if (allele2 == 1) //-1
      return ('2');
    else if (allele2 == ADO) //-?
      return ('-');
    else if (allele2 == DELETION) //--
      return ('-');
    else
      return ('N');
  } else
    return ('N');
}

/********************* WhichNucWhichNucChar ************************/
/* Returns integer representation for character nucleotides */

int WhichNucChar(char nucleotide) {
  if (nucleotide == 'A')
    return (A);
  else if (nucleotide == 'C')
    return (C);
  else if (nucleotide == 'G')
    return (G);
  else if (nucleotide == 'T')
    return (T);
  else if (nucleotide == '?')
    return (ADO);
  else if (nucleotide == '-')
    return (DELETION);
  else {
    fprintf(stderr, "\nERROR in WhicNucChar: nucleotide = %c\n", nucleotide);
    exit(-1);
  }
}

/********************* WhichNuc ************************/
/* Returns character representation for nucleotides */

char WhichNuc(int nucleotide) {
  if (nucleotide == A)
    return ('A');
  else if (nucleotide == C)
    return ('C');
  else if (nucleotide == G)
    return ('G');
  else if (nucleotide == T)
    return ('T');
  else if (nucleotide == ADO)
    return ('?');
  else if (nucleotide == DELETION)
    return ('-');
  else
    return ('N');
}

/********************* WhichMut ************************/
/* Returns character representation for binary data */

char WhichMut(int state) {
  if (state == 0)
    return ('0');
  else if (state == 1)
    return ('1');
  else if (state == ADO)
    return ('?');
  else if (state == DELETION)
    return ('-');
  else
    return ('N');
}

/***************************** PrintHeader *******************************/
/* Prints program header */

static void PrintHeader(FILE *fp) {
  fprintf(fp, "________________________________________________________________"
              "______\n\n");
  fprintf(fp, "Cell coalescent simulation - %s", PROGRAM_NAME);
  fprintf(fp, "  %s", VERSION_NUMBER);
  fprintf(fp, "\n(c) 2019 David Posada - dposada@uvigo.es");
  fprintf(fp, "\n______________________________________________________________"
              "________\n\n");
}

/***************************** PrintDate *******************************/
/* Prints current date and time  */

static void PrintDate(FILE *fp) {
  time_t now;
  char *date;

  now = time(NULL);
  date = ctime(&now);
  fprintf(fp, "%s", date);
}

/***************************** PrintRunInformation
 * *******************************/
/* Prints a summary of run settings */

static void PrintRunInformation(FILE *fp) {
  int i;

  fprintf(fp, "\n\nRun settings\n------------");
  fprintf(fp, "\n[Parameter file = %s]\n", parameterFile);
  if (noisy > 1)
    fprintf(fp, "\n[Assumptions in brackets]\n");
  fprintf(fp, "\n Seed                                         =   %-3ld",
          originalSeed);
  fprintf(fp, "\n Number replicate data sets                   =   %-3d",
          numDataSets);
  fprintf(fp, "\n Number of cells                              =   %-3d",
          numCells);
  if (doUserGenome == YES)
    fprintf(fp, "\n Number of sites in the user genome           =   %-3d",
            numSites);
  else
    fprintf(fp, "\n Number of genomic sites                      =   %-3d",
            numSites);

  fprintf(fp, "\n\nDemographics");
  fprintf(fp, "\n Effective population size                    =   %-3d", N);
  if (doExponential == YES)
    fprintf(fp, "\n Exponential growth rate                      =   %2.1e",
            growthRate);
  else if (doDemographics == YES) {
    fprintf(fp, "\n Period   Nbegin     Nend    Duration    Growth rate");
    for (i = 1; i <= numPeriods; i++)
      fprintf(fp, "\n   %d    %7d   %7d     %7d       %+2.1e", i, Nbegin[i],
              Nend[i], cumDuration[i] - cumDuration[i - 1], periodGrowth[i]);
  }
  if (doUserTree == YES)
    fprintf(fp, "\n\nUser tree file                                =   %s",
            userTreeFile);
  else {
    fprintf(fp, "\n\nCoalescent");
    // fprintf (fp, "\n Mean number of coalescence events            = %3.2f",
    // meanNumCA);
    fprintf(fp, "\n Mean time to the MRCA (2N generations)       =   %3.2f",
            meanTMRCA);
    if (noisy > 1)
      fprintf(fp, "\n Exp time to the MRCA [constant Ne]           =   %3.2f",
              expTMRCA);
    if (numDataSets > 1) {
      fprintf(fp, "\n Variance time MRCA                           =   %3.2f",
              varTMRCA);
      if (noisy > 1)
        fprintf(fp, "\n Exp variance time MRCA [constant Ne]         =   %3.2f",
                expVarTMRCA);
    }

    if (rateVarAmongLineages == YES)
      fprintf(fp, "\n Rate variation among branches (alpha)        =   %-3.2f",
              alphaBranches);
    else
      fprintf(fp, "\n No rate variation among branches");
  }

  fprintf(fp, "\n\nUser-defined branches");
  fprintf(fp, "\n Transforming branch length ratio             =   %2.1e",
          transformingBranchLengthRatio);
  fprintf(fp, "\n Healthy tip branch length ratio              =   %2.1e",
          healthyTipBranchLengthRatio);

  if (doUserGenome == YES)
    fprintf(fp, "\n\nUser genome file                              =   %s",
            userGenomeFile);

  if (doSimulateData == YES) {
    fprintf(fp, "\n\nMutation models");
    fprintf(fp, "\n Character alphabet                           =   ");
    if (alphabet == BINARY)
      fprintf(fp, "0/1");
    else if (alphabet == DNA)
      fprintf(fp, "DNA");

    if (doGeneticSignatures == YES) {
      fprintf(fp, "\n Number of trinucleotide genetic signatures   =   %d",
              numUserSignatures);
      for (i = 0; i < numUserSignatures; i++)
        fprintf(fp,
                "\n  signature %02d weight                         =   %4.2f",
                signatureID[i], signatureWeight[i]);
    } else if (propAltModelSites == 0)
      fprintf(fp, "\n All sites evolving under the ISM diploid");

    else if (propAltModelSites > 0) {
      fprintf(fp, "\n Proportion of ISM diploid sites              =   %-3.2f",
              1.0 - propAltModelSites);
      fprintf(fp, "\n Proportion of non-ISM diploid sites          =   %-3.2f",
              propAltModelSites);
      fprintf(fp, "\n Non-ISM diploid mutation model               =   ");
      if (altModel == 0) /* ISM haploid model */
        fprintf(fp, "ISM haploid");
      else if (altModel == 1) /* Markov (Mk2) model */
        fprintf(fp, "Markov 0/1");
      else if (altModel == 2) /* finite DNA model */
      {
        fprintf(fp, "Markov finite DNA: ");
        if (doJC == YES)
          fprintf(fp, "JC");
        else if (doHKY == YES)
          fprintf(fp, "HKY");
        else if (doGTR == YES)
          fprintf(fp, "GTR");
        else if (doGTnR == YES)
          fprintf(fp, "GTnR");
      }

      if (rateVarAmongSites == YES)
        fprintf(fp,
                "\n Rate variation among sites (alpha)           =   %-3.2f",
                alphaSites);
      else
        fprintf(fp, "\n Equal rates among sites");
    }

    if (alphabet == DNA && propAltModelSites > 0 && altModel == 2 &&
        doGeneticSignatures == NO) {
      fprintf(fp,
              "\n Base frequencies (ACGT)                      =   %3.2f %3.2f "
              "%3.2f %3.2f",
              freq[0], freq[1], freq[2], freq[3]);
      if (doHKY == YES)
        fprintf(fp,
                "\n Transition/transversion ratio                =   %3.2f  "
                "(kappa = %3.2f)",
                titv, kappa);
      if (thereIsMij == YES) {
        fprintf(fp,
                "\n Mutation rate matrix                         =   %3.2f "
                "%3.2f %3.2f %3.2f",
                Mij[0][0], Mij[0][1], Mij[0][2], Mij[0][3]);
        fprintf(fp,
                "\n                                                  %3.2f "
                "%3.2f %3.2f %3.2f",
                Mij[1][0], Mij[1][1], Mij[1][2], Mij[1][3]);
        fprintf(fp,
                "\n                                                  %3.2f "
                "%3.2f %3.2f %3.2f",
                Mij[2][0], Mij[2][1], Mij[2][2], Mij[2][3]);
        fprintf(fp,
                "\n                                                  %3.2f "
                "%3.2f %3.2f %3.2f",
                Mij[3][0], Mij[3][1], Mij[3][2], Mij[3][3]);
      }
    }
    if (SNPrate > 0)
      fprintf(fp, "\n SNPrate                                      =   %2.1e",
              SNPrate);

    fprintf(fp, "\n\nSNV rates");
    if (doSimulateFixedNumMutations == YES)
      fprintf(fp, "\n Fixed number of mutations                    =   %-3d",
              numFixedMutations);
    else {
      fprintf(fp, "\n Mutation rate                                =   %2.1e",
              mutationRate);
      fprintf(fp, "\n Population mutation parameter (4Nu2L)        =   %3.2f",
              theta);
      if (noisy > 1)
        fprintf(fp, "\n Exp num mutation events [Ne=cte,ISMh]        =   %3.2f",
                expNumMU);
    }
    fprintf(fp, "\n Mean number of mutation events               =   %3.2f",
            meanNumMU);
    if (numDataSets > 1) {
      if (doSimulateFixedNumMutations == NO) {
        if (noisy > 1)
          fprintf(fp,
                  "\n Exp variance num mut events [Ne=cte,ISM]     =   %3.2f",
                  expVarNumMU);
        fprintf(fp, "\n Variance number of mutation events           =   %3.2f",
                varNumMU);
      }
    }
    fprintf(fp, "\n Mean number of SNVs                          =   %3.2f",
            meanNumSNVs);
    if (numDataSets > 1) {
      fprintf(fp, "\n Variance number of SNVs                      =   %3.2f",
              varNumSNVs);
      fprintf(fp, "\n Proportion of reps without SNV sites         =   %3.2f",
              zeroSNVs / (double)numDataSets);
    }

    if (deletionRate > 0) {
      fprintf(fp, "\n\nDeletion rates");
      fprintf(fp, "\n Deletion rate                                =   %2.1e",
              deletionRate);
      fprintf(fp, "\n Mean number of deletion events               =   %3.2f",
              meanNumDEL);
      if (numDataSets > 1)
        fprintf(fp, "\n Variance number of deletion events           =   %3.2f",
                varNumDEL);
    }

    if (CNLOHrate > 0) {
      fprintf(fp, "\n\nCopy-neutral LOH rates");
      fprintf(fp, "\n CN-LOH rate                                  =   %2.1e",
              CNLOHrate);
      fprintf(fp, "\n Mean number of CN-LOH events                 =   %3.2f",
              meanNumCNLOH);
      if (numDataSets > 1)
        fprintf(fp, "\n Variance number of CN-LOH events             =   %3.2f",
                varNumCNLOH);
    }

    if (genotypingError > 0 || doNGS == YES)
      fprintf(fp, "\n\nNGS");
    else
      fprintf(fp,
              "\n\nNGS genotype errors/read counts are not being simulated");

    if (genotypingError > 0) {
      fprintf(fp, "\n Genotype error                               =   %2.1e",
              genotypingError);
      // fprintf (fp, "\n Exp number of FP SNVs due to genotype errors =
      // %3.2f", (1 - pow(1-genotypingError, numCells+1)) * numSites);
    }
    if (ADOrate > 0) {
      fprintf(fp, "\n Alellic dropout                              =   %2.1e",
              ADOrate);
      if (ADOvarAmongSites == YES)
        fprintf(fp, "\n  ADO variation among sites (alpha)           =   %2.1e",
                alphaADOsites);
      else
        fprintf(fp, "\n  ADO is constant among sites");
      if (ADOvarAmongCells == YES)
        fprintf(fp, "\n  ADO variation among cells (alpha)           =   %2.1e",
                alphaADOcells);
      else
        fprintf(fp, "\n  ADO is constant among cells");
    }

    if (doNGS == YES) {
      fprintf(fp, "\n Sequencing coverage                          =   %fX",
              allelicCovMean);
      if (rateVarCoverage == YES)
        fprintf(fp,
                "\n Coverage dispersion (alpha)                  =   %-3.2f",
                alphaCoverage);
      else
        fprintf(fp, "\n  [coverage follows a Poisson distribution]");
      fprintf(fp, "\n Maternal allelic imbalance                   =   %-3.2f",
              allelicImbalance);
      fprintf(fp, "\n Haploid coverage                             =   %-3.2f",
              haploidCoverageReduction);
      if (meanAmplificationError == 0)
        fprintf(fp, "\n Amplification error                          =   %2.1e",
                meanAmplificationError);
      else {
        fprintf(fp, "\n Amplification error");
        if (simulateOnlyTwoTemplates == YES)
          fprintf(fp, " (2-template model)");
        else
          fprintf(fp, " (4-template model)");
        fprintf(fp, "\n  Mean                                        =   %2.1e",
                meanAmplificationError);
        fprintf(fp, "\n  Variance                                    =   %2.1e",
                varAmplificationError);
      }
      fprintf(fp, "\n Sequencing error                             =   %2.1e",
              sequencingError);
      fprintf (fp,"\n-B: doublet rate per cell; mean, var =  %2.1e, %2.1e", meanDoubletRate, varDoubletRate);

      if (thereIsEij == YES) {
        fprintf(fp,
                "\n NGS error rate matrix                        =   %3.2f "
                "%3.2f %3.2f %3.2f",
                Eij[0][0], Eij[0][1], Eij[0][2], Eij[0][3]);
        fprintf(fp,
                "\n                                                  %3.2f "
                "%3.2f %3.2f %3.2f",
                Eij[1][0], Eij[1][1], Eij[1][2], Eij[1][3]);
        fprintf(fp,
                "\n                                                  %3.2f "
                "%3.2f %3.2f %3.2f",
                Eij[2][0], Eij[2][1], Eij[2][2], Eij[2][3]);
        fprintf(fp,
                "\n                                                  %3.2f "
                "%3.2f %3.2f %3.2f",
                Eij[3][0], Eij[3][1], Eij[3][2], Eij[3][3]);
      }

      // fprintf (fp, "\n Probability of untrue genotype calls         =
      // %3.2f\n (see documentation)", cumCountMLgenotypeErrors/numDataSets);
    }
  }
}

/* added code */
/***************************** PrintRunInformation
 * *******************************/
/* Prints a summary of run settings */
void PrintRunInformationExtd(FILE *fp) {
  fprintf(fp, "\n\nRun settings\n------------");
  fprintf(fp, "\n[Parameter file = %s]\n", parameterFile);
  if (noisy > 1)
    fprintf(fp, "\n[Assumptions in brackets]\n");
  fprintf(fp, "\n Seed                                         =   %-3ld",
          originalSeed);
  fprintf(fp, "\n Number replicate data sets                   =   %-3d",
          numDataSets);
  fprintf(fp, "\n Number of cells                              =   %-3d",
          numCells);
  if (doUserGenome == YES)
    fprintf(fp, "\n Number of sites in the user genome           =   %-3d",
            numSites);
  else
    fprintf(fp, "\n Number of genomic sites                      =   %-3d",
            numSites);

  fprintf(fp, "\n\nDemographics");
  fprintf(fp, "\n Effective population size                    =   %-3d", N);
  if (doExponential == YES)
    fprintf(fp, "\n Exponential growth rate                      =   %2.1e",
            growthRate);

  if (doUserTree == YES)
    fprintf(fp, "\n\nUser tree file                                =   %s",
            userTreeFile);
  else {
    fprintf(fp, "\n\nCoalescent");
    // fprintf (fp, "\n Mean number of coalescence events            = %3.2f",
    // meanNumCA);
    fprintf(fp, "\n Mean time to the MRCA (2N generations)       =   %3.2f",
            meanTMRCA);
    if (noisy > 1)
      fprintf(fp, "\n Exp time to the MRCA [constant Ne]           =   %3.2f",
              expTMRCA);
    if (numDataSets > 1) {
      fprintf(fp, "\n Variance time MRCA                           =   %3.2f",
              varTMRCA);
      if (noisy > 1)
        fprintf(fp, "\n Exp variance time MRCA [constant Ne]         =   %3.2f",
                expVarTMRCA);
    }

    if (rateVarAmongLineages == YES)
      fprintf(fp, "\n Rate variation among branches (alpha)        =   %-3.2f",
              alphaBranches);
    else
      fprintf(fp, "\n No rate variation among branches");
  }

  fprintf(fp, "\n\nUser-defined branches");
  fprintf(fp, "\n Transforming branch length ratio             =   %2.1e",
          transformingBranchLengthRatio);
  fprintf(fp, "\n Healthy tip branch length ratio              =   %2.1e",
          healthyTipBranchLengthRatio);

  if (doUserGenome == YES)
    fprintf(fp, "\n\nUser genome file                              =   %s",
            userGenomeFile);

  fprintf(fp, "\n\nMutation models");
  fprintf(fp, "\n Character alphabet                           =   DNA");

  fprintf(fp, "\n All sites evolving under the finite-sites "
              "mutation-deletion-insertion model");

  if (rateVarAmongSites == YES)
    fprintf(fp, "\n Rate variation among sites (alpha)           =   %-3.2f",
            alphaSites);
  else
    fprintf(fp, "\n Equal rates among sites");

  fprintf(fp,
          "\n Base frequencies (ACGT)                      =   %3.2f %3.2f "
          "%3.2f %3.2f",
          freq[0], freq[1], freq[2], freq[3]);

  fprintf(fp, "\n\nSNV rates");
  fprintf(fp, "\n Mutation rate                                =   %2.1e",
          mutationRate);
  fprintf(fp, "\n Mean number of mutation events               =   %3.2f",
          meanNumMU);
  if (numDataSets > 1)
    fprintf(fp, "\n Variance number of mutation events           =   %3.2f",
            varNumMU);
  if (deletionRate > 0) {
    fprintf(fp, "\n Deletion rate                                =   %2.1e",
            deletionRate);
    fprintf(fp, "\n Mean number of deletion events               =   %3.2f",
            meanNumDEL);
    if (numDataSets > 1)
      fprintf(fp, "\n Variance number of deletion events           =   %3.2f",
              varNumDEL);
  }
  if (insertionRate > 0) {
    fprintf(fp, "\n Insertion rate                               =   %2.1e",
            insertionRate);
    fprintf(fp, "\n Mean number of insertion events              =   %3.2f",
            meanNumINS);
    if (numDataSets > 1)
      fprintf(fp, "\n Variance number of insertion events          =   %3.2f",
              varNumINS);
  }

  fprintf(fp, "\n Mean number of SNVs                          =   %3.2f",
          meanNumSNVs);
  if (numDataSets > 1) {
    fprintf(fp, "\n Variance number of SNVs                      =   %3.2f",
            varNumSNVs);
    fprintf(fp, "\n Proportion of reps without SNV sites         =   %3.2f",
            zeroSNVs / (double)numDataSets);
  }

  fprintf(fp, "\n\nNGS");

  if (ADOrate > 0) {
    fprintf(fp, "\n Alellic dropout                              =   %2.1e",
            ADOrate);
    if (ADOtype == 1) {
      fprintf(fp,
              "\n Alellic dropout type                         =   single ADO");
    } else if (ADOtype == 2) {
      fprintf(fp,
              "\n Alellic dropout type                         =   locus ADO");
    }
    if (ADOvarAmongSites == YES)
      fprintf(fp, "\n  ADO variation among sites (alpha)           =   %2.1e",
              alphaADOsites);
    else
      fprintf(fp, "\n  ADO is constant among sites");
    if (ADOvarAmongCells == YES)
      fprintf(fp, "\n  ADO variation among cells (alpha)           =   %2.1e",
              alphaADOcells);
    else
      fprintf(fp, "\n  ADO is constant among cells");
  }

  fprintf(fp, "\n Sequencing coverage                          =   %fX",
          allelicCovMean);
  fprintf(fp, "\n Coverage dispersion                          =   %-3.2f",
          allelicRawVarMean);

  fprintf(fp, "\n Amplification error");
  if (fixedAmplificationError) {
    fprintf(fp, "\n  Fixed amplification error                   =   %2.1e",
            meanAmplificationError);
  } else {
    fprintf(fp, " (4-template model)");
    fprintf(fp, "\n  Mean                                        =   %2.1e",
            meanAmplificationError);
    fprintf(fp, "\n  Variance                                    =   %2.1e",
            varAmplificationError);
  }

  fprintf(fp, "\n Sequencing error                             =   %2.1e",
          sequencingError);

  fprintf(fp, "\n Sequencing coverage size factor");
  fprintf(fp, "\n  Mean                                        =   %2.1e",
          meanSizeFactor);
  fprintf(fp, "\n  Variance                                    =   %2.1e",
          varSizeFactor);

  fprintf(fp, "\n Overdispersion for each situation");
  for (int i = 0; i < 2; ++i) {
    fprintf(fp, "\n  Situation %d                                 =  %-3.2f",
            i + 1, shapeController[i]);
  }

  if (thereIsEij == YES) {
    fprintf(fp,
            "\n NGS error rate matrix                        =   %3.2f "
            "%3.2f %3.2f %3.2f",
            Eij[0][0], Eij[0][1], Eij[0][2], Eij[0][3]);
    fprintf(fp,
            "\n                                                  %3.2f "
            "%3.2f %3.2f %3.2f",
            Eij[1][0], Eij[1][1], Eij[1][2], Eij[1][3]);
    fprintf(fp,
            "\n                                                  %3.2f "
            "%3.2f %3.2f %3.2f",
            Eij[2][0], Eij[2][1], Eij[2][2], Eij[2][3]);
    fprintf(fp,
            "\n                                                  %3.2f "
            "%3.2f %3.2f %3.2f",
            Eij[3][0], Eij[3][1], Eij[3][2], Eij[3][3]);
  }
}
/* end added code */

/***************** PrintCommandLine **********************/
/*	Prints the settings used in the simulations */
static void PrintCommandLine(FILE *fp, int argc, char **argv) {
  int i;

  if (readingParameterFile == NO)
    fprintf(fp, "Command line arguments (some might not affect) = %s",
            CommandLine);
  else {
    fprintf(fp, "Parameter file [\"%s\"] arguments (some might not affect) = ",
            parameterFile);

    /* Coalescent */
    fprintf(fp, " -n%d -s%d -l%d -e%d -g%2.1e -h%d", numDataSets, numCells,
            numSites, N, growthRate, numPeriods);
    for (i = 1; i <= numPeriods; i++)
      fprintf(fp, " %d %d %d", Nbegin[i], Nend[i],
              cumDuration[i] - cumDuration[i - 1]);

    /* Post-coalescent */
    fprintf(fp, " -k%2.1e", transformingBranchLengthRatio);
    fprintf(fp, " -q%2.1e", healthyTipBranchLengthRatio);
    fprintf(fp, " -i%6.4f", alphaBranches);

    /* Mutation model */
    fprintf(fp, " -b%d", alphabet);
    fprintf(fp, " -c%2.1e", SNPrate);
    fprintf(fp, " -u%2.1e", mutationRate);
    fprintf(fp, " -d%2.1e", deletionRate);
    fprintf(fp, " -H%2.1e", CNLOHrate);
    fprintf(fp, " -j%d", numFixedMutations);
    fprintf(fp, " -S%d", numUserSignatures);
    for (i = 0; i < numUserSignatures; i++)
      fprintf(fp, " %d %4.2f", signatureID[i], signatureWeight[i]);
    fprintf(fp, " -m%d", altModel);
    fprintf(fp, " -p%6.4f", propAltModelSites);
    fprintf(fp, " -w%6.4f", nonISMRelMutRate);
    fprintf(fp, " -f%3.2f %3.2f %3.2f %3.2f", freq[0], freq[1], freq[2],
            freq[3]);
    fprintf(fp, " -t%6.4f", titv);
    fprintf(fp, " -a%6.4f", alphaSites);
    fprintf(fp,
            " -r%3.2f %3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f "
            "%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f",
            Mij[0][0], Mij[0][1], Mij[0][2], Mij[0][3], Mij[1][0], Mij[1][1],
            Mij[1][2], Mij[1][3], Mij[2][0], Mij[2][1], Mij[2][2], Mij[2][3],
            Mij[3][0], Mij[3][1], Mij[3][2], Mij[3][3]);

    /* NGS */
    fprintf(fp, " -G%2.1e", genotypingError);
    fprintf(fp, " -C%d", coverage);
    fprintf(fp, " -V%6.4f", alphaCoverage);
    fprintf(fp, " -I%2.1e", allelicImbalance);
    fprintf(fp, " -D%2.1e", ADOrate);
    fprintf(fp, " -P%2.1e", alphaADOsites);
    fprintf(fp, " -Q%2.1e", alphaADOcells);
    fprintf(fp, " -R%2.1e", haploidCoverageReduction);
    fprintf(fp, " -A%2.1e %2.1e %d", meanAmplificationError,
            varAmplificationError, simulateOnlyTwoTemplates);
    fprintf(fp, " -E%2.1e", sequencingError);
    fprintf (fp," -B%2.1e,%2.1e", meanDoubletRate, varDoubletRate);
    fprintf(fp,
            " -X%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f "
            "%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f",
            Eij[0][0], Eij[0][1], Eij[0][2], Eij[0][3], Eij[1][0], Eij[1][1],
            Eij[1][2], Eij[1][3], Eij[2][0], Eij[2][1], Eij[2][2], Eij[2][3],
            Eij[3][0], Eij[3][1], Eij[3][2], Eij[3][3]);

    /* Output */
    if (doPrintSNVgenotypes == YES)
      fprintf(fp, " -%d", 1);
    if (doPrintSNVhaplotypes == YES)
      fprintf(fp, " -%d", 2);
    if (doPrintFullGenotypes == YES)
      fprintf(fp, " -%d", 3);
    if (doPrintFullHaplotypes == YES)
      fprintf(fp, " -%d", 4);
    if (doPrintAncestors == YES)
      fprintf(fp, " -%d", 5);
    if (doPrintTree == YES)
      fprintf(fp, " -%d", 6);
    if (doPrintTimes == YES)
      fprintf(fp, " -%d", 7);
    if (doPrintCATG == YES)
      fprintf(fp, " -%d", 8);
    if (doPrintTrueHaplotypes == YES)
      fprintf(fp, " -%d", 9);
    if (doPrintSeparateReplicates == YES)
      fprintf(fp, " -%c", 'v');
    if (doPrintIUPAChaplotypes == YES)
      fprintf(fp, " -%c", 'x');
    fprintf(fp, " -o%s", resultsDir);
    if (doUserTree == YES)
      fprintf(fp, " -T%s", userTreeFile);
    if (doUserGenome == YES)
      fprintf(fp, " -U%s", userGenomeFile);

    /* Other */
    if (doTumorNames == YES)
      fprintf(fp, " -%c", 'W');
    fprintf(fp, " -y%d -#%ld", noisy, originalSeed);
  }
}

/***************************** PrintUsage *******************************/
/* Prints a short description of program usage  */
void PrintUsage(FILE *fp) {
  PrintHeader(fp);
  fprintf(fp,
          "\n%s generates a coalescent tree and simulates a sample of diploid "
          "genomes from somatic cells (no recombination), together with a "
          "(healthy) cell as outgroup.",
          PROGRAM_NAME);
  fprintf(fp,
          "\n\nUsage: %s [-FTXT -n# -s# -l# -e# -g# -h# (# # #) -k# -q# -i# "
          "-b# -c# -u# -d# -H# -j# -S# (# #) -m# -p# -w# -f# # # # -t# -a# -r# "
          "# # # # # # # # # # # # # # # -G# -C# -V# -I# -D# -P# -Q# -R# -A# # "
          "# -E# -B# -X# # # # # # # # # # # # # # # # -1 -2 -3 -4 -5 -6 -7 -8 "
          "-9 -v -x -oTXT -TTXT -UTXT -W# -y# -## -?]",
          PROGRAM_NAME_LOWERCASE);
  fprintf(fp, "\n");

  fprintf(fp, "\n-F: user-defined parameter file (e.g. -Fmyparameters)");
  fprintf(fp, "\n-n: number of replicates (e.g. -n1000)");
  fprintf(fp, "\n-s: sample size (# cells) (e.g. -s8)");
  fprintf(fp, "\n-l: number of sites (e.g. -l500)");
  fprintf(fp, "\n-e: effective population size (e.g. -e1000)");
  fprintf(fp, "\n-g: exponential growth rate (e.g. -g1e-5)");
  fprintf(fp,
          "\n-h: number of demographic periods followed by effective  "
          "population size at the beginning and at the end of the period, ");
  fprintf(fp, "\n    and the duration of the period in generations. (e.g. -d1 "
              "100 200 30000) (e.g. -d2 100 100 40000 1000 2000 20000)");
  fprintf(fp, "\n-k: transforming branch length (e.g. -u1e-2)");
  fprintf(fp, "\n-q: healthy tip branch length (e.g. -q1e-6)");
  fprintf(fp, "\n-i: shape of the gamma distribution for rate variation among "
              "lineages (e.g. -i0.5)");
  fprintf(fp, "\n-b: alphabet [0:binary 1:DNA] (e.g. -b0)");
  fprintf(fp, "\n-c: germline SNP rate (e.g. -c1e-5)");
  fprintf(fp, "\n-u: somatic mutation rate (e.g. -u1e-6)");
  fprintf(fp, "\n-d: deletion rate (e.g. -d1e-5)");
  fprintf(fp, "\n-H: CNLOH rate (e.g. -H1e-5)");
  fprintf(fp, "\n-j: fixed number of mutations (e.g. -j100)");
  fprintf(fp, "\n-S: trinucleotide genetic signatures (e.g. -S2 3 0.8 13 0.2)");
  fprintf(fp, "\n-m: alternative mutation model [0:ISMhap 1:Mk2 2:finiteDNA] "
              "[default mutation model is ISM diploid] (e.g. -m2)");
  fprintf(fp, "\n-p: proportion of alternative model sites (e.g. -p0.1)");
  fprintf(fp,
          "\n-w: alternative/default model relative mutation rate (e.g. -w1)");
  fprintf(fp, "\n-f: nucleotide base frequencies (e.g. -f0.4 0.3 0.2 0.1)");
  fprintf(fp, "\n-t: transition/transversion ratio (e.g. -t2)");
  fprintf(fp, "\n-a: shape of the gamma distribution for rate variation among "
              "sites (e.g. -a0.2)");
  fprintf(fp, "\n-r: mutation matrix ACGT x ACGT (e.g. -r0 1 2 3 1 0 4 5 2 4 0 "
              "1 3 5 1 0)");
  fprintf(fp, "\n-G: genotyping error [no read counts] (e.g. -G0.01)");
  fprintf(fp, "\n-C: sequencing coverage [read counts] (e.g. -C60)");
  fprintf(fp, "\n-V: sequencing coverage overdispersion (e.g. -V5)");
  fprintf(fp, "\n-I: maternal allelic imbalance (e.g. -I0.5)");
  fprintf(fp, "\n-D: allelic dropout (e.g. -D0.1)");
  fprintf(fp, "\n-P: allelic dropout variation among sites (e.g. -P1)");
  fprintf(fp, "\n-Q: allelic dropout variation among cells (e.g. -Q1)");
  fprintf(fp, "\n-R: ADO/deletion haploid coverage reduction (e.g. -R0.5)");
  fprintf(fp, "\n-A: amplification error (e.g. -A0.1 0.01 0)");
  fprintf(fp, "\n-E: sequencing error (e.g. -E0.001)");
  fprintf(fp, "\n-B: doublet rate per cell (e.g. -B0.1)");
  fprintf(fp, "\n-X: error matrix ACGT x ACGT (e.g. -X0 1 1 1 1 0 1 1 1 1 0 1 "
              "1 1 1 0)");

  fprintf(fp, "\n-1: print SNV genotypes to a file (e.g. -1)");
  fprintf(fp, "\n-2: print SNV haplotypes to a file (e.g. -2)");
  fprintf(fp, "\n-3: print full genotypes to a file (e.g. -3)");
  fprintf(fp, "\n-4: print full haplotypes to a file (e.g. -4)");
  fprintf(fp, "\n-5: print ancestral genotypes (e.g. -5)");
  fprintf(fp, "\n-6: print trees to a file (e.g. -6)");
  fprintf(fp, "\n-7: print times to a file (e.g. -7)");
  fprintf(fp, "\n-8: print read counts in CATG format (e.g. -8)");
  fprintf(fp, "\n-9: print true SNV haplotypes to a file (e.g. -9)");
  fprintf(fp, "\n-v: print replicates in individual folders (e.g. -v)");
  fprintf(fp, "\n-x: print consensus/IUPAC haplotypes (e.g. -x)");
  fprintf(fp, "\n-o: results folder name (e.g. -oresultsFolder)");
  fprintf(fp, "\n-T: user-tree file (e.g. -Tusertree)");
  fprintf(fp, "\n-U: user genome (e.g. -Uusergenome)");
  fprintf(fp, "\n-W: use tumor nomenclature for cells (e.g. -W)");
  fprintf(fp, "\n-y: noisy (e.g. -y1)");
  fprintf(fp, "\n      = 0: does not print anything");
  fprintf(fp, "\n      = 1:  + simulation summary");
  fprintf(fp, "\n      = 2:  + replicate information");
  fprintf(fp, "\n      = 3: + calculation status and event information");
  fprintf(fp, "\n-#: seed (e.g. -#37864287)");
  fprintf(fp, "\n-?: Print help");

  fprintf(fp, "\n\nDefaults: ");
  PrintDefaults(fp);
  fprintf(fp, "\n--------------------------------------------------------------"
              "------------------------------------------\n");
  fprintf(fp, "\n...now exiting to system...\n\n");
  exit(-1);
}

/***************** PrintDefaults **********************/
/*	Prints the default settings */

static void PrintDefaults(FILE *fp) {
  int i;

  fprintf(fp, "\n-F: parameter file =  %s", "parameters");

  /* Coalescent */
  fprintf(fp, "\n-n: number of replicates =  %d", numDataSets);
  fprintf(fp, "\n-s: sample size (# cells) =  %d", numCells);
  fprintf(fp, "\n-l: number of sites =  %d", numSites);
  fprintf(fp, "\n-e: effective population size =  %d", N);
  fprintf(fp, "\n-g: exponential growth rate =  %2.1e", growthRate);
  fprintf(fp, "\n-h: number of demographic periods = %d", numPeriods);
  for (i = 1; i <= numPeriods; i++)
    fprintf(fp, "\n  period %d =  %d %d %d", i, Nbegin[i], Nend[i],
            cumDuration[i] - cumDuration[i - 1]);

  /* Post-coalescent */
  fprintf(fp, "\n-k: transforming branch length ratio =  %2.1e",
          transformingBranchLengthRatio);
  fprintf(fp, "\n-q: healthy tip branch length ratio =  %2.1e",
          healthyTipBranchLengthRatio);
  fprintf(fp,
          "\n-i: shape of the gamma distribution for rate variation among "
          "lineages =  %6.4f",
          alphaBranches);

  /* Mutation model */
  fprintf(fp, "\n-b: alphabet [0:binary 1:DNA] =  %d", alphabet);
  fprintf(fp, "\n-c: germline SNP rate =  %2.1e", SNPrate);
  fprintf(fp, "\n-u: mutation rate =  %2.1e", mutationRate);
  fprintf(fp, "\n-d: deletion rate =  %2.1e", deletionRate);
  fprintf(fp, "\n-H: CNLOH rate =  %2.1e", CNLOHrate);
  fprintf(fp, "\n-j: fixed number of mutations =  %d", numFixedMutations);
  fprintf(fp, "\n-S: number of trinucleotide genetic signature =  %d",
          numUserSignatures);
  fprintf(fp, "\n-m: alternative mutation model =  %d", altModel);
  fprintf(fp, "\n-p: proportion of alternative model sites =  %6.4f",
          propAltModelSites);
  fprintf(fp, "\n-w: alternative/default model relative mutation rate =  %6.4f",
          nonISMRelMutRate);
  fprintf(fp, "\n-f: nucleotide base frequencies = f%3.2f %3.2f %3.2f %3.2f",
          freq[0], freq[1], freq[2], freq[3]);
  fprintf(fp, "\n-t: transition/transversion ratio =  %6.4f", titv);
  fprintf(fp,
          "\n-a: shape of the gamma distribution for rate variation among "
          "sites =  %6.4f",
          alphaSites);
  fprintf(fp,
          "\n-r: mutation matrix ACGT x ACGT = %3.2f %3.2f %3.2f %3.2f  %3.2f "
          "%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f",
          Mij[0][0], Mij[0][1], Mij[0][2], Mij[0][3], Mij[1][0], Mij[1][1],
          Mij[1][2], Mij[1][3], Mij[2][0], Mij[2][1], Mij[2][2], Mij[2][3],
          Mij[3][0], Mij[3][1], Mij[3][2], Mij[3][3]);

  /* NGS */
  fprintf(fp, "\n-G: genotyping error =  %2.1e", genotypingError);
  fprintf(fp, "\n-C: sequencing coverage =  %d", coverage);
  fprintf(fp, "\n-V: sequencing coverage overdispersion =  %6.4f",
          alphaCoverage);
  fprintf(fp, "\n-I: allelic imbalance =  %2.1e", allelicImbalance);
  fprintf(fp, "\n-D: allelic dropout =  %2.1e", ADOrate);
  fprintf(fp, "\n-P: allelic dropout variation among sites =  %2.1e",
          ADOvarAmongSites);
  fprintf(fp, "\n-Q: allelic dropout cariation among cells=  %2.1e",
          ADOvarAmongCells);
  fprintf(fp, "\n-R: haploid coverage reduction =  %2.1e",
          haploidCoverageReduction);
  fprintf(
      fp,
      "\n-A: amplification error, var, 2-template model =  %2.1e, %2.1e, %d)",
      meanAmplificationError, varAmplificationError, simulateOnlyTwoTemplates);
  fprintf(fp, "\n-E: sequencing error =  %2.1e", sequencingError);
  fprintf (fp,"\n-B: doublet rate per cell; mean, var =  %2.1e, %2.1e", meanDoubletRate, varDoubletRate);
  fprintf(fp,
          "\n-X: error matrix ACGT x ACGT = %3.2f %3.2f %3.2f %3.2f  %3.2f "
          "%3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f",
          Eij[0][0], Eij[0][1], Eij[0][2], Eij[0][3], Eij[1][0], Eij[1][1],
          Eij[1][2], Eij[1][3], Eij[2][0], Eij[2][1], Eij[2][2], Eij[2][3],
          Eij[3][0], Eij[3][1], Eij[3][2], Eij[3][3]);

  /* Output */
  fprintf(fp, "\n-1: print SNV genotypes to a file =  %d", doPrintSNVgenotypes);
  fprintf(fp, "\n-2: print SNV haplotypes to a file =  %d",
          doPrintSNVhaplotypes);
  fprintf(fp, "\n-3: print full genotypes to a file =  %d",
          doPrintFullGenotypes);
  fprintf(fp, "\n-4: print full haplotypes to a file =  %d",
          doPrintFullHaplotypes);
  fprintf(fp, "\n-5: print ancestral genotypes =  %d", doPrintAncestors);
  fprintf(fp, "\n-6: print trees to a file =  %d", doPrintTree);
  fprintf(fp, "\n-7: print times to a file =  %d", doPrintTimes);
  fprintf(fp, "\n-8: print read counts in CATG format =  %d", doPrintCATG);
  fprintf(fp, "\n-9: print true haplotypes to a file =  %d",
          doPrintTrueHaplotypes);
  fprintf(fp, "\n-v: print replicates in individual folders =  %d",
          doPrintSeparateReplicates);
  fprintf(fp, "\n-x: print consensus/IUPAC haplotypes =  %d",
          doPrintIUPAChaplotypes);
  if (strlen(userTreeFile) == 0)
    strcpy(userTreeFile, "none");
  fprintf(fp, "\n-T: user tree file name =  %s", userTreeFile);
  if (strlen(userGenomeFile) == 0)
    strcpy(userGenomeFile, "none");
  fprintf(fp, "\n-U: user genome file name =  %s", userGenomeFile);
  if (strlen(resultsDir) == 0)
    strcpy(resultsDir, "Results");
  fprintf(fp, "\n-o: results folder name =  %s", resultsDir);

  /* Other */
  fprintf(fp, "\n-W: use tumor nomenclature for cells =  %d", doTumorNames);
  fprintf(fp, "\n-y: noisy = %d", noisy);
  fprintf(fp, "\n-#: seed = %ld", originalSeed);
}

/***************** Index ***************/
/* Returns index for a given node */

int Index(TreeNode *p) { return (p == NULL) ? -1 : p->index + 1; }

/***************** Label ***************/
/* Returns label for a given node */

int Label(TreeNode *p) {
  return (p->anc == NULL && p->left == NULL && p->right == NULL) ? -1
                                                                 : p->label + 1;
}

/************************************* ChooseTrinucMutation
 * **********************************************/
/* Chooses uniformy a random trinucleotide mutation (out of 96) for a given
 * signature 3D matrix with 96 types */

int ChooseTrinucMutation(double ***prob, RandomDevice *rd) {
  int i, j = 0, k = 0;
  double ran, cumProb;

  cumProb = 0;
  ran = RandUniform(rd);
  for (i = 0; i < 4; i++)
    for (j = 0; j < 6; j++)
      for (k = 0; k < 4; k++) {
        cumProb += prob[i][j][k];
        if (ran <= cumProb)
          break;
      }
  return trimut(i, j, k);
}

/************************************* ChooseUniformState
 * **********************************************/
/* Chooses uniformy a random state according to a vector of state
 * probabilities
 */

int ChooseUniformState(double *prob, RandomDevice *rd) {
  int chosenState;
  double ran, cumProb;

  chosenState = 0;
  cumProb = prob[chosenState];
  ran = RandUniform(rd);

  while (ran > cumProb)
    cumProb += prob[++chosenState];

  return chosenState;
}

/***************************** RandomUniform
 * **********************************/
/*  It returns a random uniform variate in range 0..1. It is described in
 Park, S. K. and K. W. Miller.  1988.  Random number generators: good
 ones are hard to find.  Communications of the ACM, 31(10):1192-1201.
 */

// double RandomUniform(long int *seed) {
//   long int lo, hi, test;
//
//   hi = (*seed) / 127773;
//   lo = (*seed) % 127773;
//   test = 16807 * lo - 2836 * hi;
//   if (test > 0)
//     *seed = test;
//   else
//     *seed = test + 2147483647;
//   return (double)(*seed) / (double)2147483647;
// }

/********************** RandomUniformTo ****************************/
/* it returns random uniform in range 0...max-1          */

// static int RandomUniformTo(int max, long int *seed) {
//   double rd;
//   rd = RandomUniform(seed);
//   return (floor(rd * max));
// }

/********************* RandomPoisson ********************/
/* Generates a random number from a Poisson distibution with
 mean lambda.
 */

// int RandomPoisson(double lambda, long int *seed) {
//   int poissonNumber;
//   double sum;
//
//   sum = 0;
//   poissonNumber = -1;
//
//   while (sum <= 1.0) {
//     sum += RandomExponential(lambda, seed);
//     poissonNumber++;
//   }
//
//   return poissonNumber;
// }

/********************* RandomExponential ********************/
/* Generates a random number from a Exponential distibution with
 mean lambda.
 */

// double RandomExponential(double lambda, long int *seed) {
//
//   double exponentialNumber, U;
//
//   do
//     U = RandomUniform(seed);
//   while (U == 0);
//
//   exponentialNumber = -log(U) / lambda;
//
//   return exponentialNumber;
// }

/********************* RandomBinomial ********************/
/* Generates a random number from a Binomial distibution with sucess
   probabilty p and n trials using the direct method (sum of Bernoulli
   variables).
 */

// int RandomBinomial(double prob, int numTrials, long int *seed) {
//   int i, sum;
//   sum = 0;
//
//   for (i = 0; i < numTrials; i++) {
//     if (RandomUniform(seed) < prob)
//       sum++;
//   }
//   return sum;
// }

/********************* RandomNegativeBinomial ********************/
/*
 *    Random variates from the negative binomial distribution.
 *
 *  NOTES
 *
 *    x = the number of failures before the n-th success
 *
 *  REFERENCE
 *
 *    Devroye, L. (1986).
 *    Non-Uniform Random Variate Generation.
 *    New York:Springer-Verlag.  Pages 488 and 543.
 *
 *  METHOD
 *
 *    Generate lambda as gamma with shape parameter "dispersion" (aka size)
and scale
 *    parameter "mean/dispersion".  Return a Poisson deviate with mean lambda.

**** NOTE: Extracted from rnbinom.c R code:
        rpois(rgamma(size, (1 - prob) / prob));
        rpois(rgamma(size, mu / size));

The negative binomial distribution with dispersion = n and prob = p has
density

p(x) = Gamma(x+n)/(Gamma(n) x!) p^n (1-p)^x

for x = 0, 1, 2, ..., n > 0 and 0 < p <= 1.

A negative binomial distribution can arise as a mixture of Poisson
distributions with mean distributed as a Γ (pgamma) distribution with scale
parameter (1 - prob)/prob and shape parameter dispersion. In this model prob =
scale/(1+scale), and the mean is dispersion * (1 - prob)/prob. The variance in
this parametrization is n (1-p)/p^2.

The alternative parameterization, often used in ecology, and the one used
here, is by the mean mu, and the dispersion parameter, where prob =
dispersion/(dispersion+mu). The variance is mu + mu^2/dispersion in this
parametrization.
*/

// int RandomNegativeBinomial(double mean, double dispersion, long int *seed) {
//   int poissonRand;
//   double gammaRand;
//
//   /* the RandomGamma function here has mean 1, so we need to scale it
//    * ourselves
//    */
//   do {
//     gammaRand = mean / dispersion * RandomGamma(dispersion, seed);
//   } while (fabs(gammaRand - 0) < 1e-16);
//   poissonRand = RandomPoisson(gammaRand, seed);
//
//   return poissonRand;
// }

/**************************** RandomBeta *************************/
/*	Generates a beta random number given mean and variance

An algorithm for generating beta variates B(α,β) is to generate X/(X + Y),
where X is a gamma variate with parameters (α, 1) and Y is an independent
gamma variate with parameters (β, 1).[53] From Numerical Recipes 3rd Edition:
The Art of Scientific Computing. 2007. Press et al. ISBN-13: 978-0521880688

note that the variance has to be > 0
*/

// double RandomBeta(double mean, double var, long int *seed) {
//   double sample_size, shape1, shape2, gamma1, gamma2, randBeta;
//
//   /* assuming variance < mean (1-mean) */
//   sample_size = (mean * (1.0 - mean) / var) - 1.0;
//
//   shape1 = mean * sample_size;
//   shape2 = (1.0 - mean) * sample_size;
//
//   gamma1 = RandomGamma(shape1, seed);
//   gamma2 = RandomGamma(shape2, seed);
//   randBeta = gamma1 / (gamma1 + gamma2);
//
//   return randBeta;
// }

/* Gamma functions are taken from Ziheng Yang's PAML package.
 http://abacus.gene.ucl.ac.uk/*/

/**************************** RandomGamma *************************/
/*	Generates a gamma number using routines in Ziheng's
        Yang tools.h in PAML

        Random standard gamma (Mean=Var=s,  with shape par=s, scale par=1)
        r^(s-1)*exp(-r)

        J. Dagpunar (1988) Principles of random variate generation,
        Clarendon Press, Oxford

        Calling rndgamma1() if s<1 or rndgamma2() if s>1 or exponential if s=1
 */

// double RandomGamma(double shape, long int *seed) {
//   double gammaNumber = 0;
//
//   if (shape <= 0)
//     fprintf(stderr,
//             "!!! ERROR: problems with gamma variable generation, shape < 0");
//   else if (shape < 1)
//     gammaNumber = RandomGamma1(shape, seed);
//   else if (shape > 1)
//     gammaNumber = RandomGamma2(shape, seed);
//   else
//     gammaNumber = -log(RandomUniform(seed));
//
//   return (gammaNumber);
// }

/*************** RandomGamma1 ***************/
// double RandomGamma1(double s, long int *seed) {
//   /* Random standard gamma for s<1
//    switching method
//    */
//   double r, x = 0.0, small = 1e-37, w;
//   static double a, p, uf, ss = 10.0, d;
//
//   if (s != ss) {
//     a = 1.0 - s;
//     p = a / (a + s * exp(-a));
//     uf = p * pow(small / a, s);
//     d = a * log(a);
//     ss = s;
//   }
//   for (;;) {
//     r = RandomUniform(seed);
//     if (r > p) {
//       x = a - log((1.0 - r) / (1.0 - p));
//       w = a * log(x) - d; /* this was with comma in line above before
//       270917*/
//     } else if (r > uf) {
//       x = a * pow(r / p, 1 / s);
//       w = x; /* this was with comma in line above before 270917*/
//     } else
//       return (0.0);
//     r = RandomUniform(seed);
//     if (1.0 - r <= w && r > 0.0)
//       if (r * (w + 1.0) >= 1.0 || -log(r) <= w)
//         continue;
//     break;
//   }
//   return (x);
// }

/*************** RandomGamma2 ***************/
// double RandomGamma2(double s, long int *seed) {
//   /* Random standard gamma for s>1
//    Best's (1978) t distribution method
//    */
//   double r, d, f, g, x;
//   static double b, h, ss = 0;
//
//   if (s != ss) {
//     b = s - 1.0;
//     h = sqrt(3.0 * s - 0.75);
//     ss = s;
//   }
//   for (;;) {
//     r = RandomUniform(seed);
//     g = r - r * r;
//     f = (r - 0.5) * h / sqrt(g);
//     x = b + f;
//     if (x <= 0.0)
//       continue;
//     r = RandomUniform(seed);
//     d = 64 * r * r * g * g * g;
//     if (d * x < x - 2.0 * f * f || log(d) < 2 * (b * log(x / b) - f))
//       break;
//   }
//   return (x);
// }

/************************ ReadUntil **************************/
/* Utility to skip comments between brackets. Borrowed from Andy
 Rambaut and Nick Grassly */

void ReadUntil(FILE *fv, char stopChar, char *what) {
  char ch;

  ch = fgetc(fv);
  while (!feof(fv) && ch != stopChar)
    ch = fgetc(fv);

  if (feof(fv) || ch != stopChar) {
    fprintf(stderr, "%s missing", what);
    exit(0);
  }
}

/************************ Unlink_callback **************************/
/* Call back to remove the passed path .
From
http://stackoverflow.com/questions/5467725/how-to-delete-a-directory-and-its-contents-in-posix-c
*/
int Unlink_callback(const char *fpath, const struct stat *sb, int typeflag,
                    struct FTW *ftwbuf) {
  int rv = remove(fpath);

  if (rv)
    perror(fpath);

  return rv;
}

/************************ RemoveDir **************************/
/* Removes directory and its contents
From
http://stackoverflow.com/questions/5467725/how-to-delete-a-directory-and-its-contents-in-posix-c
*/
int RemoveDir(char *path) {
  return nftw(path, Unlink_callback, 64, FTW_DEPTH | FTW_PHYS);
}

/************************ CheckMatrixSymmetry **************************/
/* Checks whether a given matrix is symmetric */

int CheckMatrixSymmetry(double matrix[4][4]) {
  int i, j;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      if (matrix[i][j] != matrix[j][i])
        return NO;
  return YES;
}

/******************** ReadParametersFromCommandLine
 * **************************/
/*
USED IN ORDER
F n s l e g h k q i b u d H j S m p w c f t a r G C V A E D P Q I R B M 1 2 3
4 5 6 7 8 9 v x o T U  W y  #

USED
a b c d e f g h i j k l m n o p q r s t u v w x y
A B C D E F G H I       M     P Q R S T U V W X
1 2 3 4 5 6 7 8 9 #
*/

static void ReadParametersFromCommandLine(int argc, char **argv) {
  int i, j;
  char flag, ch, *cat;
  double sumPi, sum;
  float argument;
  double argdouble;

  if (argc > 0) {
    cat = (char *)calloc(300, sizeof(char));
    CommandLine = (char *)calloc(MAX_LINE, sizeof(char));

    for (i = 0; i < argc; i++) {
      sprintf(cat, " %s", argv[i]);
      strcat(CommandLine, cat);
    }
    free(cat);
    cat = NULL;
  }

  for (i = 1; i < argc; i++) {
    argv[i]++;
    flag = *argv[i];
    argv[i]++;
    argument = -9999;

    switch (flag) {
    /* specifying a specific parameter file */
    case 'F':
      ch = *argv[i];
      if (!isspace(ch))
        strcpy(parameterFile, argv[i]);
      i = argc;
      doSpecificParameterFile = YES;
      break;
    case 'n':
      argument = atof(argv[i]);
      numDataSets = (int)argument;
      if (numDataSets < 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad number of replicates (%d)\n\n",
                numDataSets);
        PrintUsage(stderr);
      }
      break;
    case 's':
      argument = atof(argv[i]);
      numCells = (int)argument;
      if (numCells < 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad sample size (%d)\n\n",
                numCells);
        PrintUsage(stderr);
      }
      break;
    case 'l':
      argument = atof(argv[i]);
      numSites = (int)argument;
      if (numSites < 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad sequence length (%d)\n\n",
                numSites);
        PrintUsage(stderr);
      }
      break;
    case 'e':
      argument = atof(argv[i]);
      N = (int)argument;
      if (N < 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad effective population size (%d)\n\n",
                N);
        PrintUsage(stderr);
      }
      break;
    case 'g':
      growthRate = atof(argv[i]);
      if (growthRate < 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad growth rate (%f)\n\n",
                growthRate);
        PrintUsage(stderr);
      }
      if (growthRate != 0) {
        doExponential = YES;
        if (doDemographics == YES) {
          fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot have both exponential "
                          "(-g) and  demographics (-h)\n\n");
          exit(-1);
        }
      }
      break;
    case 'h':
      argument = atof(argv[i]);
      numPeriods = (int)argument;
      if (numPeriods > 0)
        doDemographics = YES;
      if (doDemographics == YES && doExponential == YES) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot have both demographics "
                        "periods (-h) and other demographics (-d)\n\n");
        exit(-1);
      }
      if (numPeriods <= 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad number of periods (%d)\n\n",
                numPeriods);
        PrintUsage(stderr);
      }
      Nbegin = (int *)calloc(numPeriods + 1, sizeof(int));
      Nend = (int *)calloc(numPeriods + 1, sizeof(int));
      cumDuration = (int *)calloc(numPeriods + 1, sizeof(int));
      periodGrowth = (double *)calloc(numPeriods + 1, sizeof(double));
      if (Nbegin == NULL || Nend == NULL || cumDuration == NULL) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Could not allocate demographic vectors "
                "(%lu)\n",
                numPeriods * sizeof(int));
        exit(-1);
      }
      for (j = 1; j <= numPeriods; j++) {
        argument = atof(argv[++i]);
        Nbegin[j] = (int)argument;
        argument = atof(argv[++i]);
        Nend[j] = (int)argument;
        argument = atof(argv[++i]);
        cumDuration[j] = (int)argument + cumDuration[j - 1];
      }
      break;
    case 'k':
      transformingBranchLengthRatio = atof(argv[i]);
      if (transformingBranchLengthRatio < 0) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad transforming branch length (%f)\n\n",
            transformingBranchLengthRatio);
        PrintUsage(stderr);
      }
      break;
    case 'q':
      healthyTipBranchLengthRatio = atof(argv[i]);
      if (healthyTipBranchLengthRatio < 0) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad healthy tip branch length (%f)\n\n",
                healthyTipBranchLengthRatio);
        PrintUsage(stderr);
      }
      break;
    case 'i':
      alphaBranches = atof(argv[i]);
      if (alphaBranches <= 0) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad gamma alphaBranches shape (%f)\n\n",
                alphaBranches);
        PrintUsage(stderr);
      }
      rateVarAmongLineages = YES;
      break;
    case 'b':
      argument = atof(argv[i]);
      alphabet = (int)argument;
      if (alphabet < 0 || alphabet > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad alphabet (%d)\n\n",
                alphabet);
        PrintUsage(stderr);
      }
      break;
    case 'u':
      mutationRate = atof(argv[i]);
      if (mutationRate < 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mutation rate (%f)\n\n",
                mutationRate);
        PrintUsage(stderr);
      }
      break;
    case 'd':
      deletionRate = atof(argv[i]);
      if (deletionRate < 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad deletion rate (%f)\n\n",
                deletionRate);
        PrintUsage(stderr);
      }
      break;
    case 'H':
      CNLOHrate = atof(argv[i]);
      if (CNLOHrate < 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad CNLOH rate (%f)\n\n",
                CNLOHrate);
        PrintUsage(stderr);
      }
      break;
    case 'j':
      argument = atof(argv[i]);
      numFixedMutations = (int)argument;
      if (numFixedMutations < 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad number of fixed mutations (%d)\n\n",
                numFixedMutations);
        PrintUsage(stderr);
      }
      doSimulateFixedNumMutations = YES;
      if (propAltModelSites > 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: You cannot specify a fixed "
                        "number of mutations if there is any non-ISM site. Set "
                        "the proportion of non-ISM sites to zero\n\n");
        PrintUsage(stderr);
      }
      break;
    case 'S':
      argument = atof(argv[i]);
      numUserSignatures = (int)argument;
      if (numUserSignatures < 1 || numUserSignatures > 30) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad  number of genetic signatures "
                "(%d), it has to be a single number between 1 and 30\n\n",
                numUserSignatures);
        PrintUsage(stderr);
      }
      doGeneticSignatures = YES;

      if (alphabet == BINARY) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Trinucleotide signatures "
                        "cannot be simulated with the binary alphabet");
        PrintUsage(stderr);
      }

      signatureWeight = (double *)calloc(numUserSignatures, sizeof(double));
      if (!signatureWeight) {
        fprintf(stderr, "Could not allocate the signatureWeight array\n");
        exit(-1);
      }

      signatureID = (int *)calloc(numUserSignatures, sizeof(int));
      if (!signatureID) {
        fprintf(stderr, "Could not allocate the signatureID array\n");
        exit(-1);
      }
      sum = 0;
      for (j = 0; j < numUserSignatures; j++) {
        argument = atof(argv[++i]);
        signatureID[j] = (int)argument;
        argument = atof(argv[++i]);
        signatureWeight[j] = argument;
        if (signatureWeight[j] <= 0) {
          fprintf(stderr, "\n!!! PARAMETER ERROR: Bad signatureWeight (%f)\n\n",
                  signatureWeight[j]);
          PrintUsage(stderr);
        }
        sum += signatureWeight[j];
      }
      if (sum != 1.0)
        for (j = 0; j < numUserSignatures; j++)
          signatureWeight[j] /= sum;
      break;
    case 'm':
      argument = atof(argv[i]);
      altModel = (int)argument;
      if (altModel < 0 || altModel > 2) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad alternative mutation model (%d)\n\n",
            altModel);
        PrintUsage(stderr);
      }
      if (alphabet == DNA && propAltModelSites > 0) {
        if (altModel == Mk) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: The DNA alphabet and the model (%d) "
                  "specified are incompatible",
                  (int)argument);
          PrintUsage(stderr);
        }
      } else if (alphabet == BINARY && propAltModelSites > 0) {
        if (altModel == finiteDNA) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: The binary alphabet and the model "
                  "(%d) specified are incompatible",
                  (int)argument);
          PrintUsage(stderr);
        }
      }
      break;
    case 'p':
      propAltModelSites = atof(argv[i]);
      if (propAltModelSites < 0 || propAltModelSites > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad proportion of alternative model "
                "sites (%f). It has to be between 0 and 1\n\n",
                propAltModelSites);
        PrintUsage(stderr);
      }
      if (propAltModelSites > 0 && altModel != ISMhap &&
          doSimulateFixedNumMutations == YES) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: You cannot specify a "
                        "proportion of non-ISM "
                        " sites bigger than zero if the number of mutations "
                        "is fixed\n\n");
        PrintUsage(stderr);
      }
      if (alphabet == DNA && propAltModelSites > 0) {
        if (altModel == Mk) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: The DNA alphabet and the alt model "
                  "(%d) specified are incompatible",
                  (int)altModel);
          PrintUsage(stderr);
        }
      } else if (alphabet == BINARY && propAltModelSites > 0) {
        if (altModel == finiteDNA) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: The binary alphabet and the alt "
                  "model (%d) specified are incompatible",
                  (int)altModel);
          PrintUsage(stderr);
        }
      }
      break;
    case 'w':
      nonISMRelMutRate = atof(argv[i]);
      if (nonISMRelMutRate < 0 || nonISMRelMutRate > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad relative nonISM/ISM mutation rate "
                "(%f)\n\n",
                nonISMRelMutRate);
        PrintUsage(stderr);
      }
      break;
    case 'c':
      SNPrate = atof(argv[i]);
      if (SNPrate < 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad germline SNP rate (%f)\n\n",
                SNPrate);
        PrintUsage(stderr);
      }
      break;
    case 'f':
      freq[0] = atof(argv[i]);
      freq[1] = atof(argv[++i]);
      freq[2] = atof(argv[++i]);
      freq[3] = atof(argv[++i]);

      if ((freq[0] == freq[1]) && (freq[0] == freq[2]) &&
          (freq[0] == freq[3])) /* TODO fix with loop */
        equalBaseFreq = YES;
      else
        equalBaseFreq = NO;

      sumPi = freq[0] + freq[1] + freq[2] + freq[3];
      if (sumPi != 1.0) {
        freq[0] /= sumPi;
        freq[1] /= sumPi;
        freq[2] /= sumPi;
        freq[3] /= sumPi;
      }
      break;
    case 't':
      titv = atof(argv[i]);
      if (titv < 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad ti/tv (%f)\n\n", titv);
        PrintUsage(stderr);
      } else {
        doJC = NO;
        doHKY = YES;
        doGTR = NO;
        doGTnR = NO;
      }
      if (thereIsMij == YES) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Cannot specify a mutation matrix (GTR "
                "model) and a ti/tv (HKY model) at the same time\n\n");
        PrintUsage(stderr);
      }
      break;
    case 'a':
      alphaSites = atof(argv[i]);
      if (alphaSites < 0) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad gamma alphaSites shape (%f)\n\n",
                alphaSites);
        PrintUsage(stderr);
      }
      rateVarAmongSites = YES;
      break;
    case 'r':
      if (doHKY == YES) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Cannot specify a mutation matrix (GTR "
                "model) and a ti/tv (HKY model) at the same time\n\n");
        PrintUsage(stderr);
      }
      /*	AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT */
      Mij[0][0] = atof(argv[i]);
      Mij[0][1] = atof(argv[++i]);
      Mij[0][2] = atof(argv[++i]);
      Mij[0][3] = atof(argv[++i]);
      Mij[1][0] = atof(argv[++i]);
      Mij[1][1] = atof(argv[++i]);
      Mij[1][2] = atof(argv[++i]);
      Mij[1][3] = atof(argv[++i]);
      Mij[2][0] = atof(argv[++i]);
      Mij[2][1] = atof(argv[++i]);
      Mij[2][2] = atof(argv[++i]);
      Mij[2][3] = atof(argv[++i]);
      Mij[3][0] = atof(argv[++i]);
      Mij[3][1] = atof(argv[++i]);
      Mij[3][2] = atof(argv[++i]);
      Mij[3][3] = atof(argv[++i]);
      if (Mij[0][0] != 0 || Mij[1][1] != 0 || Mij[2][2] != 0 ||
          Mij[3][3] != 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad general rate matrix: "
                        "diagonals should be 0 \n\n");
        PrintUsage(stderr);
      }
      thereIsMij = YES;
      if (CheckMatrixSymmetry(Mij) == YES) {
        doJC = NO;
        doHKY = NO;
        doGTR = YES;
        doGTnR = NO;
      } else {
        doJC = NO;
        doHKY = NO;
        doGTR = NO;
        doGTnR = YES;
      }
      break;
    case 'G':
      genotypingError = atof(argv[i]);
      if (genotypingError < 0 || genotypingError > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad genotyping error (%f)\n\n",
                genotypingError);
        PrintUsage(stderr);
      }
      if (genotypingError > 0 && doNGS == YES) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot specify a coverage "
                        "larger than 0, which implies read count generation, "
                        "and a genotyping error at the same time\n\n");
        PrintUsage(stderr);
      }
      break;
    case 'C':
      argument = atof(argv[i]);
      coverage = (int)argument;
      if (coverage < 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad coverage (%d)\n\n",
                coverage);
        PrintUsage(stderr);
      }
      if (coverage > 0)
        doNGS = YES;
      if (genotypingError > 0 && doNGS == YES) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot specify a coverage "
                        "larger than 0, which implies read count generation, "
                        "and a genotyping error at the same time\n\n");
        PrintUsage(stderr);
      }
      break;
    case 'V':
      alphaCoverage = atof(argv[i]);
      if (alphaCoverage <= 0) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad coverage dispersion (%f)\n\n",
                alphaCoverage);
        PrintUsage(stderr);
      }
      rateVarCoverage = YES;
      break;
    case 'A':
      meanAmplificationError = atof(argv[i]);
      varAmplificationError = atof(argv[++i]);
      simulateOnlyTwoTemplates = (int)atof(argv[++i]);
      if (meanAmplificationError < 0 || meanAmplificationError > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad mean amplification error (%f)\n\n",
                meanAmplificationError);
        PrintUsage(stderr);
      }
      if (varAmplificationError <= 0 ||
          (meanAmplificationError > 0 &&
           varAmplificationError >=
               (meanAmplificationError * (1.0 - meanAmplificationError)))) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad variance amplification error (%f); "
                "it has to be > 0 and < mean*(1-mean)\n\n",
                varAmplificationError);
        PrintUsage(stderr);
      }
      if (simulateOnlyTwoTemplates < 0 || simulateOnlyTwoTemplates > 1) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); "
            "it has to be 0 (assume 2 templates) or 1 (assume 4 templates)",
            simulateOnlyTwoTemplates);
        PrintUsage(stderr);
      }
      break;
    case 'E':
      sequencingError = atof(argv[i]);
      if (sequencingError < 0 || sequencingError > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad sequencing error (%f)\n\n",
                sequencingError);
        PrintUsage(stderr);
      }
      break;
    case 'D':
      ADOrate = atof(argv[i]);
      if (ADOrate < 0 || ADOrate > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad allelic dropout (%f)\n\n",
                ADOrate);
        PrintUsage(stderr);
      }
      break;
    case 'P':
      alphaADOsites = atof(argv[i]);
      ADOvarAmongSites = YES;
      if (alphaADOsites < 0 || alphaADOsites > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad allelic dropout heterogeneity "
                "among sites (%f)\n\n",
                alphaADOsites);
        PrintUsage(stderr);
      }
      break;
    case 'Q':
      alphaADOcells = atof(argv[i]);
      ADOvarAmongCells = YES;
      if (alphaADOcells < 0 || alphaADOcells > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad allelic dropout heterogeneity "
                "among cells (%f)\n\n",
                alphaADOcells);
        PrintUsage(stderr);
      }
      break;
    case 'I':
      allelicImbalance = atof(argv[i]);
      if (allelicImbalance < 0 || allelicImbalance > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad allelicImbalance (%f)\n\n",
                haploidCoverageReduction);
        PrintUsage(stderr);
      }
      break;
    case 'R':
      haploidCoverageReduction = atof(argv[i]);
      if (haploidCoverageReduction < 0 || haploidCoverageReduction > 1) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad haploid coverage reduction (%f)\n\n",
            haploidCoverageReduction);
        PrintUsage(stderr);
      }
      break;
    case 'B':
      meanDoubletRate = atof(argv[i]);
      varDoubletRate = atof(argv[++i]);
      if (meanDoubletRate < 0 || meanDoubletRate > 1)
      {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mean doublet rate (%f)\n\n", meanDoubletRate);
        PrintUsage(stderr);
      }
      if (varDoubletRate <= 0 || (meanDoubletRate > 0 && varDoubletRate >= (meanDoubletRate * (1.0 - meanDoubletRate))))
      {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad variance doublet rate (%f); it has to be > 0 and < mean*(1-mean)\n\n", varDoubletRate);
        PrintUsage(stderr);
      }
      break;
    case 'X':
      /*	AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT */
      Eij[0][0] = atof(argv[i]);
      Eij[0][1] = atof(argv[++i]);
      Eij[0][2] = atof(argv[++i]);
      Eij[0][3] = atof(argv[++i]);
      Eij[1][0] = atof(argv[++i]);
      Eij[1][1] = atof(argv[++i]);
      Eij[1][2] = atof(argv[++i]);
      Eij[1][3] = atof(argv[++i]);
      Eij[2][0] = atof(argv[++i]);
      Eij[2][1] = atof(argv[++i]);
      Eij[2][2] = atof(argv[++i]);
      Eij[2][3] = atof(argv[++i]);
      Eij[3][0] = atof(argv[++i]);
      Eij[3][1] = atof(argv[++i]);
      Eij[3][2] = atof(argv[++i]);
      Eij[3][3] = atof(argv[++i]);
      if (Eij[0][0] != 0 || Eij[1][1] != 0 || Eij[2][2] != 0 ||
          Eij[3][3] != 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad error matrix: diagonals "
                        "should be 0 \n\n");
        PrintUsage(stderr);
      }
      thereIsEij = YES;
      break;
    case '1':
      doPrintSNVgenotypes = YES;
      break;
    case '2':
      doPrintSNVhaplotypes = YES;
      break;
    case '3':
      doPrintFullGenotypes = YES;
      break;
    case '4':
      doPrintFullHaplotypes = YES;
      break;
    case '5':
      doPrintAncestors = YES;
      break;
    case '6':
      doPrintTree = YES;
      break;
    case '7':
      doPrintTimes = YES;
      break;
    case '8':
      doPrintCATG = YES;
      if (doNGS == NO) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot print the CATG format "
                        "when coverage is <= 0\n\n");
        PrintUsage(stderr);
      }
      break;
    case '9':
      doPrintTrueHaplotypes = YES;
      break;
    case 'v':
      doPrintSeparateReplicates = YES;
      break;
    case 'x':
      doPrintIUPAChaplotypes = YES;
      break;
    case 'o':
      ch = *argv[i];
      if (!isspace(ch))
        strcpy(resultsDir, argv[i]);
      break;
    case 'T':
      ch = *argv[i];
      if (!isspace(ch))
        strcpy(userTreeFile, argv[i]);
      doUserTree = YES;
      break;
    case 'U':
      ch = *argv[i];
      if (!isspace(ch))
        strcpy(userGenomeFile, argv[i]);
      doUserGenome = YES;
      break;
    /*case '0':
doSimulateData = NO;
break;*/
    case 'W':
      doTumorNames = YES;
      break;
    case 'y':
      argument = atof(argv[i]);
      noisy = (int)argument;
      if (noisy < 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad noisy value (%d)\n\n",
                noisy);
        PrintUsage(stderr);
      }
      break;
    /*case 'z':
        argument = atof(argv[i]);
        numNodes = (int) argument;
        if (numNodes<1)
            {
            fprintf (stderr, "\n!!! PARAMETER ERROR: Bad number of nodes
       (%d)\n\n", numNodes); PrintUsage(stderr);
            }
        break;*/
    case '#':
      argdouble = atof(argv[i]);
      userSeed = (long int)argdouble;
      if (userSeed < 0) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad user seed for random number "
                "(%ld)\n\n",
                userSeed);
        PrintUsage(stderr);
      }
      break;
    case '?':
      PrintUsage(stderr);
      break;
    default:
      fprintf(stderr, "\n!!! PARAMETER ERROR: Incorrect parameter: %c\n\n",
              flag);
      PrintUsage(stderr);
      break;
    }
  }
}

/***************************** ReadParametersFromFile
 * *******************************/
/*
 USED IN ORDER
 n s l e g h k q i b u d H j S m p w c f t a r G C V A E D P Q I R B M 1 2 3 4
 5 6 7 8 9 v x o T U  W y  #

 USED
 a b c d e f g h i j k l m n o p q r s t u v w x y
 A B C D E F G H I       M     P Q R S T U V W X
 1 2 3 4 5 6 7 8 9 #

 ADDED
 z (defined after d in order)
 J K L N O Y Z
 0 -
 */

void ReadParametersFromFile() {
  int j;
  char ch;
  double sumPi, sum;
  float argument;
  double argdouble;

  if (feof(stdin)) {
    fprintf(stderr, "FILE ERROR: Unable to read parameters from stdin\n");
    exit(-1);
  }

  ch = fgetc(stdin);
  while (isspace(ch))
    ch = fgetc(stdin);
  while (ch == '[') {
    ReadUntil(stdin, ']', "closing bracket");
    ch = fgetc(stdin);
    while (isspace(ch))
      ch = fgetc(stdin);
  }

  while (!feof(stdin)) {
    argument = 0;

    switch (ch) {
      /* added code */
    case '-':
      if (fscanf(stdin, "%d", &numThreads) != 1 || numThreads < 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad number of threads (%d) \n\n",
                numThreads);
        PrintUsage(stderr);
      }
      AllocateThreads(numThreads);
      break;
      /* end added code */
    case 'n':
      /* added code */
      if (fscanf(stdin, "%d %d", &numDataSets, &fixedTree) != 2 ||
          numDataSets < 1 || !(fixedTree == YES || fixedTree == NO)) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad number of replicates (%d) or fixed "
                "tree usage (%d) \n\n",
                numDataSets, fixedTree);
        PrintUsage(stderr);
      }
      /* end added code */

      /* original code */
      /*
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad number of replicates (%d)\n\n",
                (int)argument);
        PrintUsage(stderr);
      }
      numDataSets = (int)argument;
       */
      break;
    case 's':
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad sample size (%d)\n\n",
                (int)argument);
        PrintUsage(stderr);
      }
      numCells = (int)argument;
      break;
    case 'l':
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad sequence length (%d)\n\n",
                (int)argument);
        PrintUsage(stderr);
      }
      numSites = (int)argument;
      break;
    case 'e':
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad effective population size (%d)\n\n",
                (int)argument);
        PrintUsage(stderr);
      }
      N = (int)argument;
      break;
    case 'g':
      if (fscanf(stdin, "%lf", &growthRate) != 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad exponential growth rate (%f)\n\n",
                growthRate);
        PrintUsage(stderr);
      }
      if (growthRate != 0) {
        doExponential = YES;
        if (doDemographics == YES) {
          fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot have both exponential "
                          "(-b) and other demographics(-d)\n\n");
          exit(-1);
        }
      }
      break;
    case 'h':
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad number of periods (%d)\n\n",
                (int)argument);
        PrintUsage(stderr);
      }
      numPeriods = (int)argument;
      if (numPeriods > 0)
        doDemographics = YES;
      if (doDemographics == YES && doExponential == YES) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot have both demographics "
                        "periods (-d) and other demographics (-b)\n\n");
        exit(-1);
      }
      Nbegin = (int *)calloc(numPeriods + 1, sizeof(int));
      Nend = (int *)calloc(numPeriods + 1, sizeof(int));
      cumDuration = (int *)calloc(numPeriods + 1, sizeof(int));
      periodGrowth = (double *)calloc(numPeriods + 1, sizeof(double));
      if (Nbegin == NULL || Nend == NULL || cumDuration == NULL) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Could not allocate demographic vectors "
                "(%lu)\n",
                numPeriods * sizeof(int));
        exit(-1);
      }
      for (j = 1; j <= numPeriods; j++) {
        fscanf(stdin, "%f", &argument);
        Nbegin[j] = (int)argument;
        fscanf(stdin, "%f", &argument);
        Nend[j] = (int)argument;
        fscanf(stdin, "%f", &argument);
        cumDuration[j] = (int)argument + cumDuration[j - 1];
      }
      break;
    case 'k':
      if (fscanf(stdin, "%lf", &transformingBranchLengthRatio) != 1 ||
          transformingBranchLengthRatio < 0) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad transforming branch length (%f)\n\n",
            transformingBranchLengthRatio);
        PrintUsage(stderr);
      }
      break;
    case 'q':
      if (fscanf(stdin, "%lf", &healthyTipBranchLengthRatio) != 1 ||
          healthyTipBranchLengthRatio < 0) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad healthy tip branch length (%f)\n\n",
                healthyTipBranchLengthRatio);
        PrintUsage(stderr);
      }
      break;
    case 'i':
      if (fscanf(stdin, "%lf", &alphaBranches) != 1 || alphaBranches <= 0) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad alphaBranches shape (%f)\n\n",
                alphaBranches);
        PrintUsage(stderr);
      }
      rateVarAmongLineages = YES;
      break;
    case 'b':
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 0 || argument > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad alphabet (%d)\n\n",
                (int)argument);
        PrintUsage(stderr);
      }
      alphabet = (int)argument;
      break;
    case 'u':
      if (fscanf(stdin, "%lf", &mutationRate) != 1 || mutationRate < 0 ||
          mutationRate > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mutation rate (%f)\n\n",
                mutationRate);
        PrintUsage(stderr);
      }
      break;
    case 'd':
      if (fscanf(stdin, "%lf", &deletionRate) != 1 || deletionRate < 0 ||
          deletionRate > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad deletion rate (%f)\n\n",
                deletionRate);
        PrintUsage(stderr);
      }
      break;
    /* added code */
    case 'z':
      if (fscanf(stdin, "%lf", &insertionRate) != 1 || insertionRate < 0 ||
          insertionRate > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad copy and insertion rate (%f)\n\n",
                insertionRate);
        PrintUsage(stderr);
      }
      break;
    /* end added code */
    case 'H':
      if (fscanf(stdin, "%lf", &CNLOHrate) != 1 || CNLOHrate < 0 ||
          CNLOHrate > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad CNLOH rate (%f)\n\n",
                CNLOHrate);
        PrintUsage(stderr);
      }
      break;
    case 'j':
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad number of mutations (%d)\n\n",
                (int)argument);
        PrintUsage(stderr);
      }
      numFixedMutations = (int)argument;
      doSimulateFixedNumMutations = YES;
      if (propAltModelSites > 0 && altModel != ISMhap) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: You cannot specify a fixed number of "
                "mutations if there is any non-ISM  site. Set the proportion "
                "of non-ISM diploid sites to zero\n\n");
        PrintUsage(stderr);
      }
      break;
    case 'S':
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 1 ||
          argument > 30) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad number of genetic signatures (%d), "
                "it has to be a single number between 1 and 30\n\n",
                (int)argument);
        PrintUsage(stderr);
      }
      numUserSignatures = (int)argument;
      doGeneticSignatures = YES;

      if (alphabet == BINARY) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Trinucleotide signatures "
                        "cannot be simulated with the binary alphabet");
        PrintUsage(stderr);
      }

      signatureWeight = (double *)calloc(numUserSignatures, sizeof(double));
      if (!signatureWeight) {
        fprintf(stderr, "Could not allocate the signatureWeight array\n");
        exit(-1);
      }

      signatureID = (int *)calloc(numUserSignatures, sizeof(int));
      if (!signatureID) {
        fprintf(stderr, "Could not allocate the signatureID array\n");
        exit(-1);
      }

      sum = 0;
      for (j = 0; j < numUserSignatures; j++) {
        fscanf(stdin, "%f", &argument);
        signatureID[j] = (int)argument;
        fscanf(stdin, "%lf", &signatureWeight[j]);
        if (signatureWeight[j] <= 0) {
          fprintf(stderr, "\n!!! PARAMETER ERROR: Bad signatureWeight (%f)\n\n",
                  signatureWeight[j]);
          PrintUsage(stderr);
        }
        sum += signatureWeight[j];
      }
      if (sum != 1.0)
        for (j = 0; j < numUserSignatures; j++)
          signatureWeight[j] /= sum;
      break;
    case 'm':
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 0 || argument > 7) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad alternative mutation model (%d)\n\n",
            (int)argument);
        PrintUsage(stderr);
      }
      altModel = (int)argument;
      if (alphabet == DNA && propAltModelSites > 0) {
        if (altModel == Mk) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: The DNA alphabet and the model (%d) "
                  "specified are incompatible",
                  (int)argument);
          PrintUsage(stderr);
        }
      } else if (alphabet == BINARY && propAltModelSites > 0) {
        if (altModel == finiteDNA) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: The binary alphabet and the model "
                  "(%d) specified are incompatible",
                  (int)argument);
          PrintUsage(stderr);
        }
      }

      /* added code */
      if (altModel >= FiniteMu) {
        if (alphabet == BINARY && propAltModelSites > 0) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: The binary alphabet and the model "
                  "(%d) specified are incompatible",
                  (int)argument);
          PrintUsage(stderr);
        }

        if (altModel == 6) {
          constrained = YES;
          altModel = FiniteMuDel;
        } else if (altModel == 7) {
          constrained = YES;
          altModel = FiniteMuInDel;
        }

        if (altModel == FiniteMu) {
          InitRateMatrixMD(mutationRate, 0);
        } else if (altModel == FiniteMuDel) {
          InitRateMatrixMD(mutationRate, deletionRate);
        } else if (altModel == FiniteMuInDel) {
          InitRateMatrixMDI(mutationRate, deletionRate, insertionRate);
        }
      }
      /* end added code */

      break;
    case 'p':
      if (fscanf(stdin, "%lf", &propAltModelSites) != 1 ||
          propAltModelSites < 0 || propAltModelSites > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad proportion of alternative model "
                "sites (%f)\n\n",
                propAltModelSites);
        PrintUsage(stderr);
      }
      if (propAltModelSites > 0 && doSimulateFixedNumMutations == YES &&
          altModel != ISMhap) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: You cannot specify a "
                        "proportion of non-ISM diploid sites bigger than zero "
                        "if the number of mutations is fixed\n\n");
        PrintUsage(stderr);
      }
      if (alphabet == DNA && propAltModelSites > 0) {
        if (altModel == Mk) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: The DNA alphabet and the alt model "
                  "(%d) specified are incompatible",
                  (int)altModel);
          PrintUsage(stderr);
        }
      } else if (alphabet == BINARY && propAltModelSites > 0) {
        if (altModel == finiteDNA) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: The binary alphabet and the alt "
                  "model (%d) specified are incompatible",
                  (int)altModel);
          PrintUsage(stderr);
        }
      }
      break;
    case 'w':
      if (fscanf(stdin, "%lf", &nonISMRelMutRate) != 1 ||
          nonISMRelMutRate < 0 || nonISMRelMutRate > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad relative non-ISM/ISM mutation rate "
                "(%f)\n\n",
                nonISMRelMutRate);
        PrintUsage(stderr);
      }
      break;
    case 'c':
      if (fscanf(stdin, "%lf", &SNPrate) != 1 || SNPrate < 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad germline SNP rate (%f)\n\n",
                SNPrate);
        PrintUsage(stderr);
      }
      break;
    case 'f':
      if (fscanf(stdin, "%lf %lf %lf %lf", &freq[0], &freq[1], &freq[2],
                 &freq[3]) != 4) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad base frequencies (%lf %lf %lf %lf)\n",
            freq[0], freq[1], freq[2], freq[3]);
        PrintUsage(stderr);
      } else if (freq[0] == freq[1] == freq[2] == freq[3])
        equalBaseFreq = YES;
      else
        equalBaseFreq = NO;
      sumPi = freq[0] + freq[1] + freq[2] + freq[3];
      if (sumPi != 1.0) {
        freq[0] /= sumPi;
        freq[1] /= sumPi;
        freq[2] /= sumPi;
        freq[3] /= sumPi;
      }
      break;
    case 't':
      if (fscanf(stdin, "%lf", &titv) != 1 || titv <= 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad ti/tv (%f)\n\n", titv);
        PrintUsage(stderr);
      } else {
        doJC = NO;
        doHKY = YES;
        doGTR = NO;
        doGTnR = NO;
      }
      if (thereIsMij == YES) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Cannot specify a mutation matrix (GTR "
                "model) and a ti/tv (HKY model) at the same time\n\n");
        PrintUsage(stderr);
      }
      break;
    case 'a':
      if (fscanf(stdin, "%lf", &alphaSites) != 1 || alphaSites <= 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad alphaSites shape (%f)\n\n",
                alphaSites);
        PrintUsage(stderr);
      }
      rateVarAmongSites = YES;
      break;
    case 'r':
      if (doHKY == YES) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Cannot specify a mutation matrix (GTR "
                "model) and a ti/tv (HKY model) at the same time\n\n");
        PrintUsage(stderr);
      }
      if (fscanf(stdin,
                 "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf "
                 "%lf",
                 &Mij[0][0], &Mij[0][1], &Mij[0][2], &Mij[0][3], &Mij[1][0],
                 &Mij[1][1], &Mij[1][2], &Mij[1][3], &Mij[2][0], &Mij[2][1],
                 &Mij[2][2], &Mij[2][3], &Mij[3][0], &Mij[3][1], &Mij[3][2],
                 &Mij[3][3]) != 16) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad general rate matrix (-rx x x x x x x "
            "x x x x) (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT)\n\n");
        PrintUsage(stderr);
      }

      if (Mij[0][0] != 0 || Mij[1][1] != 0 || Mij[2][2] != 0 ||
          Mij[3][3] != 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad general rate matrix: "
                        "diagonals should be 0 \n\n");
        PrintUsage(stderr);
      }
      thereIsMij = YES;

      if (CheckMatrixSymmetry(Mij) == YES) {
        doJC = NO;
        doHKY = NO;
        doGTR = YES;
        doGTnR = NO;
      } else {
        doJC = NO;
        doHKY = NO;
        doGTR = NO;
        doGTnR = YES;
      }

      /* added code */
      if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu) {
        doJC = NO;
        doHKY = NO;
        doGTR = NO;
        doGTnR = NO;
        Mij[0][0] = Mij[1][1] = Mij[2][2] = Mij[3][3] = 0.0;
        Mij[0][1] = Mij[0][2] = Mij[0][3] = 0.0;
        Mij[1][0] = Mij[1][2] = Mij[1][3] = 0.0;
        Mij[2][0] = Mij[2][1] = Mij[2][3] = 0.0;
        Mij[3][0] = Mij[3][1] = Mij[3][2] = 0.0;

        fprintf(stderr, "\nWARNING! INVALID PARAMETER: Substitution rate "
                        "matrix has already been defined: "
                        "user-defined matrix will be overridden. \n\n");
      }
      /* end added code */

      break;
    case 'G':
      if (fscanf(stdin, "%lf", &genotypingError) != 1 || genotypingError < 0 ||
          genotypingError > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad genotyping error (%f)\n\n",
                genotypingError);
        PrintUsage(stderr);
      }
      if (genotypingError > 0 && doNGS == YES) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot specify a coverage "
                        "larger than 0, which implies read count generation, "
                        "and a genotyping error at the same time\n\n");
        PrintUsage(stderr);
      }
      break;
      /* added code */
    case 'O':
      if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu) {
        if (fscanf(stdin, "%d", &seqCovCellCoalStyle) != 1 || (seqCovCellCoalStyle != 0 && seqCovCellCoalStyle != 1)) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: Bad sequencing coverage simulation style (%d)\n\n",
                  seqCovCellCoalStyle);
          PrintUsage(stderr);
        }
      } else
      {
        fprintf(stderr,
                  "\n!!! PARAMETER ERROR: Unknown argument\n\n");
        PrintUsage(stderr);
      }
      break;
      /* end added code */
    case 'C':
      /* added code */
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        if (fscanf(stdin, "%f", &argument) != 1 || argument < 1) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: Bad sequencing coverage (%d)\n\n",
                  (int)coverage);
          PrintUsage(stderr);
        }
        coverage = (int)argument;
        if (coverage > 0)
          doNGS = YES;
      } else {
        if (fscanf(stdin, "%lf %lf", &allelicCovMean, &allelicCovVar) != 2 ||
            allelicCovMean < 1 || allelicCovVar < 0) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: Bad allelic sequencing "
                  "coverage mean (%f) or variance (%f)\n\n",
                  allelicCovMean, allelicCovVar);
          PrintUsage(stderr);
        }

        if (allelicCovMean > 0)
          doNGS = YES;
      }

      if (genotypingError > 0 && doNGS == YES) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot specify a coverage "
                        "larger than 0, which implies read count generation, "
                        "and a genotyping error at the same time\n\n");
        PrintUsage(stderr);
      }
      /* end added code */

      /* original code */
      /*if (fscanf(stdin, "%f", &argument) != 1 || argument < 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad sequencing coverage (%d)\n\n",
                (int)coverage);
        PrintUsage(stderr);
      }
      coverage = (int)argument;
      if (coverage > 0)
        doNGS = YES;
      if (genotypingError > 0 && doNGS == YES) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot specify a coverage "
                        "larger than 0, which implies read count generation, "
                        "and a genotyping error at the same time\n\n");
        PrintUsage(stderr);
      }*/
      break;
    case 'V':
      /* added code */
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        if (fscanf(stdin, "%lf", &alphaCoverage) != 1 || alphaCoverage <= 0) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: Bad coverage dispersion (%f)\n\n",
                  alphaCoverage);
          PrintUsage(stderr);
        }
        rateVarCoverage = YES;
      } else {
        if (fscanf(stdin, "%lf %lf", &allelicRawVarMean, &allelicRawVarVar) !=
                2 ||
            allelicRawVarMean <= 0 || allelicRawVarVar < 0) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: Bad allelic sequencing raw "
                  "variance mean (%f) or variance (%f)\n\n",
                  allelicRawVarMean, allelicRawVarVar);
          PrintUsage(stderr);
        }
      }
      /* end added code */

      /* original code
      if (fscanf(stdin, "%lf", &alphaCoverage) != 1 || alphaCoverage <= 0) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad coverage dispersion (%f)\n\n",
                alphaCoverage);
        PrintUsage(stderr);
      }
      rateVarCoverage = YES;
       */
      break;
      /* added code */
    case 'K':
      if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu) {
        if (fscanf(stdin, "%lf %lf", &meanSizeFactor, &varSizeFactor) != 2) {
          fprintf(stderr, "\n!!! PARAMETER ERROR: Expecting two parameters for "
                          "size factors.\n\n");
          PrintUsage(stderr);
        }
        if (meanSizeFactor < 0) {
          fprintf(stderr,
                  "\n!!! PARAMETER ERROR: Bad mean value for size factors "
                  "(%lf).\n\n",
                  meanSizeFactor);
          PrintUsage(stderr);
        }
        if (varSizeFactor < 0) {
          fprintf(
              stderr,
              "\n!!! PARAMETER ERROR: Bad variance for size factors (%lf).\n\n",
              varSizeFactor);
          PrintUsage(stderr);
        }
      } else {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Size factors should be defined under "
                "the finite mutation insertion deletion model (m3)!\n\n");
        PrintUsage(stderr);
      }
      break;
    case 'N':
      if (altModel == FiniteMuDel || altModel == FiniteMu || altModel == FiniteMuInDel) {
        if (fscanf(stdin, "%lf %lf", &shapeController[0],
                   &shapeController[1]) != 2) {
          fprintf(stderr, "\n!!! PARAMETER ERROR: Expecting two parameters for "
                          "shape controllers.\n\n");
          PrintUsage(stderr);
        }
      }
//      else if (altModel == FiniteMuInDel) {
//        fprintf(stderr, "\n!!! PARAMETER ERROR: Undefined behavior.\n\n");
//        PrintUsage(stderr);
//      }
      else {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Shape controllers should work with "
                "m3, m4, or m5.\n\n");
        PrintUsage(stderr);
      }
      break;
      //  end added code
    case 'A':
      /* added code */
      if (!(altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu)) {
        if (fscanf(stdin, "%lf %lf %d", &meanAmplificationError,
                   &varAmplificationError, &simulateOnlyTwoTemplates) != 3) {
          fprintf(
              stderr,
              "\n!!! PARAMETER ERROR: Bad mean/var/model amplification error "
              "(%f ; %f ; model=%d)\n\n",
              meanAmplificationError, varAmplificationError,
              simulateOnlyTwoTemplates);
          PrintUsage(stderr);
        }
        if (meanAmplificationError < 0 || meanAmplificationError > 1) {
          fprintf(
              stderr,
              "\n!!! PARAMETER ERROR: Bad mean amplification error (%f)\n\n",
              meanAmplificationError);
          PrintUsage(stderr);
        }
        if (varAmplificationError <= 0 ||
            (meanAmplificationError > 0 &&
             varAmplificationError >=
                 (meanAmplificationError * (1.0 - meanAmplificationError)))) {
          fprintf(
              stderr,
              "\n!!! PARAMETER ERROR: Bad variance amplification error (%f); "
              "it has to be >0 and < mean*(1-mean)\n\n",
              varAmplificationError);
          PrintUsage(stderr);
        }
        if (simulateOnlyTwoTemplates != 0 && simulateOnlyTwoTemplates != 1) {
          fprintf(
              stderr,
              "\n!!! PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); "
              "it has to be 0 (assume 4 templates) or 1 (assume 2 templates)",
              simulateOnlyTwoTemplates);
          PrintUsage(stderr);
        }
      } else {
        if (fscanf(stdin, "%lf %lf %d", &meanAmplificationError,
                   &varAmplificationError, &fixedAmplificationError) != 3) {
          fprintf(
              stderr,
              "\n!!! PARAMETER ERROR: Bad mean/var/fixed amplification error "
              "(%f ; %f ; model=%d)\n\n",
              meanAmplificationError, varAmplificationError,
              fixedAmplificationError);
          PrintUsage(stderr);
        }
        if (meanAmplificationError < 0 || meanAmplificationError > 1) {
          fprintf(
              stderr,
              "\n!!! PARAMETER ERROR: Bad mean amplification error (%f)\n\n",
              meanAmplificationError);
          PrintUsage(stderr);
        }
        if (varAmplificationError <= 0 ||
            (meanAmplificationError > 0 &&
             varAmplificationError >=
                 (meanAmplificationError * (1.0 - meanAmplificationError)))) {
          fprintf(
              stderr,
              "\n!!! PARAMETER ERROR: Bad variance amplification error (%f); "
              "it has to be >0 and < mean*(1-mean)\n\n",
              varAmplificationError);
          PrintUsage(stderr);
        }
        if (fixedAmplificationError != 0 && fixedAmplificationError != 1) {
          fprintf(
              stderr,
              "\n!!! PARAMETER ERROR: Bad fixedAmplificationError error (%d); "
              "it has to be 0 (assume site-wise amplification error) or 1 "
              "(assume fixed amplification error); using 4 templates.",
              fixedAmplificationError);
          PrintUsage(stderr);
        }
      }
      /* end added code */

      /* original code */
      /*
      if (fscanf(stdin, "%lf %lf %d", &meanAmplificationError,
                 &varAmplificationError, &simulateOnlyTwoTemplates) != 3) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad mean/var/model amplification error "
                "(%f ; %f ; model=%d)\n\n",
                meanAmplificationError, varAmplificationError,
                simulateOnlyTwoTemplates);
        PrintUsage(stderr);
      }
      if (meanAmplificationError < 0 || meanAmplificationError > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad mean amplification error (%f)\n\n",
                meanAmplificationError);
        PrintUsage(stderr);
      }
      if (varAmplificationError <= 0 ||
          (meanAmplificationError > 0 &&
           varAmplificationError >=
               (meanAmplificationError * (1.0 - meanAmplificationError)))) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad variance amplification error (%f); "
                "it has to be >0 and < mean*(1-mean)\n\n",
                varAmplificationError);
        PrintUsage(stderr);
      }
      if (simulateOnlyTwoTemplates != 0 && simulateOnlyTwoTemplates != 1) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); "
            "it has to be 0 (assume 4 templates) or 1 (assume 2 templates)",
            simulateOnlyTwoTemplates);
        PrintUsage(stderr);
      }
       */
      break;
    case 'E':
      if (fscanf(stdin, "%lf", &sequencingError) != 1 || sequencingError < 0 ||
          sequencingError > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad sequencing error (%f)\n\n",
                sequencingError);
        PrintUsage(stderr);
      }
      break;
    case 'D':
      if (fscanf(stdin, "%lf", &ADOrate) != 1 || ADOrate < 0 || ADOrate > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad allelic dropout rate (%f)\n\n",
                ADOrate);
        PrintUsage(stderr);
      }
      break;

    /* added code */
    case 'Y':
      if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
          altModel == FiniteMu) {
        doPrintHealthyTip = YES;
      } else {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: print healthy tip flag should be defined "
            "under the finite mutation insertion deletion model "
            "(m3, m4, or m5)!\n\n");
        PrintUsage(stderr);
      }
      break;
    case 'Z':
      if (fscanf(stdin, "%d", &ADOtype) == 1) {
        if (altModel == FiniteMuInDel || altModel == FiniteMuDel ||
            altModel == FiniteMu) {
          if (ADOtype < 0 || ADOtype > 2) {
            fprintf(stderr,
                    "\n!!! PARAMETER ERROR: Bad allelic dropout type (%d)\n\n",
                    ADOtype);
            PrintUsage(stderr);
          }
        } else {
          fprintf(stderr, "\n!!! PARAMETER ERROR: Invaid definition of Z.\n\n");
          PrintUsage(stderr);
        }
      }
      break;
    case 'P':
      if (fscanf(stdin, "%lf", &alphaADOsites) != 1 || alphaADOsites < 0 ||
          alphaADOsites > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad allelic dropout hetereogeneity "
                "among sites (%f)\n\n",
                alphaADOsites);
        PrintUsage(stderr);
      }
      ADOvarAmongSites = YES;
      break;
    case 'Q':
      if (fscanf(stdin, "%lf", &alphaADOcells) != 1 || alphaADOcells < 0 ||
          alphaADOcells > 1) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad allelic dropout hetereogeneity "
                "among cells (%f)\n\n",
                alphaADOcells);
        PrintUsage(stderr);
      }
      ADOvarAmongCells = YES;
      break;
    case 'I':
      if (fscanf(stdin, "%lf", &allelicImbalance) != 1 ||
          allelicImbalance < 0 || allelicImbalance > 1) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad allelic imbalance (%f)\n\n",
                allelicImbalance);
        PrintUsage(stderr);
      }
      break;
    case 'R':
      if (fscanf(stdin, "%lf", &haploidCoverageReduction) != 1 ||
          haploidCoverageReduction < 0 || haploidCoverageReduction > 1) {
        fprintf(
            stderr,
            "\n!!! PARAMETER ERROR: Bad haploid coverage reduction (%f)\n\n",
            haploidCoverageReduction);
        PrintUsage(stderr);
      }
      break;
    case 'B':
      if (fscanf(stdin, "%lf %lf", &meanDoubletRate, &varDoubletRate) != 2)
      {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mean/var/ doublet rate (%f ; %f)\n\n", meanDoubletRate, varDoubletRate);
        PrintUsage(stderr);
      }
      if (meanDoubletRate < 0 || meanDoubletRate > 1)
      {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mean doublet rate (%f)\n\n", meanDoubletRate);
        PrintUsage(stderr);
      }
      if (varDoubletRate <= 0 || (meanDoubletRate > 0 && varDoubletRate >= (meanDoubletRate * (1.0 - meanDoubletRate))))
      {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad variance doublet rate error (%f); it has to be >0 and < mean*(1-mean)\n\n", varDoubletRate);
        PrintUsage(stderr);
      }
      break;
    case 'X':
      if (fscanf(stdin,
                 "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf "
                 "%lf",
                 &Eij[0][0], &Eij[0][1], &Eij[0][2], &Eij[0][3], &Eij[1][0],
                 &Eij[1][1], &Eij[1][2], &Eij[1][3], &Eij[2][0], &Eij[2][1],
                 &Eij[2][2], &Eij[2][3], &Eij[3][0], &Eij[3][1], &Eij[3][2],
                 &Eij[3][3]) != 16) {
        fprintf(stderr,
                "\n!!! PARAMETER ERROR: Bad  error matrix (-rx x x x x x x x x "
                "x x) (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT)\n\n");
        PrintUsage(stderr);
      }

      if (Eij[0][0] != 0 || Eij[1][1] != 0 || Eij[2][2] != 0 ||
          Eij[3][3] != 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad error matrix: diagonals "
                        "should be 0 \n\n");
        PrintUsage(stderr);
      }
      thereIsEij = YES;
      break;
      /* added code */
    case '0':
      doPrintFullNucleotideReads = YES;
      break;
      /* end added code */
    case '1':
      doPrintSNVgenotypes = YES;
      break;
    case '2':
      doPrintSNVhaplotypes = YES;
      break;
    case '3':
      doPrintFullGenotypes = YES;
      /* added code */
      if (numSites >= 1000) {
        fprintf(stderr, "\nWARNING! INVALID PARAMETER: number of sites is more "
                        "than 1000. Full genotype will not be printed in order "
                        "to save disk space. \n\n");
        doPrintFullGenotypes = NO;
      }
      /* end added code */
      break;
    case '4':
      doPrintFullHaplotypes = YES;
      break;
    case '5':
      doPrintAncestors = YES;
      break;
    case '6':
      doPrintTree = YES;
      break;
    case '7':
      doPrintTimes = YES;
      break;
    case '8':
      doPrintCATG = YES;
      if (doNGS == NO) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Cannot print the CATG format "
                        "when reads are not simulated (coverage >0)\n\n");
        PrintUsage(stderr);
      }
      break;
    case '9':
      doPrintTrueHaplotypes = YES;
      break;
    case 'v':
      doPrintSeparateReplicates = YES;
      break;
    case 'x':
      doPrintIUPAChaplotypes = YES;
      break;
    case 'o':
      ch = fgetc(stdin);
      if (!isspace(ch)) {
        j = 0;
        do {
          resultsDir[j] = ch;
          j++;
          ch = fgetc(stdin);
        } while (!isspace(ch));
        resultsDir[j] = '\0';
      }
      break;
    case 'T':
      ch = fgetc(stdin);
      if (!isspace(ch)) {
        j = 0;
        do {
          userTreeFile[j] = ch;
          j++;
          ch = fgetc(stdin);
        } while (!isspace(ch));
        userTreeFile[j] = '\0';
      }
      doUserTree = YES;
      break;
    case 'U':
      ch = fgetc(stdin);
      if (!isspace(ch)) {
        j = 0;
        do {
          userGenomeFile[j] = ch;
          j++;
          ch = fgetc(stdin);
        } while (!isspace(ch));
        userGenomeFile[j] = '\0';
      }
      doUserGenome = YES;
      break;
    /*case '0':
doSimulateData = NO;
    break;*/
    case 'W':
      doTumorNames = YES;
      break;
    case 'y':
      /* original code */
      // if (fscanf(stdin, "%f", &argument) != 1 || argument < 0 || argument >
      // 3)

      /* added code */
      if (fscanf(stdin, "%f", &argument) != 1 || argument < 0 || argument > 4)
      /* end added code */
      {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad noisy value (%d)\n\n",
                (int)argument);
        PrintUsage(stderr);
      }
      noisy = (int)argument;
      break;
    case '*':
      if (fscanf(stdin, "%d", readLengths + 1) != 1 || readLengths[1] <= 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad read length (%d)\n\n",
                readLengths[1]);
        PrintUsage(stderr);
      }

      SetReadLengths(numSites);

      break;
    case '#':
      if (fscanf(stdin, "%lf", &argdouble) != 1 || argument < 0) {
        fprintf(stderr, "\n!!! PARAMETER ERROR: Bad seed (%ld)\n\n",
                (long int)argdouble);
        PrintUsage(stderr);
      }
      userSeed = (unsigned long)argdouble;
      break;
    case '?':
      PrintUsage(stderr);
      break;
    default:
      fprintf(stderr, "\n!!! PARAMETER ERROR: Incorrect parameter: %c\n\n", ch);
      PrintUsage(stderr);
      break;
    }
    ch = fgetc(stdin);
    while (isspace(ch) && !feof(stdin))
      ch = fgetc(stdin);
    while (ch == '[') {
      ReadUntil(stdin, ']', "closing bracket");
      ch = fgetc(stdin);
      while (isspace(ch))
        ch = fgetc(stdin);
    }
  }
}
