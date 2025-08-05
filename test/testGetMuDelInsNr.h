//
// Created by Senbai Kang on 11/25/19.
//

#ifndef MYSIMULATOR_TESTGETMUDELINSNR_H
#define MYSIMULATOR_TESTGETMUDELINSNR_H

struct HashMapIntToIntArr *map = NULL;
DataStr *dataPa = NULL, *dataCh = NULL;
unsigned int evolutionaryEvents[5];

void AllocateData() {
  dataPa = (DataStr *)calloc(2, sizeof(DataStr));
  dataPa->alleles = NULL;
  dataPa->pastEvolutionaryEvents = NULL;
  dataCh = (DataStr *)calloc(2, sizeof(DataStr));
  dataCh->alleles = NULL;
  dataCh->pastEvolutionaryEvents = NULL;
}

void FreeData() {
  free(dataCh);
  dataCh = NULL;
  free(dataPa);
  dataPa = NULL;
}

void GetMuDelInsNrCore(int referenceNuc, int maxCopies) {
  int overlapNr, overlap[maxCopies];

  {
    for (int i = 0; i < 5; i++)
      evolutionaryEvents[i] = 0;
  }

  if (dataCh->numExistingAlleles != dataPa->numExistingAlleles) {
    if (dataPa->numExistingAlleles == 3) {
      evolutionaryEvents[3]++;

      if (dataCh->numExistingAlleles == 2) {
        if (dataCh->alleles[0] < 0) {
          // e.g., A/AA -> -/AA
          overlapNr = NumOverlap(dataPa->alleles, dataCh->alleles, maxCopies, 1,
                                 1, overlap);

          CountMutationNum(dataCh->alleles, 1, maxCopies, referenceNuc,
                           dataPa->alleles, overlap, evolutionaryEvents);
        } else {
          overlapNr = NumOverlap(dataPa->alleles, dataCh->alleles, maxCopies, 0,
                                 0, overlap);
          if (overlapNr == 2) {
            if (dataPa->alleles[0] != dataCh->alleles[0] &&
                dataPa->alleles[0] != dataCh->alleles[1]) {
              // e.g., A/CG -> C/G
              // Ambiguity will raise.
              // allSites[site].numMutations++;
              // numMU++;
              evolutionaryEvents[0]++;
              if (ValueInArray(referenceNuc, dataCh->alleles, 2) == 1)
                evolutionaryEvents[2]++;
            }
          } else {
            CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                             dataPa->alleles, overlap, evolutionaryEvents);
          }
        }
      } else if (dataCh->numExistingAlleles == 1) {
        evolutionaryEvents[3]++;

        overlapNr = NumOverlap(dataPa->alleles, dataCh->alleles, maxCopies, 0,
                               0, overlap);
        if (overlapNr == 0) {
          CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                           dataPa->alleles, overlap, evolutionaryEvents);
        }
      } else {
        evolutionaryEvents[3] += 2;
      }
    } else if (dataPa->numExistingAlleles == 2) {
      if (dataCh->numExistingAlleles == 3) {
        if (dataPa->alleles[0] >= 0) {
          evolutionaryEvents[4]++;

          overlapNr = NumOverlap(dataPa->alleles, dataCh->alleles, maxCopies, 0,
                                 0, overlap);
          if (dataPa->alleles[0] == dataPa->alleles[1]) {
            // e.g., G/G
            if (overlapNr == 2) {
              if (dataCh->alleles[0] != dataCh->alleles[1] ||
                  dataCh->alleles[0] != dataCh->alleles[2] ||
                  dataCh->alleles[1] != dataCh->alleles[2]) {
                CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                                 dataPa->alleles, overlap, evolutionaryEvents);
              }
            } else if (overlapNr == 1) {
              if (dataCh->alleles[1] == dataCh->alleles[2]) {
                // e.g., G/G -> G/AA
                CountMutationNum(dataCh->alleles, 0, 2, referenceNuc,
                                 dataPa->alleles, overlap, evolutionaryEvents);
              } else {
                // e.g., G/G -> A/CG
                CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                                 dataPa->alleles, overlap, evolutionaryEvents);
              }
            } else {
              if (dataCh->alleles[1] == dataCh->alleles[2]) {
                // e.g., G/G -> A/CC
                CountMutationNum(dataCh->alleles, 0, 2, referenceNuc,
                                 dataPa->alleles, overlap, evolutionaryEvents);
              } else {
                // e.g., G/G -> A/CT
                CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                                 dataPa->alleles, overlap, evolutionaryEvents);
              }
            }
          } else {
            // e.g., C/T
            if (overlapNr == 2) {
              if (dataCh->alleles[1] != dataCh->alleles[2]) {
                if ((dataCh->alleles[0] != dataPa->alleles[0]) &&
                    (dataCh->alleles[0] != dataPa->alleles[1])) {
                  // e.g., C/T -> A/CT
                  // Ambiguity could raise.
                  evolutionaryEvents[0] += 2;
                  if (referenceNuc == dataCh->alleles[0])
                    evolutionaryEvents[1] += 1;
                  if (referenceNuc == dataCh->alleles[1] ||
                      referenceNuc == dataCh->alleles[2])
                    evolutionaryEvents[2] += 1;
                } else {
                  // e.g., C/T -> C/AT
                  CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                                   dataPa->alleles, overlap,
                                   evolutionaryEvents);
                }
              }
            } else if (overlapNr == 1) {
              if (dataCh->alleles[1] != dataCh->alleles[2]) {
                // e.g., C/T -> C/AG
                CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                                 dataPa->alleles, overlap, evolutionaryEvents);
              } else {
                // e.g., C/T -> A/TT
                CountMutationNum(dataCh->alleles, 0, 2, referenceNuc,
                                 dataPa->alleles, overlap, evolutionaryEvents);
              }
            } else {
              if (dataCh->alleles[1] != dataCh->alleles[2]) {
                // e.g., C/T -> G/AG
                CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                                 dataPa->alleles, overlap, evolutionaryEvents);
              } else {
                // e.g., C/T -> G/AA
                CountMutationNum(dataCh->alleles, 0, 2, referenceNuc,
                                 dataPa->alleles, overlap, evolutionaryEvents);
              }
            }
          }
        }
      } else if (dataCh->numExistingAlleles == 1) {
        evolutionaryEvents[3]++;

        overlapNr = NumOverlap(dataPa->alleles, dataCh->alleles, maxCopies, 0,
                               0, overlap);
        if (overlapNr == 0) {
          // e.g., C/T -> A/-
          CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                           dataPa->alleles, overlap, evolutionaryEvents);
        }
      } else {
        evolutionaryEvents[3] += 2;
      }
    } else if (dataPa->numExistingAlleles == 1) {
      // e.g., C/-
      if (dataCh->numExistingAlleles == 2 && dataCh->alleles[0] == DELETED) {
        evolutionaryEvents[4]++;

        overlapNr = NumOverlap(dataPa->alleles, dataCh->alleles, maxCopies, 0,
                               0, overlap);
        if (overlapNr == 0) {
          if (dataCh->alleles[1] == dataCh->alleles[2]) {
            // e.g., C/- -> -/AA
            CountMutationNum(dataCh->alleles, 1, 2, referenceNuc,
                             dataPa->alleles, overlap, evolutionaryEvents);
          } else {
            // e.g., C/- -> -/GT
            CountMutationNum(dataCh->alleles, 1, maxCopies, referenceNuc,
                             dataPa->alleles, overlap, evolutionaryEvents);
          }
        } else {
          // e.g., C/- -> -/CT
          if (dataCh->alleles[1] != dataCh->alleles[2]) {
            CountMutationNum(dataCh->alleles, 1, maxCopies, referenceNuc,
                             dataPa->alleles, overlap, evolutionaryEvents);
          }
        }
      } else if (dataCh->numExistingAlleles == 0) {
        evolutionaryEvents[3]++;
      }
    }
  } else {
    if (dataPa->numExistingAlleles == 3) {
      if (dataPa->alleles[0] != dataCh->alleles[0]) {
        // e.g., A/CC -> G/??
        evolutionaryEvents[0]++;
        if (referenceNuc == dataCh->alleles[0])
          evolutionaryEvents[1]++;
      }

      overlapNr = NumOverlap(dataPa->alleles, dataCh->alleles, maxCopies, 1, 1,
                             overlap);
      CountMutationNum(dataCh->alleles, 1, maxCopies, referenceNuc,
                       dataPa->alleles, overlap, evolutionaryEvents);
    } else if (dataPa->numExistingAlleles == 1) {
      CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                       dataPa->alleles, overlap, evolutionaryEvents);
    } else {
      if (dataPa->alleles[0] != DELETED && dataCh->alleles[0] == DELETED) {
        evolutionaryEvents[3]++;
        evolutionaryEvents[4]++;

        overlapNr = NumOverlap(dataPa->alleles, dataCh->alleles, maxCopies, 0,
                               1, overlap);
        if (overlapNr == 2) {
          if (dataPa->alleles[0] != dataPa->alleles[1]) {
            // e.g., G/T -> -/GT
            // Ambiguity might raise.
            evolutionaryEvents[0]++;
            if (ValueInArray(referenceNuc, dataCh->alleles, 3) == 1)
              evolutionaryEvents[2]++;
          }
        } else if (overlapNr == 1) {
          if (dataCh->alleles[1] != dataCh->alleles[2]) {
            // e.g., G/T -> -/AT
            CountMutationNum(dataCh->alleles, 1, maxCopies, referenceNuc,
                             dataPa->alleles, overlap, evolutionaryEvents);
          }
        } else {
          if (dataCh->alleles[1] == dataCh->alleles[2]) {
            // e.g., G/T -> -/AA
            CountMutationNum(dataCh->alleles, 1, 2, referenceNuc,
                             dataPa->alleles, overlap, evolutionaryEvents);
          } else {
            // e.g., G/T -> -/AC
            CountMutationNum(dataCh->alleles, 1, maxCopies, referenceNuc,
                             dataPa->alleles, overlap, evolutionaryEvents);
          }
        }
      } else if (!(dataPa->alleles[0] == DELETED &&
                   dataCh->alleles[0] != DELETED)) {
        overlapNr = NumOverlap(dataPa->alleles, dataCh->alleles, maxCopies, 0,
                               0, overlap);
        CountMutationNum(dataCh->alleles, 0, maxCopies, referenceNuc,
                         dataPa->alleles, overlap, evolutionaryEvents);
      }
    }
  }
}

int GetCodeMD(int code) {
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

int GetCodeMDI(int code) {
  if (code < 40) {
    return 3;
  } else if (code < 60 && code >= 40) {
    return 2;
  } else if (code < 64 && code >= 60) {
    return 1;
  } else if (code == 64) {
    return 0;
  }
  return -1;
}

void SetDataStrMD(DataStr *dataStr, int code) {
  dataStr->code = code;
  dataStr->numExistingAlleles = GetCodeMD(code);
  dataStr->alleles = GetValIntToIntArr(map, dataStr->code);
}

void SetDataStrMDI(DataStr *dataStr, int code) {
  dataStr->code = code;
  dataStr->numExistingAlleles = GetCodeMDI(code);
  dataStr->alleles = GetValIntToIntArr(map, dataStr->code);
}

#endif // MYSIMULATOR_TESTGETMUDELINSNR_H