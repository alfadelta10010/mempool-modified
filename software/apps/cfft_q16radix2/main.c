// Copyright 2022 ETH Zurich and University of Bologna.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

// Author: Marco Bertuletti, ETH Zurich

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Mempool runtime libraries */
#include "encoding.h"
#include "printf.h"
#include "runtime.h"
#include "synchronization.h"
#include "xpulp/builtins_v2.h"

/* CFFT mempool libraries */
#include "kernel/mempool_radix2_cfft_q16p.h"
#include "kernel/mempool_radix2_cfft_q16s.h"
#include "tables/mempool_radix2_cfft_q16_BitRevIndexTable.h"
#include "tables/mempool_radix2_cfft_q16_twiddleCoef.h"

#define N_CSAMPLES 256
#define N_RSAMPLES 2 * N_CSAMPLES
//#define VERBOSE
#define TEST_256
#define PARALLEL
#define SINGLE

/* CFFT mempool data */
int16_t pSrc[N_RSAMPLES] __attribute__((aligned(N_CSAMPLES), section(".l1")));

void initialize_vector(int16_t *pSrc, uint32_t N_el) {
  int lower = SHRT_MIN, upper = SHRT_MAX;
  srand((unsigned)1);
  for (uint32_t i = 0; i < N_el; i++) {
    pSrc[i] = (int16_t)((rand() % (upper - lower + 1)) + lower);
  }
}

int main() {

  uint32_t core_id = mempool_get_core_id();
  uint32_t num_cores = mempool_get_core_count();
  mempool_barrier_init(core_id);

  /* SINGLE-CORE */

#ifdef SINGLE

  if (core_id == 0) {

    printf("On the run...\n");
    initialize_vector(pSrc, N_RSAMPLES);

#ifdef TEST_16
    mempool_start_benchmark();
    mempool_radix2_cfft_q16s((uint16_t)16, twiddleCoef_16_q16,
                             BitRevIndexTable_fixed_16, pSrc,
                             BITREVINDEXTABLE_FIXED_16_TABLE_LENGTH, 0, 0);
    mempool_stop_benchmark();

#elif defined(TEST_32)
    mempool_start_benchmark();
    mempool_radix2_cfft_q16s((uint16_t)32, twiddleCoef_32_q16,
                             BitRevIndexTable_fixed_32, pSrc,
                             BITREVINDEXTABLE_FIXED_32_TABLE_LENGTH, 0, 0);
    mempool_stop_benchmark();
#elif defined(TEST_64)
    mempool_start_benchmark();
    mempool_radix2_cfft_q16s((uint16_t)64, twiddleCoef_64_q16,
                             BitRevIndexTable_fixed_64, pSrc,
                             BITREVINDEXTABLE_FIXED_64_TABLE_LENGTH, 0, 0);
    mempool_stop_benchmark();
#elif defined(TEST_128)
    mempool_start_benchmark();
    mempool_radix2_cfft_q16s((uint16_t)128, twiddleCoef_128_q16,
                             BitRevIndexTable_fixed_128, pSrc,
                             BITREVINDEXTABLE_FIXED_128_TABLE_LENGTH, 0, 0);
    mempool_stop_benchmark();
#elif defined(TEST_256)
    mempool_start_benchmark();
    mempool_radix2_cfft_q16s((uint16_t)256, twiddleCoef_256_q16,
                             BitRevIndexTable_fixed_256, pSrc,
                             BITREVINDEXTABLE_FIXED_256_TABLE_LENGTH, 0, 0);
    mempool_stop_benchmark();
#elif defined(TEST_512)
    mempool_start_benchmark();
    mempool_radix2_cfft_q16s((uint16_t)512, twiddleCoef_512_q16,
                             BitRevIndexTable_fixed_512, pSrc,
                             BITREVINDEXTABLE_FIXED_512_TABLE_LENGTH, 0, 0);
    mempool_stop_benchmark();
#elif defined(TEST_1024)
    mempool_start_benchmark();
    mempool_radix2_cfft_q16s((uint16_t)1024, twiddleCoef_1024_q16,
                             BitRevIndexTable_fixed_1024, pSrc,
                             BITREVINDEXTABLE_FIXED_1024_TABLE_LENGTH, 0, 0);
    mempool_stop_benchmark();
#elif defined(TEST_2048)
    mempool_start_benchmark();
    mempool_radix2_cfft_q16s((uint16_t)2048, twiddleCoef_2048_q16,
                             BitRevIndexTable_fixed_2048, pSrc,
                             BITREVINDEXTABLE_FIXED_2048_TABLE_LENGTH, 0, 0);
    mempool_stop_benchmark();
#elif defined(TEST_4096)
    mempool_start_benchmark();
    mempool_radix2_cfft_q16s((uint16_t)4096, twiddleCoef_4096_q16,
                             BitRevIndexTable_fixed_4096, pSrc,
                             BITREVINDEXTABLE_FIXED_4096_TABLE_LENGTH, 0, 0);
    mempool_stop_benchmark();
#endif

#ifdef VERBOSE
    for (uint32_t i = 0; i < N_RSAMPLES; i += 2) {
      printf("{%6d;%6d } \n", pSrc[i], pSrc[i + 1]);
    }
#endif
    printf("Done SINGLE!\n");
  }
  mempool_barrier(num_cores);

#endif

  /* PARALLEL-CORE */

#ifdef PARALLEL

  if (core_id == 0) {
    printf("On the run...\n");
    initialize_vector(pSrc, N_RSAMPLES);
  }
  mempool_barrier(num_cores);

#ifdef TEST_16
  mempool_start_benchmark();
  mempool_radix2_cfft_q16p(
      (uint16_t)16, twiddleCoef_16_q16, BitRevIndexTable_fixed_16, pSrc,
      BITREVINDEXTABLE_FIXED_16_TABLE_LENGTH, 0, 0, num_cores);
  mempool_stop_benchmark();

#elif defined(TEST_32)
  mempool_start_benchmark();
  mempool_radix2_cfft_q16p(
      (uint16_t)32, twiddleCoef_32_q16, BitRevIndexTable_fixed_32, pSrc,
      BITREVINDEXTABLE_FIXED_32_TABLE_LENGTH, 0, 0, num_cores);
  mempool_stop_benchmark();
#elif defined(TEST_64)
  mempool_start_benchmark();
  mempool_radix2_cfft_q16p(
      (uint16_t)64, twiddleCoef_64_q16, BitRevIndexTable_fixed_64, pSrc,
      BITREVINDEXTABLE_FIXED_64_TABLE_LENGTH, 0, 0, num_cores);
  mempool_stop_benchmark();
#elif defined(TEST_128)
  mempool_start_benchmark();
  mempool_radix2_cfft_q16p(
      (uint16_t)128, twiddleCoef_128_q16, BitRevIndexTable_fixed_128, pSrc,
      BITREVINDEXTABLE_FIXED_128_TABLE_LENGTH, 0, 0, num_cores);
  mempool_stop_benchmark();
#elif defined(TEST_256)
  mempool_start_benchmark();
  mempool_radix2_cfft_q16p(
      (uint16_t)256, twiddleCoef_256_q16, BitRevIndexTable_fixed_256, pSrc,
      BITREVINDEXTABLE_FIXED_256_TABLE_LENGTH, 0, 0, num_cores);
  mempool_stop_benchmark();
#elif defined(TEST_512)
  mempool_start_benchmark();
  mempool_radix2_cfft_q16p(
      (uint16_t)512, twiddleCoef_512_q16, BitRevIndexTable_fixed_512, pSrc,
      BITREVINDEXTABLE_FIXED_512_TABLE_LENGTH, 0, 0, num_cores);
  mempool_stop_benchmark();
#elif defined(TEST_1024)
  mempool_start_benchmark();
  mempool_radix2_cfft_q16p(
      (uint16_t)1024, twiddleCoef_1024_q16, BitRevIndexTable_fixed_1024, pSrc,
      BITREVINDEXTABLE_FIXED_1024_TABLE_LENGTH, 0, 0, num_cores);
  mempool_stop_benchmark();
#elif defined(TEST_2048)
  mempool_start_benchmark();
  mempool_radix2_cfft_q16p(
      (uint16_t)2048, twiddleCoef_2048_q16, BitRevIndexTable_fixed_2048, pSrc,
      BITREVINDEXTABLE_FIXED_2048_TABLE_LENGTH, 0, 0, num_cores);
  mempool_stop_benchmark();
#elif defined(TEST_4096)
  mempool_start_benchmark();
  mempool_radix2_cfft_q16p(
      (uint16_t)4096, twiddleCoef_4096_q16, BitRevIndexTable_fixed_4096, pSrc,
      BITREVINDEXTABLE_FIXED_4096_TABLE_LENGTH, 0, 0, num_cores);
  mempool_stop_benchmark();
#endif

  if (core_id == 0) {
#ifdef VERBOSE
    for (uint32_t i = 0; i < N_RSAMPLES; i += 2) {
      printf("{%6d;%6d } \n", pSrc[i], pSrc[i + 1]);
    }
#endif
    printf("Done PARALLEL!\n");
  }
  mempool_barrier(num_cores);

#endif

  return 0;
}
