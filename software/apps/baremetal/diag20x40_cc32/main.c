#include <stdint.h>
#include "encoding.h"
#include <stdio.h>
#include <stdbool.h>
#include "printf.h"
#include "runtime.h"
#include "synchronization.h"
#include <string.h>

// DEFINITIONS

#define MAX_SEQ_LENGTH 200
#define MATCH_SCORE 2
#define MISMATCH_SCORE 1
#define GAP_SCORE 1

// CORE_COUNT
#define cc 32

// TIMERS
int global_compute_time = 0;
int global_sync_time = 0;

// Cell update counter
int global_cell_updates = 0;

// compute does not include indices calculation duration
uint32_t compute_start_times[cc] = {0};
uint32_t compute_end_times[cc] = {0};

// SEQUENCE DATA
const char seq1[] = "GATGGGCGTGCCAAAGGAGATCGGACGAAATGGATTTCGA";
const char seq2[] = "TGCGGGCCAAGATAGGCGGT";
int seq1_len = strlen(seq1);
int seq2_len = strlen(seq2);

// DP matrix of size m*n not (m+1)(n+1)
int16_t dp[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH] __attribute__((section(".l1")));

// Function to compute score
int16_t score(char a, char b) {
  if (a == b) return MATCH_SCORE;
  return MISMATCH_SCORE;
}

// function to find min
int16_t minimum(int16_t x, int16_t y) {
  return (x < y) ? x : y;
}

// function to find max
int16_t maximum(int16_t x, int16_t y) {
  return (x > y) ? x : y;
}

// Unified max calculation function
int16_t max3(int16_t a, int16_t b, int16_t c) {
  return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

// computes length of diagonal given index
int diagonal_length_cal(int diagonal_index) {
  int num_diagonal_elements;
  int min_dim = minimum(seq1_len, seq2_len);
  int max_dim = maximum(seq1_len, seq2_len);
  if (diagonal_index < min_dim) {
    num_diagonal_elements = diagonal_index + 1;
  } else if (diagonal_index < max_dim) {
    num_diagonal_elements = min_dim;
  } else {
    num_diagonal_elements = seq1_len + seq2_len - 1 - diagonal_index;
  }
  return num_diagonal_elements;
}

// Function to compute a segment of a diagonal
void compute_diagonal_segment(int diagonal_index, int start, int end) {
    // Directly compute elements within the specified segment of the diagonal
    int count = 0;

    // Calculate the exact starting (i, j) coordinates for the first element
    int initial_i, initial_j;
    if (diagonal_index < seq2_len) {
        initial_i = diagonal_index;
        initial_j = 0;
    } else {
        initial_i = seq1_len - 1;
        initial_j = diagonal_index - seq1_len + 1;
    }

    // Move to the correct starting point in the segment
    while (count < start) {
        initial_i--;
        initial_j++;
        count++;
    }

    // Compute only the required segment
    while (count < end) {
        // Ensure we're within sequence bounds
        if (initial_i >= 0 && initial_j < seq2_len) {
            int16_t match, delete, insert, max_score;

            // Simplified condition handling (similar to your existing logic)
            if (initial_i == 0 && initial_j == 0) {
                match = score(seq1[initial_i], seq2[initial_j]);
                max_score = max3(match, GAP_SCORE, GAP_SCORE);
            } else if (initial_i == 0) {
                match = score(seq1[initial_i], seq2[initial_j]);
                max_score = max3(match, GAP_SCORE, dp[initial_i][initial_j-1] + GAP_SCORE);
            } else if (initial_j == 0) {
                match = score(seq1[initial_i], seq2[initial_j]);
                max_score = max3(match, dp[initial_i-1][initial_j] + GAP_SCORE, GAP_SCORE);
            } else {
                match = dp[initial_i-1][initial_j-1] + score(seq1[initial_i], seq2[initial_j]);
                delete = dp[initial_i-1][initial_j] + GAP_SCORE;
                insert = dp[initial_i][initial_j-1] + GAP_SCORE;
                max_score = max3(match, delete, insert);
            }

            // Ensure non-negative, simplified
            dp[initial_i][initial_j] = (max_score > 0) ? max_score : 0;
        }

        // Move to next element in diagonal
        initial_i--;
        initial_j++;
        count++;
    }
}

void smith_waterman(uint32_t core_id) {
  int total_diagonals = seq1_len + seq2_len - 1;
  int roof_start=mempool_get_timer();
  for (int diagonal_index = 0; diagonal_index < total_diagonals; diagonal_index++) {
    int num_diagonal_elements = diagonal_length_cal(diagonal_index);

    int active_cores = minimum(num_diagonal_elements, cc);
    int elements_per_core = (num_diagonal_elements + active_cores - 1) / active_cores;

    compute_start_times[core_id] = 0;
    compute_end_times[core_id] = 0;

    if (core_id < active_cores) {
      int start = core_id * elements_per_core;
      int end = (core_id == active_cores - 1) ? num_diagonal_elements : (core_id + 1) * elements_per_core;

      compute_start_times[core_id] = mempool_get_timer();
      compute_diagonal_segment(diagonal_index, start, end);
      compute_end_times[core_id] = mempool_get_timer();
    }

    mempool_barrier(cc);
    //
    // if (core_id == 0) {
    //   uint32_t earliest_start = UINT32_MAX;
    //   uint32_t latest_end = 0;
    //
    //   for (int i = 0; i < active_cores; i++) {
    //     if (compute_start_times[i] < earliest_start) {
    //       earliest_start = compute_start_times[i];
    //     }
    //     if (compute_end_times[i] > latest_end) {
    //       latest_end = compute_end_times[i];
    //     }
    //   }
    //
    //   uint32_t diagonal_compute_time = latest_end - earliest_start;
    //   uint32_t max_individual_compute = 0;
    //   for (int i = 0; i < active_cores; i++) {
    //     uint32_t core_compute_time = compute_end_times[i] - compute_start_times[i];
    //     if (core_compute_time > max_individual_compute) {
    //       max_individual_compute = core_compute_time;
    //     }
    //   }
    //
    //   uint32_t diagonal_sync_time = diagonal_compute_time - max_individual_compute;
    //
    //   global_compute_time += max_individual_compute;
    //   global_sync_time += diagonal_sync_time;
    //
    //   printf("Diagonal %d - Compute time: %u, Sync time: %u\n", diagonal_index, max_individual_compute, diagonal_sync_time);
    // }
     mempool_barrier(cc);
  }
  int roof_end=mempool_get_timer();
  if(core_id==0){
    printf("the roofstart and roof end is %d\n%d\n",roof_start,roof_end);
  }

}

int main() {
  uint32_t core_id = mempool_get_core_id();

  mempool_barrier_init(core_id);

  if (core_id >= cc) {
    while (1);
  }

  if (core_id == 0) {
    memset(dp, 0, sizeof(dp));

    printf("Sequence 1: %s (Length: %d)\n", seq1, seq1_len);
    printf("Sequence 2: %s (Length: %d)\n", seq2, seq2_len);
  }

  mempool_barrier(cc);

  uint32_t start_cycles = mempool_get_timer();
  smith_waterman(core_id);
  uint32_t end_cycles = mempool_get_timer();
  uint32_t cycle_count = end_cycles - start_cycles;

  mempool_barrier(cc);

  if (core_id == 0) {
    printf("The cycles taken for the complete smith_waterman kernel: %d\n", cycle_count);
    printf("End-Start: =%d-%d\n", end_cycles, start_cycles);
    printf("Cycles taken for compute: %d\n", global_compute_time);
    printf("Cycles taken for synchronization: %d\n", global_sync_time);
    printf("Total cell updates made: %d\n", global_cell_updates);
  }

  mempool_barrier(cc);

  return 0;
}
