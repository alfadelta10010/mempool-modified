#include <stdint.h>
#include "encoding.h"
#include <stdio.h>
#include <stdbool.h>
#include "printf.h"
#include "runtime.h"
#include "synchronization.h"
#include <string.h>

// DEFINITIONS
#define MAX_SEQ_LENGTH 100
#define MATCH_SCORE 2
#define MISMATCH_SCORE -1
#define GAP_SCORE -2

// CORE_COUNT - adjust as needed
#define cc 4

// TIMERS
uint32_t total_compute_time = 0;
uint32_t total_sync_time = 0;

// SEQUENCE DATA
const char seq1[] = "CCGGGGGGATCACCGCGGGCCTCCGCCGGTAGCTTGCGAT";
const char seq2[] = "AATCAATGAACACCGCATACCCTATTGTTGAAAGCGCTACGATCGTGATTTACTATCGTGCTCGACTGCAATGAGCCTAACGAAAAGTCGTGCTAGAACT";

int seq1_len;
int seq2_len;

// DP matrix - note: we need (m+1)x(n+1) for proper Smith-Waterman
int16_t dp[MAX_SEQ_LENGTH + 1][MAX_SEQ_LENGTH + 1] __attribute__((section(".l1")));

// Function to compute score
int16_t score(char a, char b) {
    if (a == b) return MATCH_SCORE;
    return MISMATCH_SCORE;
}

// Function to find minimum
int16_t minimum(int16_t x, int16_t y) {
    return (x < y) ? x : y;
}

// Function to find maximum
int16_t maximum(int16_t x, int16_t y) {
    return (x > y) ? x : y;
}

// Function to compute diagonal length
int diagonal_length(int diag_idx, int m, int n) {
    if (diag_idx <= minimum(m, n)) {
        return diag_idx;
    } else if (diag_idx <= maximum(m, n)) {
        return minimum(m, n);
    } else {
        return m + n - diag_idx;
    }
}

// Function to get (i,j) coordinates for k-th element in diagonal
void get_diagonal_coords(int diag_idx, int k, int m, int n, int *i, int *j) {
    if (diag_idx <= n) {
        *i = k + 1;
        *j = diag_idx - k;
    } else {
        *i = k + diag_idx - n + 1;
        *j = n - k;
    }
}

// Compute a segment of diagonal elements
void compute_diagonal_segment(int diagonal_idx, int start_elem, int end_elem, int m, int n) {
    for (int elem = start_elem; elem < end_elem; elem++) {
        int i, j;
        get_diagonal_coords(diagonal_idx, elem, m, n, &i, &j);
        
        // Boundary check
        if (i > m || j < 1 || j > n) continue;
        
        int16_t match_score, delete_score, insert_score, max_score;
        
        if (i == 1 && j == 1) {
            // First cell
            match_score = score(seq1[i-1], seq2[j-1]);
            delete_score = 0;
            insert_score = 0;
        } else if (i == 1) {
            // First row
            match_score = score(seq1[i-1], seq2[j-1]);
            delete_score = 0;
            insert_score = dp[i][j-1] + GAP_SCORE;
        } else if (j == 1) {
            // First column
            match_score = score(seq1[i-1], seq2[j-1]);
            delete_score = dp[i-1][j] + GAP_SCORE;
            insert_score = 0;
        } else {
            // General case
            match_score = dp[i-1][j-1] + score(seq1[i-1], seq2[j-1]);
            delete_score = dp[i-1][j] + GAP_SCORE;
            insert_score = dp[i][j-1] + GAP_SCORE;
        }
        
        // Find maximum score
        max_score = maximum(match_score, maximum(delete_score, insert_score));
        max_score = maximum(max_score, 0); // Smith-Waterman: no negative scores
        
        dp[i][j] = max_score;
    }
}

void smith_waterman_parallel(uint32_t core_id, uint32_t num_cores) {
    int m = seq1_len;
    int n = seq2_len;
    int total_diagonals = m + n;
    
    uint32_t core_compute_time = 0;
    
    for (int diag_idx = 1; diag_idx <= total_diagonals; diag_idx++) {
        int diag_len = diagonal_length(diag_idx, m, n);
        
        if (diag_len == 0) continue;
        
        // Distribute elements among active cores
        int active_cores = minimum(diag_len, num_cores);
        
        if (core_id < active_cores) {
            int elements_per_core = (diag_len + active_cores - 1) / active_cores;
            int start_elem = core_id * elements_per_core;
            int end_elem = minimum(start_elem + elements_per_core, diag_len);
            
            uint32_t start_time = mempool_get_timer();
            compute_diagonal_segment(diag_idx, start_elem, end_elem, m, n);
            uint32_t end_time = mempool_get_timer();
            
            core_compute_time += (end_time - start_time);
        }
        
        // Synchronize all cores after each diagonal
        mempool_barrier(num_cores);
    }
    
    // Aggregate timing information
    if (core_id == 0) {
        total_compute_time = core_compute_time;
        printf("Core 0 compute time: %u cycles\n", core_compute_time);
    }
}

// Function to find the maximum score in the DP matrix
int16_t find_max_score() {
    int16_t max_score = 0;
    for (int i = 1; i <= seq1_len; i++) {
        for (int j = 1; j <= seq2_len; j++) {
            if (dp[i][j] > max_score) {
                max_score = dp[i][j];
            }
        }
    }
    return max_score;
}

// Function to print a portion of the DP matrix for debugging
void print_dp_matrix_portion(int max_i, int max_j) {
    printf("\nDP Matrix (first %dx%d):\n", max_i, max_j);
    for (int i = 0; i <= max_i && i <= seq1_len; i++) {
        for (int j = 0; j <= max_j && j <= seq2_len; j++) {
            printf("%3d ", dp[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main() {
    uint32_t core_id = mempool_get_core_id();
    uint32_t num_cores = mempool_get_core_count();
    
    // Initialize synchronization
    mempool_barrier_init(core_id);
    
    // Sleep unused cores
    if (core_id >= cc) {
        while(1) {
            mempool_wfi();
        }
    }
    
    // Initialize sequences and DP matrix
    if (core_id == 0) {
        seq1_len = strlen(seq1);
        seq2_len = strlen(seq2);
        
        printf("Initializing Smith-Waterman algorithm\n");
        printf("Sequence 1: %s (Length: %d)\n", seq1, seq1_len);
        printf("Sequence 2: %s (Length: %d)\n", seq2, seq2_len);
        printf("Using %d cores\n", cc);
        
        // Initialize DP matrix to zero
        memset(dp, 0, sizeof(dp));
        
        printf("DP matrix initialized\n");
    }
    
    mempool_barrier(cc);
    
    // Start timing
    uint32_t start_cycles = mempool_get_timer();
    
    // Execute Smith-Waterman algorithm
    smith_waterman_parallel(core_id, cc);
    
    uint32_t end_cycles = mempool_get_timer();
    
    mempool_barrier(cc);
    
    // Output results
    if (core_id == 0) {
        uint32_t total_cycles = end_cycles - start_cycles;
        
        printf("\n=== RESULTS ===\n");
        printf("Total execution time: %u cycles\n", total_cycles);
        printf("Compute time: %u cycles\n", total_compute_time);
        
        // Find and print maximum alignment score
        int16_t max_score = find_max_score();
        printf("Maximum alignment score: %d\n", max_score);
        
        // Print a portion of the DP matrix for verification
        print_dp_matrix_portion(minimum(10, seq1_len), minimum(10, seq2_len));
    }
    
    mempool_barrier(cc);
    return 0;
}
