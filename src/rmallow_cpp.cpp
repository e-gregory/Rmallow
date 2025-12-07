// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

//' Compute pairwise comparison information for Kendall's distance (C++ version)
//'
//' Fast C++ implementation that computes pairwise comparisons for Kendall distance.
//' Returns 0 for increases, 1 for decreases, and NA for ties.
//'
//' @param r NumericMatrix of rankings (each row is a ranking)
//' @param inds IntegerMatrix of column index pairs (2 x n_pairs), 1-indexed
//' @return NumericMatrix of comparison info (n_obs x n_pairs)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix kendall_info_cpp(NumericMatrix r, IntegerMatrix inds) {
    int n_obs = r.nrow();
    int n_pairs = inds.ncol();
    
    NumericMatrix infos(n_obs, n_pairs);
    
    // Fill with NA initially
    std::fill(infos.begin(), infos.end(), NA_REAL);
    
    // Process each pair comparison
    for (int j = 0; j < n_pairs; j++) {
        int col1 = inds(0, j) - 1;  // Convert to 0-indexed
        int col2 = inds(1, j) - 1;
        
        for (int i = 0; i < n_obs; i++) {
            double diff = r(i, col1) - r(i, col2);
            if (diff > 0) {
                infos(i, j) = 1.0;  // Decrease
            } else if (diff < 0) {
                infos(i, j) = 0.0;  // Increase
            }
            // NA remains for ties (diff == 0)
        }
    }
    
    return infos;
}

//' Generate pairwise column indices
//'
//' Fast C++ implementation to generate all (i, j) pairs where i < j.
//'
//' @param n Number of columns
//' @return IntegerMatrix of size 2 x (n choose 2)
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix combn_pairs_cpp(int n) {
    int n_pairs = n * (n - 1) / 2;
    IntegerMatrix inds(2, n_pairs);
    
    int k = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = i + 1; j <= n; j++) {
            inds(0, k) = i;
            inds(1, k) = j;
            k++;
        }
    }
    
    return inds;
}

//' Compute Kendall distances between rankings and reference sequences (C++ version)
//'
//' Fast C++ implementation that computes all pairwise Kendall distances.
//'
//' @param data_info NumericMatrix of Kendall info for data (n_obs x n_pairs)
//' @param seqs_info NumericMatrix of Kendall info for sequences (n_seqs x n_pairs)
//' @return NumericMatrix of distances (n_obs x n_seqs)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix all_kendall_cpp(NumericMatrix data_info, NumericMatrix seqs_info) {
    int n_obs = data_info.nrow();
    int n_seqs = seqs_info.nrow();
    int n_pairs = data_info.ncol();
    
    NumericMatrix dists(n_obs, n_seqs);
    
    for (int i = 0; i < n_obs; i++) {
        for (int j = 0; j < n_seqs; j++) {
            double dist = 0.0;
            for (int k = 0; k < n_pairs; k++) {
                double d = data_info(i, k);
                double s = seqs_info(j, k);
                // Skip if either is NA (tie)
                if (!NumericMatrix::is_na(d) && !NumericMatrix::is_na(s)) {
                    dist += std::abs(d - s);
                }
            }
            dists(i, j) = dist;
        }
    }
    
    return dists;
}

//' Compute distances to canonical ordering (C++ version)
//'
//' Counts the number of inversions (1s) in Kendall info matrix.
//'
//' @param infos NumericMatrix of Kendall info
//' @return IntegerVector of distances
//' @keywords internal
// [[Rcpp::export]]
IntegerVector all_seq_dists_cpp(NumericMatrix infos) {
    int n_obs = infos.nrow();
    int n_pairs = infos.ncol();
    
    IntegerVector dists(n_obs);
    
    for (int i = 0; i < n_obs; i++) {
        int count = 0;
        for (int j = 0; j < n_pairs; j++) {
            double val = infos(i, j);
            if (!NumericMatrix::is_na(val) && val == 1.0) {
                count++;
            }
        }
        dists[i] = count;
    }
    
    return dists;
}

//' Compute weighted sums for UpdateR (C++ version)
//'
//' Computes weighted sums of membership probabilities for 0s and 1s
//' in the Kendall info matrix.
//'
//' @param infos NumericMatrix of Kendall info (n_obs x n_pairs)
//' @param z NumericMatrix of membership probabilities (n_obs x G)
//' @return List with zero_sums and one_sums matrices (G x n_pairs)
//' @keywords internal
// [[Rcpp::export]]
List update_r_sums_cpp(NumericMatrix infos, NumericMatrix z) {
    int n_obs = infos.nrow();
    int n_pairs = infos.ncol();
    int G = z.ncol();
    
    NumericMatrix zero_sums(G, n_pairs);
    NumericMatrix one_sums(G, n_pairs);
    
    // Initialize to zero
    std::fill(zero_sums.begin(), zero_sums.end(), 0.0);
    std::fill(one_sums.begin(), one_sums.end(), 0.0);
    
    // Process each pair comparison
    for (int j = 0; j < n_pairs; j++) {
        for (int i = 0; i < n_obs; i++) {
            double val = infos(i, j);
            if (!NumericMatrix::is_na(val)) {
                if (val == 0.0) {
                    // Add membership probs for this obs to zero_sums
                    for (int g = 0; g < G; g++) {
                        zero_sums(g, j) += z(i, g);
                    }
                } else if (val == 1.0) {
                    // Add membership probs for this obs to one_sums
                    for (int g = 0; g < G; g++) {
                        one_sums(g, j) += z(i, g);
                    }
                }
            }
        }
    }
    
    return List::create(
        Named("zero_sums") = zero_sums,
        Named("one_sums") = one_sums
    );
}

//' E-Step computation (C++ version)
//'
//' Computes unnormalized log probabilities and applies softmax normalization.
//'
//' @param all_dists NumericMatrix of distances (n_obs x G)
//' @param lambda NumericVector of spread parameters (length G)
//' @param C NumericVector of normalizing constants (length G)
//' @param p NumericVector of cluster proportions (length G)
//' @return NumericMatrix of membership probabilities (n_obs x G)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix estep_cpp(NumericMatrix all_dists, NumericVector lambda, 
                        NumericVector C, NumericVector p) {
    int N = all_dists.nrow();
    int G = all_dists.ncol();
    
    NumericMatrix z(N, G);
    
    // Precompute log(C * p)
    NumericVector log_cp(G);
    for (int g = 0; g < G; g++) {
        log_cp[g] = std::log(C[g] * p[g]);
    }
    
    for (int i = 0; i < N; i++) {
        // Compute log probabilities
        double max_log_prob = R_NegInf;
        for (int g = 0; g < G; g++) {
            double log_prob = -lambda[g] * all_dists(i, g) + log_cp[g];
            z(i, g) = log_prob;
            if (log_prob > max_log_prob) {
                max_log_prob = log_prob;
            }
        }
        
        // Subtract max for numerical stability and exponentiate
        double sum_exp = 0.0;
        for (int g = 0; g < G; g++) {
            z(i, g) = std::exp(z(i, g) - max_log_prob);
            sum_exp += z(i, g);
        }
        
        // Normalize
        if (sum_exp > 0) {
            for (int g = 0; g < G; g++) {
                z(i, g) /= sum_exp;
            }
        } else {
            // Handle edge case
            for (int g = 0; g < G; g++) {
                z(i, g) = 1.0 / G;
            }
        }
    }
    
    return z;
}

//' NextTable computation (C++ version)
//'
//' Computes distance distribution for (N+1)! from N! space using recurrence.
//'
//' @param last_table NumericVector of distance counts in N! space
//' @param n_last The value of N
//' @return NumericVector of distance counts in (N+1)! space
//' @keywords internal
// [[Rcpp::export]]
NumericVector next_table_cpp(NumericVector last_table, int n_last) {
    int last_len = last_table.size();
    
    // Maximum distance in (N+1)! space is N*(N+1)/2
    int new_len = ((n_last + 1) * n_last) / 2 + 1;
    
    // Compute cumulative sum of last_table
    NumericVector cumsum_last(last_len);
    cumsum_last[0] = last_table[0];
    for (int i = 1; i < last_len; i++) {
        cumsum_last[i] = cumsum_last[i - 1] + last_table[i];
    }
    
    NumericVector result(new_len);
    
    // First n+1 entries: cumsum of first (n+1) entries
    for (int i = 0; i <= n_last && i < last_len; i++) {
        result[i] = cumsum_last[i];
    }
    
    // Middle section: sliding window sums of width (n+1)
    int middle_len = new_len - 2 * (n_last + 1);
    if (middle_len > 0) {
        for (int i = 0; i < middle_len; i++) {
            int end_idx = n_last + 1 + i;
            int start_idx = i;
            result[n_last + 1 + i] = cumsum_last[end_idx] - cumsum_last[start_idx];
        }
    }
    
    // Last n+1 entries: reverse of first part
    for (int i = 0; i <= n_last; i++) {
        result[new_len - 1 - i] = result[i];
    }
    
    return result;
}

//' SimplifySequences (C++ version)
//'
//' Transforms rankings so that tie groups use consecutive integers.
//'
//' @param rankings NumericMatrix of rankings
//' @return NumericMatrix of simplified rankings
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix simplify_sequences_cpp(NumericMatrix rankings) {
    int n_rows = rankings.nrow();
    int n_cols = rankings.ncol();
    
    NumericMatrix result(n_rows, n_cols);
    
    for (int i = 0; i < n_rows; i++) {
        // Get unique sorted values
        std::vector<double> row_vals;
        for (int j = 0; j < n_cols; j++) {
            row_vals.push_back(rankings(i, j));
        }
        
        std::vector<double> unique_vals = row_vals;
        std::sort(unique_vals.begin(), unique_vals.end());
        unique_vals.erase(std::unique(unique_vals.begin(), unique_vals.end()), unique_vals.end());
        
        // Create mapping and apply
        for (int j = 0; j < n_cols; j++) {
            auto it = std::lower_bound(unique_vals.begin(), unique_vals.end(), row_vals[j]);
            result(i, j) = std::distance(unique_vals.begin(), it) + 1;
        }
    }
    
    return result;
}

//' ConstructSeqs (C++ version)
//'
//' Reconstructs ranking sequences from their Kendall preference vectors.
//'
//' @param prefs NumericMatrix of pairwise preferences
//' @param n_abils Number of items in the ranking
//' @return List of integer vectors (rankings)
//' @keywords internal
// [[Rcpp::export]]
List construct_seqs_cpp(NumericMatrix prefs, int n_abils) {
    int n_seqs = prefs.nrow();
    
    // Precompute top positions (where each item's comparisons start)
    std::vector<int> tops(n_abils);
    tops[0] = 0;
    for (int i = 1; i < n_abils; i++) {
        tops[i] = tops[i - 1] + (n_abils - i);
    }
    
    List R(n_seqs);
    
    for (int seq_idx = 0; seq_idx < n_seqs; seq_idx++) {
        // Initialize available items
        std::vector<int> nums(n_abils);
        for (int i = 0; i < n_abils; i++) {
            nums[i] = i + 1;
        }
        
        IntegerVector result(n_abils);
        
        for (int i = 0; i < n_abils - 1; i++) {
            int start_col = tops[i];
            int end_col = (i + 1 < n_abils) ? tops[i + 1] - 1 : prefs.ncol() - 1;
            
            // Sum preferences for this position
            double sumz = 0;
            for (int c = start_col; c <= end_col; c++) {
                sumz += prefs(seq_idx, c);
            }
            
            // Position in remaining items (0-indexed)
            int pos = (int)sumz;
            int item = nums[pos];
            
            result[i] = item;
            
            // Remove used item
            nums.erase(nums.begin() + pos);
        }
        
        // Add final item
        result[n_abils - 1] = nums[0];
        
        R[seq_idx] = result;
    }
    
    return R;
}
