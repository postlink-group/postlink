// mixbayes_labels.cpp
//
// Label-alignment helpers for the S x N matrix of posterior allocations
// (values 1 = correct match, 2 = mismatch) returned by the Gibbs sampler.
// They work on the matrix in place of R code that would create several S x N
// temporaries, so that the post-processing needs little memory beyond the
// allocation matrix itself.

#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// ECR-ITERATIVE-1 for two components (see ecr_iterative_two() in
// R/mixbayes_helpers.R): each draw is either kept or swapped so as to
// minimise its disagreements with a pivot allocation, and the pivot is
// re-estimated as the modal allocation under the current permutations until
// the total cost stops decreasing. Only the columns `cols` (1-based; all
// columns when `cols` is empty) are used, without copying them. Returns TRUE
// for the draws to swap.
//' @noRd
// [[Rcpp::export]]
LogicalVector ecr_two_cpp(const IntegerMatrix& z, int maxiter, double threshold,
                          const IntegerVector& cols) {
  const R_xlen_t S = z.nrow(), Nz = z.ncol();
  const int* zp = INTEGER(z);
  const R_xlen_t N = cols.size() > 0 ? cols.size() : Nz;
  std::vector<R_xlen_t> colj(N);
  for (R_xlen_t k = 0; k < N; ++k) {
    colj[k] = cols.size() > 0 ? static_cast<R_xlen_t>(cols[k]) - 1 : k;
    if (colj[k] < 0 || colj[k] >= Nz) stop("Column index out of range.");
  }
  std::vector<char> swap(S, 0), new_swap(S, 0), pivot2(N, 0);
  std::vector<R_xlen_t> mism(S);
  double cost_prev = R_PosInf;
  for (int iter = 0; iter < maxiter; ++iter) {
    // pivot: modal label of each record under the current permutations (ties -> label 1)
    for (R_xlen_t j = 0; j < N; ++j) {
      const int* col = zp + colj[j] * S;
      R_xlen_t n2 = 0;
      for (R_xlen_t i = 0; i < S; ++i) n2 += ((col[i] == 2) != (swap[i] != 0));
      pivot2[j] = (n2 > S - n2);
    }
    // disagreements of each raw draw with the pivot
    std::fill(mism.begin(), mism.end(), 0);
    for (R_xlen_t j = 0; j < N; ++j) {
      const int* col = zp + colj[j] * S;
      const bool p2 = pivot2[j] != 0;
      for (R_xlen_t i = 0; i < S; ++i) mism[i] += ((col[i] == 2) != p2);
    }
    double cost = 0.0;
    for (R_xlen_t i = 0; i < S; ++i) {
      new_swap[i] = mism[i] > N - mism[i];          // swapping costs N - mism
      cost += static_cast<double>(std::min(mism[i], N - mism[i]));
    }
    if (cost > cost_prev) break;                    // keep the previous permutations
    swap = new_swap;
    if (cost_prev - cost <= threshold) break;
    cost_prev = cost;
  }
  LogicalVector out(S);
  for (R_xlen_t i = 0; i < S; ++i) out[i] = swap[i] != 0;
  return out;
}

// Number of cells equal to 2 in the columns `cols` (1-based; all columns when
// `cols` is empty) after exchanging the labels of the draws flagged in `swap`.
//' @noRd
// [[Rcpp::export]]
double count_label2_cpp(const IntegerMatrix& z, const LogicalVector& swap, const IntegerVector& cols) {
  const R_xlen_t S = z.nrow(), N = z.ncol();
  if (swap.size() != S) stop("`swap` must have one element per draw.");
  const int* zp = INTEGER(z);
  const R_xlen_t nc = cols.size() > 0 ? cols.size() : N;
  double count = 0.0;
  for (R_xlen_t k = 0; k < nc; ++k) {
    const R_xlen_t j = cols.size() > 0 ? static_cast<R_xlen_t>(cols[k]) - 1 : k;
    if (j < 0 || j >= N) stop("Column index out of range.");
    const int* col = zp + j * S;
    for (R_xlen_t i = 0; i < S; ++i) count += ((col[i] == 2) != (swap[i] == TRUE));
  }
  return count;
}

// Copy of z with the labels (1 <-> 2) of the draws flagged in `swap` exchanged.
//' @noRd
// [[Rcpp::export]]
IntegerMatrix swap_rows_cpp(const IntegerMatrix& z, const LogicalVector& swap) {
  const R_xlen_t S = z.nrow(), N = z.ncol();
  if (swap.size() != S) stop("`swap` must have one element per draw.");
  IntegerMatrix out(z.nrow(), z.ncol());
  const int* zp = INTEGER(z);
  int* op = INTEGER(out);
  for (R_xlen_t j = 0; j < N; ++j) {
    for (R_xlen_t i = 0; i < S; ++i) {
      const int v = zp[j * S + i];
      op[j * S + i] = (swap[i] == TRUE) ? 3 - v : v;
    }
  }
  if (z.hasAttribute("dimnames")) out.attr("dimnames") = z.attr("dimnames");
  return out;
}
