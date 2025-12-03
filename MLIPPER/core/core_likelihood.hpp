#ifndef ML_CORE_HPP
#define ML_CORE_HPP

#include <cstddef>   // size_t
#include <cstdint>   // uint*_t
#include <climits>   // UINT_MAX
#define SCALER_THRESHOLD 1e-300
#define PLL_SCALE_RATE_MAXDIFF 100
#define PLL_SCALE_THRESHOLD 1e+6
namespace Core_Likelihood
{
/* =======================================================
 * 1) Non-root internal node CLV computation
 * node_clv       : OUTPUT buffer [sites*rate_cats*states]
 * child_clv_L/R  : INPUT  buffers (same layout)
 * pmatrix_L/R    : INPUT  [rate_cats*states*states]
 * scaler         : OUTPUT site-wise scaling counters [sites] (nullable)
 * This function does not allocate or free memory.
 */
double core_edge_loglikelihood_inner_inner(
    std::size_t sites,
    std::size_t states,
    std::size_t rate_cats,
    const double * parent_clv,
    const unsigned int * parent_scaler,
    const double * child_clv,
    const unsigned int * child_scaler,
    const double * pmatrix,
    double * const * frequencies,
    const double * rate_weights,
    const unsigned int * pattern_weights,
    const double * invar_proportion,
    const int * invar_indices,
    const unsigned int * freqs_indices,
    double * persite_lnl
);

/* =======================================================
 * 2) Tip–Inner edge log-likelihood
 * parent_clv      : INPUT  [sites*rate_cats*states]
 * parent_scaler   : INPUT  [sites] (nullable)
 * tipchars        : INPUT  observed bases as ASCII ('A','C','G','T','U', etc.)
 *                   (If you use encoded 0..3 instead, rename and adjust doc.)
 * pmatrix         : INPUT  [rate_cats*states*states]
 * frequencies     : INPUT  array of π vectors (deep-const):
 *                          frequencies[freqs_indices[j]][k]
 * rate_weights    : INPUT  [rate_cats]
 * pattern_weights : INPUT  [sites]
 * invar_proportion: INPUT  nullable (indexed by freqs_indices[j])
 * invar_indices   : INPUT  nullable [sites] (-1 for ambiguous)
 * freqs_indices   : INPUT  [rate_cats] maps rate→which π vector
 * persite_lnl     : OUTPUT [sites] (nullable)
 * Returns total log-likelihood (sum over sites with weights).
 */
double core_edge_loglikelihood_tip_inner(
    std::size_t sites,
    std::size_t rate_cats,
    const double * parent_clv,                 // [sites][rate_cats][states]
    const unsigned int * parent_scaler,        // nullable, per-site
    const unsigned char * tipchars,            // child tip code
    const unsigned int  * tipmap,              // code -> bitmask over states (nullable => unambiguous)
    const double * pmatrix,                    // [rate_cats][states][states] row-major
    const double * const * frequencies,        // frequencies[mix_id][state]
    const double * rate_weights,               // len=rate_cats
    const unsigned int * pattern_weights,      // len=sites
    const double * invar_proportion,           // nullable, per-mixture
    const int * invar_indices,                 // nullable, per-site; -1 if not invariant
    const unsigned int * freqs_indices,        // len=rate_cats, maps rate->mixture
    double * persite_lnl,                      // OUTPUT or nullptr
    std::size_t states
);

/* =======================================================
 * 3) Root log-likelihood
 * root_clv        : INPUT  [sites*rate_cats*states]
 * site_scaler     : INPUT  [sites] (nullable)
 * frequencies     : INPUT  array of π vectors (deep-const)
 * rate_weights    : INPUT  [rate_cats]
 * pattern_weights : INPUT  [sites]
 * invar_proportion: INPUT  nullable (indexed by freqs_indices[j])
 * invar_indices   : INPUT  nullable [sites]
 * freqs_indices   : INPUT  [rate_cats]
 * persite_lnl     : OUTPUT [sites] (nullable)
 * Returns total log-likelihood.
 */
double core_edge_loglikelihood_root(
    std::size_t sites,
    std::size_t rate_cats,
    const double * root_clv,
    const unsigned int * site_scaler,           // nullable
    const double* frequencies,         
    const double * rate_weights,
    const unsigned int * pattern_weights,
    const double invar_proportion,            // nullable
    const int * invar_indices,                  // nullable
    //double * persite_lnl,                       // OUTPUT or nullptr
    std::size_t states
);

} // namespace mlcore

#endif // ML_CORE_HPP
