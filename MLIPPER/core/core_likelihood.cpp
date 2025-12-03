#include "core_likelihood.hpp"
#include <cmath>    
#include <cstring>  
#include <algorithm> 

double core_likelihood::core_edge_loglikelihood_inner_inner(
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
){
    std::size_t n,i,j,k;
    double logl = 0;

    const double * clvp = parent_clv;
    const double * clvc = child_clv;
    double prop_invar = 0;
    const double * pmat;
    const double * freqs = NULL;

    double terma, terma_r, termb;
    double site_lk, inv_site_lk;

    /* TODO: We need states_padded in the AVX/SSE implementations
    */
    unsigned int states_padded = states;
    unsigned int site_scalings;
    unsigned int * rate_scalings = NULL;

    /* powers of scale threshold for undoing the scaling */
    double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
    scale_factor *= PLL_SCALE_THRESHOLD;
    scale_minlh[i] = scale_factor;
    }

    for (n = 0; n < sites; ++n)
    {
        site_scalings = UINT_MAX;
        for (i = 0; i < rate_cats; ++i){
            rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
            rate_scalings[i] += (child_scaler) ? child_scaler[n*rate_cats+i] : 0;
            if (rate_scalings[i] < site_scalings)
            site_scalings = rate_scalings[i];
        }

        /* compute relative capped per-rate scalers */
        for (i = 0; i < rate_cats; ++i){
            rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings, PLL_SCALE_RATE_MAXDIFF);
        }

        pmat = pmatrix;
        terma = 0;
        for (i = 0; i < rate_cats; ++i){
            freqs = frequencies[freqs_indices[i]];
            terma_r = 0;
            for (j = 0; j < states; ++j)
            {
                termb = 0;
                for (k = 0; k < states; ++k)
                {
                termb += pmat[k] * clvc[k];
                }

                terma_r += clvp[j] * freqs[j] * termb;
                pmat += states_padded;
            }

            /* apply per-rate scalers, if necessary */
            if (rate_scalings && rate_scalings[i] > 0){
                terma_r *= scale_minlh[rate_scalings[i]-1];
            }

            /* account for invariant sites */
            prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
            if (prop_invar > 0){
                inv_site_lk = (invar_indices[n] == -1) ? 0 : freqs[invar_indices[n]];
                terma += rate_weights[i] * (terma_r * (1 - prop_invar) + inv_site_lk * prop_invar);
            }
            else{
                terma += terma_r * rate_weights[i];
            }

            clvp += states_padded;
            clvc += states_padded;
        }

        /* compute site log-likelihood and scale if necessary */
        site_lk = log(terma);
        if (site_scalings)
        site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);

        site_lk *= pattern_weights[n];

        /* store per-site log-likelihood */
        if (persite_lnl)
        persite_lnl[n] = site_lk;

        logl += site_lk;
    }

    if (rate_scalings) free(rate_scalings);

    return logl;
}
double core_likelihood::core_edge_loglikelihood_tip_inner(
    std::size_t sites,
    std::size_t states,
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
    double * persite_lnl                       // OUTPUT or nullptr
){
    double logli = 0; 
    double site_lk;
    double term_var;
    double termb, terma_r, terma;
    std::size_t i,j,k;
    double prop_invar;
    double inv_site_lk;
    const double * freqs;
    unsigned int site_scalings;
    unsigned int cstate;
    const double * pmat;
    unsigned int * rate_scalings = NULL;
    const double * clvp = parent_clv;

    /* powers of scale threshold for undoing the scaling */
    double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
        scale_factor *= PLL_SCALE_THRESHOLD;
        scale_minlh[i] = scale_factor;
    }

        
    for(i = 0; i < sites; i++){
        site_scalings = UINT_MAX;
        site_lk = 0;
        for (j = 0; j < rate_cats; ++j)
        {
            rate_scalings[j] = (parent_scaler) ? parent_scaler[i * rate_cats+j] : 0;
            if (rate_scalings[j] < site_scalings) site_scalings = rate_scalings[j];
        }

        /* compute relative capped per-rate scalers */
        for (j = 0; j < rate_cats; ++j)
        {
            rate_scalings[j] = std::min(rate_scalings[j] - site_scalings, static_cast<unsigned int>(PLL_SCALE_RATE_MAXDIFF));
        }
        pmat = pmatrix;
        terma = 0.0;
        for (j = 0; j < rate_cats; j++){
            freqs = frequencies[freqs_indices[j]];
            prop_invar = invar_proportion[freqs_indices[j]];
            terma_r = 0;
            for (j = 0; j < states; ++j)
            {
                termb = 0;
                cstate = tipmap[(unsigned int)(*tipchars)];

                for (k = 0; k < states; ++k)
                {
                    if (cstate & 1) termb += pmat[k];
                    cstate >>= 1;
                }
                terma_r += clvp[j] * freqs[j] * termb;
                pmat += states;
            }


            /* apply per-rate scalers, if necessary */
            if (rate_scalings && rate_scalings[i] > 0)
            {
                terma_r *= scale_minlh[rate_scalings[i]-1];
            }

            /* account for invariant sites */
            if (prop_invar > 0)
            {
                inv_site_lk = (invar_indices[i] == -1) ? 0 : freqs[invar_indices[i]];
                terma += rate_weights[i] * (terma_r * (1 - prop_invar) + inv_site_lk * prop_invar);
            }
            else
            {
                terma += terma_r * rate_weights[i];
            }
            clvp += states; 
        }
        
        /* compute site log-likelihood and scale if necessary */
        site_lk = log(terma);
        if (site_scalings) site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);

        site_lk *= pattern_weights[i];

        /* store per-site log-likelihood */
        if (persite_lnl)
        persite_lnl[i] = site_lk;

        logli += site_lk;

        tipchars++;
        
                
    }
    if (rate_scalings){
        free(rate_scalings);
    }
    return logli;
}


double core_likelihood::core_edge_loglikelihood_root(
    std::size_t sites,
    std::size_t states,
    std::size_t rate_cats,
    const double * root_clv,
    const unsigned int * site_scaler,           // nullable
    const double* frequencies,         
    const double * rate_weights,
    const unsigned int * pattern_weights,
    const double invar_proportion,            
    const int * invar_indices                
){
    double logli=0.0, site_lk = 0.0;
    double term,term_var;
    double inv_site_lk;
    std::size_t i,j,k;

    for(i = 0; i < sites; i++){
        term = 0;
        for(j = 0; j < rate_cats; j++){
            term_var = 0.0;
            for(k =0; k < states; k++){

                term_var += root_clv[k] * frequencies[k];
            }
            inv_site_lk = (invar_indices[i] == -1) ? 0 : frequencies[invar_indices[i]];
            if(invar_proportion > 0){
                term += (1.0 - invar_proportion) * term_var + invar_proportion * inv_site_lk;
            } else {
                term += term_var * rate_weights[j];
            }
        }
        site_lk = log(term) + site_scaler[i]*log(SCALER_THRESHOLD);
        logli += site_lk * pattern_weights[i];
        root_clv += states;
    }
    return logli;
}