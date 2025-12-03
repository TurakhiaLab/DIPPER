#pragma once
#ifndef CORE_LIKELIHOOD_CUH
#define CORE_LIKELIHOOD_CUH
#include <cstddef>

#define SCALE_FACTOR 18446744073709551616  /*  2**64 (exactly)  */
#define SCALE_THRESHOLD 5.4210108624275221700372640043e-20
#define LOG_SCALE_THRESHOLD -44.361419555836499
#define SCALE_MAX_DIFF 32  // PLL_SCALE_RATE_MAXDIFF
#define MAX_STATES 64
#define MAX_RATECATS 8

namespace core_likelihood {

    // extern __constant__ double c_frequencies[MAX_STATES];
    // extern __constant__ double c_rate_weights[MAX_RATECATS];
    // extern __constant__ double c_scale_minlh[SCALE_MAX_DIFF];
    // extern __constant__ unsigned int c_tipmap[256];

    void InitializeConstants(const double* freqs,
                           const double* rates,
                           const double* scales,
                           const unsigned* tipmap);

    struct SumAndStore {
        const unsigned* parent;
        const unsigned* child;
        unsigned*       out;

        __host__ __device__
        unsigned operator()(unsigned idx) const {
            unsigned p = parent ? parent[idx] : 0;
            unsigned c = child  ? child[idx]  : 0;
            unsigned v = p + c;
            if (out) out[idx] = v;
            return v;             
        }
    };
    struct KeyBySite {
        unsigned rate_cats;
        __host__ __device__
        unsigned operator()(unsigned idx) const {
            return idx / rate_cats;
        }
    };

    struct Param {
        std::size_t sites;
        int states;
        int rate_cats;
        double invar_proportion;
        bool per_rate_scaling;
        
        Param() : sites(0), states(0), rate_cats(0), invar_proportion(0.0), per_rate_scaling(false){};

        Param(std::size_t t_sites, int t_states, int t_rate_cats, double t_invar_proportion, bool t_per_rate_scaling){
            sites = t_sites;
            states = t_states;
            rate_cats = t_rate_cats;
            invar_proportion = t_invar_proportion;
            per_rate_scaling = t_per_rate_scaling;
        };
    };

    struct Likelihood_Root {

        double* d_per_site;   
        double* d_out_sum;

        double* d_root_clv;
        double* d_frequencies; 
        double* d_rate_weights; 
        unsigned int* d_pattern_weights; 
        unsigned int* d_site_scaler;
        int* d_invar_indices; 

        double* total_loglik;

        void Initialize(const Param& p);
        void ConstructionOnGpu(
            const double* h_root_clv,
            const double* h_frequencies,
            const double* h_rate_weights,
            const unsigned int* h_pattern_weights,
            const unsigned int* h_site_scaler,
            const int*    h_invar_indices,
            const Param& p
        );
        void ComputeLikelihood(const Param& p);     
        void PrintLikelihood() const;                
        void CleanUp();
    };

    struct Likelihood_Tip_Inner {

        double* d_per_site;   
        double* d_out_sum;

        double* d_parent_clv;
        unsigned int* d_parent_scaler; // Can be sites or sites * rate_cats
        unsigned char* d_tipchars;
        unsigned int* d_tipmap; // Maps char to state index, size 256
        double* d_p_matrix;
        double* d_frequencies; 
        double* d_rate_weights; 
        unsigned int* d_pattern_weights; 
        int*  d_invar_indices;

        unsigned int* d_scale_minlh; 
        unsigned int* d_site_scaling_min;

        double total_loglik; 

        void Initialize(const Param& p);
        void ConstructionOnGpu( 
            const double* h_parent_clv,
            const unsigned int* h_parent_site_scaler,
            const unsigned char* h_tipchars,
            const unsigned int* h_tipmap,
            const double* h_p_matrix,
            const double* h_frequencies,
            const double* h_rate_weights,
            const unsigned int* h_pattern_weights,
            const int*    h_invar_indices,
            const Param& p
        );
        void ComputeLikelihood(const Param& p);
        void PrintLikelihood() const;
        void CleanUp();
    };

    struct Likelihood_Inner_Inner {

        double* d_per_site;   
        double* d_out_sum;

        double* d_parent_clv;
        unsigned int* d_parent_scaler; // Can be sites or sites * rate_cats
        double* d_child_clv;
        unsigned int* d_child_scaler; // Can be sites or sites * rate_c
        unsigned int* d_site_scaler;
        
        double* d_p_matrix;
        double* d_frequencies; 
        double* d_rate_weights; 
        unsigned int* d_pattern_weights; 
        int*  d_invar_indices;

        unsigned int* d_scale_minlh; 
        unsigned int* d_site_scaling_min;

        double total_loglik; 

        void Initialize( Param& p);
        void ConstructionOnGpu(
            const double * h_parent_clv,
            const unsigned int * h_parent_scaler,
            const double * h_child_clv,
            const unsigned int * h_child_scaler,
            const double * h_p_matrix,
            double * h_frequencies,
            const double * h_rate_weights,
            const unsigned int * h_pattern_weights,
            const int * h_invar_indices,
            Param& p);
        void ComputeLikelihood(Param& p);
        void PrintLikelihood() const;
        void CleanUp();
    };

    double ComputeRootLogLikelihoodFromDevice(
        const double* d_root_clv,
        std::size_t sites,
        int states,
        int rate_cats,
        double invar_proportion,
        const double* h_frequencies,
        const double* h_rate_weights,
        const unsigned int* d_pattern_weights,
        const unsigned int* d_site_scaler,
        const int* d_invar_indices,
        cudaStream_t stream = 0);

} // namespace core_likelihood

#endif // CORE_LIKELIHOOD_CUH
