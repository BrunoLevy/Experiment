#ifndef EXPERIMENT_ALGO_CDT_WITH_INTERVAL
#define EXPERIMENT_ALGO_CDT_WITH_INTERVAL

#include <OGF/Experiment/common/common.h>
#include <geogram/delaunay/CDT_2d.h>

// Testing interval arithmetics
// Using https://www.codeproject.com/Articles/1040839/Interval-arithmetic-using-round-to-nearest-mode-pa
// Tested also my own implementation / simplification in algo/interval_nt.h
// Both fail with constraints_100_10.obj (copied here: PR21/bug_CDT2d)
// This one seems to work: https://accu.org/journals/overload/19/103/harris_1974/

namespace GEO {

    class GEOGRAM_API CDT2d_with_interval: public CDT2d {
    public:
        CDT2d_with_interval() {
            orient_cnt_fail_=0;
            orient_cnt_total_=0;
            incircle_cnt_fail_=0;
            incircle_cnt_total_=0;
        }
        ~CDT2d_with_interval() {
            std::cerr << "Orient   needed exact: " << orient_cnt_fail_ << " / " << orient_cnt_total_ << std::endl;
            std::cerr << "Incircle needed exact: " << incircle_cnt_fail_ << " / " << incircle_cnt_total_ << std::endl;        
        }
    protected:
        Sign orient2d(index_t i, index_t j, index_t k) const override;
        Sign incircle(index_t i,index_t j,index_t k,index_t l) const override;

        mutable index_t orient_cnt_fail_;
        mutable index_t orient_cnt_total_;
        mutable index_t incircle_cnt_fail_;
        mutable index_t incircle_cnt_total_;
    };
    
    
}

#endif
