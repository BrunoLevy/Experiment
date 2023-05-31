#ifndef EXPERIMENT_ALGO_CDT_WITH_INTERVAL
#define EXPERIMENT_ALGO_CDT_WITH_INTERVAL

#include <OGF/Experiment/common/common.h>
#include <geogram/delaunay/CDT_2d.h>

// Testing interval arithmetics
// Using https://www.codeproject.com/Articles/1040839/Interval-arithmetic-using-round-to-nearest-mode-pa
// Tested also my own implementation / simplification in algo/interval_nt.h
// Both fail with constraints_100_10.obj

namespace GEO {

    class GEOGRAM_API CDT2d_with_interval: public CDT2d {
    public:
    CDT2d_with_interval() : cnt_fail_(0) {
    }
    ~CDT2d_with_interval() {
        std::cerr << "Filter fail: " << cnt_fail_ << std::endl;
    }
    protected:
        Sign orient2d(index_t i, index_t j, index_t k) const override;
        Sign incircle(index_t i,index_t j,index_t k,index_t l) const override;

        mutable index_t cnt_fail_;
    };
    
    
}

#endif
