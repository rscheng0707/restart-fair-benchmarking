// bbob_adapter.hpp — IOHexperimenter BBOB suite behind the tfb::Problem interface
//
// WHY THIS EXISTS
// The hand-written functions in benchmarks.hpp were adequate for harness
// development but are unsuitable for the main matrix, for two reasons:
//   (1) Reviewers expect the community-standard BBOB suite with its official
//       instance transformations and target values, not hand-rolled Rastrigin.
//   (2) They are too CHEAP. A d-dimensional sphere is O(d) arithmetic. Pilot 2
//       produced runs of 20-140 ms, at which scale OpenMP overhead dominates,
//       thermal effects cannot appear, and tau cannot be estimated precisely
//       enough to fix k = floor(T/tau). BBOB applies affine/rotation transforms,
//       making evaluation O(d^2) and putting real work inside the parallel
//       region.
//
// ============================ THREAD SAFETY ============================
// ioh::problem objects are STATEFUL: each holds evaluation counters, best-so-far
// state, and transformation buffers, mutated on every call to operator().
// Calling ONE instance concurrently from several OpenMP threads is a data race.
// It will often appear to work and will silently corrupt FE accounting —
// the worst possible failure mode for this study, since FE traces are what we
// report alongside time traces.
//
// Solution: ONE PROBLEM INSTANCE PER THREAD. Instances are constructed with the
// same (function id, instance id, dimension), so they are numerically identical;
// only their internal counters are separate. FE counting is done by the harness,
// not by the ioh objects, so per-thread instances do not fragment the count.
//
// Instances are built ONCE and reused, not constructed inside the timed region:
// BBOB construction generates rotation matrices (O(d^2) or worse) and would
// otherwise be charged to the algorithm's wall-clock budget.
// =======================================================================

#pragma once
#include <memory>
#include <vector>
#include <string>
#include <stdexcept>
// Include ONLY the BBOB problem headers, not the umbrella "ioh.hpp".
// The umbrella header pulls in the CEC suites, whose problem classes register
// themselves at static-initialisation time and load external .dat resource
// files; if IOH_RESOURCE_DIR is unset the process aborts before main() does
// anything useful. Restricting includes keeps the build self-contained.
#include <iostream>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include "ioh/problem/bbob.hpp"
#include "benchmarks.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace tfb {

// Direct construction rather than ProblemRegistry::create().
// The registry instantiates every registered problem family, including the CEC
// suites, which load external .dat resource files at static-init time and abort
// if IOH's resource root is not configured. Direct construction of the BBOB
// classes avoids that dependency entirely and keeps the build self-contained.
inline std::shared_ptr<ioh::problem::BBOB> make_bbob(int fid, int iid, int dim) {
    using namespace ioh::problem::bbob;
    switch (fid) {
        case 1:  return std::make_shared<Sphere<ioh::problem::BBOB>>(iid, dim);
        case 2:  return std::make_shared<Ellipsoid<ioh::problem::BBOB>>(iid, dim);
        case 3:  return std::make_shared<Rastrigin<ioh::problem::BBOB>>(iid, dim);
        case 6:  return std::make_shared<AttractiveSector<ioh::problem::BBOB>>(iid, dim);
        case 7:  return std::make_shared<StepEllipsoid<ioh::problem::BBOB>>(iid, dim);
        case 8:  return std::make_shared<Rosenbrock<ioh::problem::BBOB>>(iid, dim);
        case 9:  return std::make_shared<RosenbrockRotated<ioh::problem::BBOB>>(iid, dim);
        case 10: return std::make_shared<EllipsoidRotated<ioh::problem::BBOB>>(iid, dim);
        case 11: return std::make_shared<Discus<ioh::problem::BBOB>>(iid, dim);
        case 12: return std::make_shared<BentCigar<ioh::problem::BBOB>>(iid, dim);
        case 13: return std::make_shared<SharpRidge<ioh::problem::BBOB>>(iid, dim);
        case 14: return std::make_shared<DifferentPowers<ioh::problem::BBOB>>(iid, dim);
        case 15: return std::make_shared<RastriginRotated<ioh::problem::BBOB>>(iid, dim);
        case 16: return std::make_shared<Weierstrass<ioh::problem::BBOB>>(iid, dim);
        case 17: return std::make_shared<Schaffers10<ioh::problem::BBOB>>(iid, dim);
        case 18: return std::make_shared<Schaffers1000<ioh::problem::BBOB>>(iid, dim);
        case 19: return std::make_shared<GriewankRosenbrock<ioh::problem::BBOB>>(iid, dim);
        case 20: return std::make_shared<Schwefel<ioh::problem::BBOB>>(iid, dim);
        case 21: return std::make_shared<Gallagher101<ioh::problem::BBOB>>(iid, dim);
        case 22: return std::make_shared<Gallagher21<ioh::problem::BBOB>>(iid, dim);
        case 23: return std::make_shared<Katsuura<ioh::problem::BBOB>>(iid, dim);
        case 24: return std::make_shared<LunacekBiRastrigin<ioh::problem::BBOB>>(iid, dim);
        default: throw std::runtime_error("BBOB: unsupported fid " + std::to_string(fid));
    }
}

// Holds one ioh problem per thread plus a scratch input buffer per thread
// (ioh takes std::vector<double>; our PSO stores a flat array, so a copy is
// needed. The buffer is preallocated to keep allocation out of the hot loop).
class BBOBProblem {
public:
    BBOBProblem(int fid, int instance, int dim, int max_threads)
        : fid_(fid), instance_(instance), dim_(dim), nthreads_(max_threads)
    {
        probs_.reserve(nthreads_);
        bufs_.assign(nthreads_, std::vector<double>(dim_));
        for (int t = 0; t < nthreads_; ++t)
            probs_.push_back(make_bbob(fid_, instance_, dim_));
        name_    = probs_[0]->meta_data().name;
        lower_   = probs_[0]->bounds().lb[0];
        upper_   = probs_[0]->bounds().ub[0];
        optimum_ = probs_[0]->optimum().y;
    }

    // Thread-safe: dispatches to this thread's own instance.
    double eval(const double* x, int n) {
#ifdef _OPENMP
        const int t = omp_get_thread_num();
#else
        const int t = 0;
#endif
        // Defensive: if OMP hands us a thread id beyond what we allocated,
        // fail loudly rather than race. Silent corruption is unacceptable here.
        if (t >= nthreads_)
            throw std::runtime_error("BBOB: thread id exceeds allocated instances; "
                                     "construct with max_threads >= omp_get_max_threads()");
        auto& b = bufs_[t];
        for (int i = 0; i < n; ++i) b[i] = x[i];
        return (*probs_[t])(b);
    }

    Problem as_problem() {
        Problem p;
        p.name    = name_;
        p.dim     = dim_;
        p.lower   = lower_;
        p.upper   = upper_;
        p.optimum = optimum_;
        p.f = [this](const double* x, int n) { return this->eval(x, n); };
        return p;
    }

    int fid()      const { return fid_; }
    int instance() const { return instance_; }
    const std::string& name() const { return name_; }

private:
    int fid_, instance_, dim_, nthreads_;
    std::string name_;
    double lower_ = -5.0, upper_ = 5.0, optimum_ = 0.0;
    std::vector<std::shared_ptr<ioh::problem::BBOB>> probs_;
    std::vector<std::vector<double>> bufs_;
};

// A landscape-spanning subset of the 24 BBOB functions.
// Chosen to cover the five official groups rather than to favour any algorithm:
//   f1  Sphere              - separable, unimodal
//   f8  Rosenbrock          - low/moderate conditioning, non-separable
//   f10 Ellipsoid (rotated) - high conditioning, unimodal
//   f15 Rastrigin (rotated) - multimodal, adequate global structure
//   f21 Gallagher 101 peaks - multimodal, weak global structure
// Report this selection criterion explicitly in the paper (§3.3); a reviewer
// will ask why these five.
inline std::vector<int> default_fids() { return {1, 8, 10, 15, 21}; }

} // namespace tfb
