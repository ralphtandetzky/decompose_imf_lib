/// @file
///
/// @author Ralph Tandetzky
/// @date 24 Mar 2014

#pragma once

#include <complex>
#include <functional>
#include <memory>
#include <vector>

namespace dimf
{

/// Contains all the parameters needed for an imf decomposition.
struct OptimizationParams
{
    std::vector<double> samples{};
    size_t swarmSize{};
    double angleDevDegs{};
    double amplitudeDev{};
    double crossOverProb{};
    double diffWeight{};
    size_t nParams{};
    double initSigmaUnits{};
    double initTauUnits{};
    double nodeDevUnits{};
    double sigmaDevUnits{};
    double tauDevUnits{};
    std::string preprocessing{};
    std::string interprocessing{};
    double xIntervalWidth{};
    std::function<
            std::vector<std::complex<double> >(
                const std::vector<double> &)> initializer{};
    std::function<void(
            const std::vector<double> & bestParams
            , double cost
            , size_t nSamples
            , size_t nIter
            , const std::vector<double> & f )> receiveBestFit{};
    std::function<bool(size_t nIter)> shallCancel{};
};

template <typename F>
void iterateMembers( OptimizationParams & params, F && f )
{
    f( params.samples        , "samples"         );
    f( params.swarmSize      , "swarmSize"       );
    f( params.angleDevDegs   , "angleDevDegs"    );
    f( params.amplitudeDev   , "amplitudeDev"    );
    f( params.crossOverProb  , "crossOverProb"   );
    f( params.diffWeight     , "diffWeight"      );
    f( params.nParams        , "nParams"         );
    f( params.initSigmaUnits , "initSigmaUnits"  );
    f( params.initTauUnits   , "initTauUnits"    );
    f( params.nodeDevUnits   , "nodeDevUnits"    );
    f( params.sigmaDevUnits  , "sigmaDevUnits"   );
    f( params.tauDevUnits    , "tauDevUnits"     );
    f( params.preprocessing  , "preprocessing"   );
    f( params.interprocessing, "interprocessing" );
    f( params.xIntervalWidth , "xIntervalWidth"  );
    f( params.initializer    , "initializer"     );
    f( params.receiveBestFit , "receiveBestFit"  );
    f( params.shallCancel    , "shallCancel"     );
}

/// @brief An active object performing imf decomposition.
///
/// All methods run asynchroneously.
class OptimizationTask
{
public:
    OptimizationTask();
    ~OptimizationTask();

    void restart( OptimizationParams params );
    void cancel();
    void continueWithNextImf();

private:
    struct Impl;
    std::unique_ptr<Impl> m;
};

} // namespace dimf
