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
            std::vector<std::complex<double>>(
                const std::vector<double> &)> initializer{};
    std::function<void(
            const std::vector<double> & bestParams
            , double cost
            , size_t nSamples
            , size_t nIter
            , const std::vector<double> & f )> receiveBestFit{};
};

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
