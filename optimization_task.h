/// @file
///
/// @author Ralph Tandetzky
/// @date 24 Mar 2014

#pragma once

namespace dimf
{

/// Contains all the parameters needed for an imf decomposition.
struct OptimizationParams
{
    ///
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
