/// @file
///
/// @author Ralph Tandetzky
/// @date 17 Feb 2014

#pragma once

#include <map>
#include <functional>
#include <string>
#include <vector>

namespace dimf
{

using ProcessingFunctions = std::map<std::string,
    std::function<std::vector<double>(
        const std::vector<double> &,
        std::vector<double>)> >;

std::vector<double> processSamples(
        std::vector<double> samples
        , const std::string & instructions
        , const ProcessingFunctions & functions
        );

ProcessingFunctions createProcessingFunctions();

/// @brief Convenience function. Behaves like
///   @code
///     processSamples( samples, instructions, createProcessingFunctions() );
///   @endcode
std::vector<double> processSamples(
        std::vector<double> samples
        , const std::string & instructions
        );

} // namespace dimf
