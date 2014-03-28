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

std::vector<double> processSamples(
        std::vector<double> samples
        , const std::string & instructions
        , const std::map<std::string,
            std::function<std::vector<double>(
                const std::vector<double> &,
                std::vector<double>)> > & functions
        );


std::map<std::string,
    std::function<std::vector<double>(
        const std::vector<double> &,
        std::vector<double>
    )>> createProcessingFunctions();

} // namespace dimf
