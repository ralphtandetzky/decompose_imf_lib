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

/// @brief Maps @c std::string to sample processing functors.
///
/// Additionally to the input samples the functors take a vector
/// of @c doubles as an argument.
using ProcessingFunctions = std::map<std::string,
    std::function<std::vector<double>(
        const std::vector<double> & params,
        std::vector<double> samples)> >;

/// @brief Parses instructions and applies them to a vector of samples.
///
/// @param instructions This string will be devided into lines. Each line
/// shall contain exactly one instruction. Empty lines will be ignored.
/// The instructions will be applied in the order of appearance. The syntax
/// of the instructions should be
///   @code
///     instruction_name [param_val_1] [param_val_2] ...
///   @endcode
/// where @c instruction_name is a chacter string without whitespaces which
/// must be a valid key in the @c functions map. The parameters (if any)
/// must be floating point values. The @c instruction_name and the parameters
/// must be separated by whitespaces but not line breaks. Subsequent
/// whitespaces will be merged into one. The parameters will be passed as a
/// vector of @c doubles to the processing function together with the vector
/// of samples. Which number of parameters is legal depends upon the
/// instruction.
std::vector<double> processSamples(
        std::vector<double> samples
        , const std::string & instructions
        , const ProcessingFunctions & functions
        );

/// @brief Constructs functors for processing of samples.
///
/// These functors can filter, clip, add noise and other things to
/// sampled data during preprocessing, interprocessing and
/// postprocessing of data. The returned map serves as input for
/// @c processSamples() which parses a string of instructions and
/// applies the instructions to a vector of samples.
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
