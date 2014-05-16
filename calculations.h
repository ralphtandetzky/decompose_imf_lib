#pragma once

#include <complex>
#include <vector>

namespace dimf
{

std::vector<std::complex<double>> groupPairsToComplex(
        const std::vector<double> & pairsOfReals );

std::vector<std::complex<double>> calculateSigmaFunctionFromSigmaSequence(
        std::vector<std::complex<double>> sigma );

std::vector<double> calculateImfFromSigmaFunction(
        const std::vector<std::complex<double>> & sigma );

std::vector<double> calculateImfFromPairsOfReals(
        const std::vector<double> & input );

double sumOfSquaresOfDifference(
        const std::vector<double> & lhs,
        const std::vector<double> & rhs );

double costFunction( const std::vector<double> & f,
                     const std::vector<double> & pairsOfReals,
                     double frequencySwingFactor );

std::vector<std::complex<double>>
    derive( std::vector<std::complex<double>> f );

double boundaryCondition(std::vector<std::complex<double>> sigma_seq
                         , double frequencySwingFactor);

std::vector<std::complex<double>>
    getInitialApproximationByInterpolatingZeros(
        const std::vector<double> & f );

/// @l \sum_{i}p_{2i}e^{-\frac{(z-p_{2i+1})^2}{2\sigma^2}}
std::vector<double>
    getSamplesFromRadialBase(
        const std::vector<double> & params, double sigma, size_t nSamples );

/// @l \sum_{i}\frac{p_{2i}}{1+e^{-\frac{z-p_{2i+1}}{\tau}}}
std::vector<double>
    getSamplesFromLogisticFunctionBase(
        const std::vector<double> & params, double tau, size_t nSamples );

std::vector<double> getSamplesFromParams(
        std::vector<double> v , size_t nSamples );

} // namespace dimf
