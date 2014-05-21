#include "decompose_imf_lib/optimization_task.h"

#include "decompose_imf_lib/calculations.h"
#include "decompose_imf_lib/processing.h"

#include "cpp_utils/math_constants.h"
#include "cpp_utils/optimize.h"
#include "cpp_utils/parallel_executor.h"
#include "cpp_utils/std_make_unique.h"

#include "qt_utils/exception_handling.h"

#include <opencv/cv.h>

namespace dimf
{

//void OptimizationTask::start( OptimizationParams params )
std::vector<std::vector<double> > runOptimization(
        const OptimizationParams & params,
        std::ostream & os )
{
    using cu::pi;

    auto done = false;
    auto currentImf = size_t{0};

    const auto processingFunctions = createProcessingFunctions();
    const auto samples = getPreprocessedSamples( params );
    auto imfs = std::vector<std::vector<double> >{};
    auto bestParamSets = std::vector<std::vector<double> >{};

    // This variable is shared between 'shallTerminate' and 'sendBestFit'.
    auto nIter = 0;

    while ( !done )
    {
        const auto oldImf = currentImf;
        os << "Optimizing IMF " << currentImf << "." << std::endl;
        auto f = samples;
        f = processSamples( f,
                params.interprocessing,
                processingFunctions );
        for ( auto i = size_t(); i < imfs.size(); ++i )
        {
            if ( i == currentImf )
                continue;
            const auto & imf = imfs[i];
            cu::subAssign( begin(f), end(f), 1., begin(imf), end(imf) );
        }
        const auto nSamples = f.size();

        // calculate an initial approximation and swarm
        CU_ASSERT_THROW( params.initializer,
                         "The algorithm for the initial "
                         "approximation has not been specified." );
        const auto initApprox = params.initializer(f);

        // calculate base with equidistant center points of logistic functions
        auto logisticBase = std::vector<std::vector<double>>{};
        auto nodes = std::vector<double>{};
        const auto factor = nSamples / params.xIntervalWidth;
        const auto initSigma = params.initSigmaUnits * factor;
        for ( auto i = size_t{0}; i < params.nParams; ++i )
        {
            nodes.push_back( (i+.5)*initApprox.size()/params.nParams );
            logisticBase.push_back( getSamplesFromLogisticFunctionBase(
                { 1., nodes.back() }, initSigma, initApprox.size() ) );
        }
        auto logisticBaseMat = cv::Mat(
                    logisticBase.front().size(), logisticBase.size(),
                    CV_64FC1, cv::Scalar::all(0) );
        for ( auto row = 0; row < logisticBaseMat.rows; ++row )
        {
            for ( auto col = 0; col < logisticBaseMat.cols; ++col )
            {
                logisticBaseMat.at<double>(row,col) =
                        logisticBase.at(col).at(row);
            }
        }

        // calculate the element closest to initApprox that lies in the
        // space spanned by the base elements. Also calculate the
        // coefficients of that minimum element with respect to the given
        // base
        auto initApproxMat = cv::Mat(
                    initApprox.size(), 1, CV_64FC1, cv::Scalar::all(0) );
        for ( auto row = 0; row < initApproxMat.rows; ++row )
        {
            initApproxMat.at<double>( row ) = initApprox.at(row).imag();
        }
        const auto invLogisticBase =
                cv::Mat{ logisticBaseMat.inv( cv::DECOMP_SVD ) };
        const auto bestApproxMat =
                cv::Mat{ invLogisticBase * initApproxMat };

        auto swarm = std::vector<std::vector<double> >( params.swarmSize );
        auto rng = std::mt19937{};
        {
            auto normal_dist = std::normal_distribution<>{};
            auto uniform = std::uniform_real_distribution<double>{-1,1};
            for ( auto & x : swarm )
            {
                const auto sigmaDev  = params.sigmaDevUnits  * factor;
                const auto tauDev    = params.tauDevUnits    * factor;
                const auto initTau   = params.initTauUnits   * factor;
                const auto nodeDev   = params.nodeDevUnits   * factor;
                const auto angleDev  = params.angleDevDegs/180*cu::pi;
                for ( auto i = size_t{0}; i < nodes.size(); ++i )
                {
                    x.push_back( params.amplitudeDev*normal_dist(rng) );
                    x.push_back( nodes[i] + nodeDev*normal_dist(rng) );
                }
                for ( auto i = size_t{0}; i < nodes.size(); ++i )
                {
                    x.push_back( bestApproxMat.at<double>( i ) +
                                 angleDev*uniform(rng) );
                    x.push_back( nodes[i] );
                }
                x.push_back( initSigma + sigmaDev*normal_dist(rng) );
                x.push_back( initTau + tauDev*normal_dist(rng) );
            }
            if ( currentImf < bestParamSets.size() )
                swarm.front() = bestParamSets.at(currentImf);
        }

        // cost function for optimization
        const auto cost = [&f, nSamples, &params]( std::vector<double> v ) -> double
        {
            return costFunction( f,
                getSamplesFromParams( std::move(v), nSamples ),
                params.freqSwingFactor );
        };

        // function which returns whether the
        // optimization algorithm shall terminate.
        const auto shallTerminate = [&]( const decltype(swarm) & )-> bool
        {
            ++nIter;
            const auto nextImf = params.howToContinue( nIter );
            if ( nextImf == ~size_t{0} )
            {
                done = true;
                return true;
            }
            if ( nextImf == currentImf )
                return false;
            currentImf = nextImf;
            return true;
        };

        // function which is called by the optimization algorithm
        // each time the best fit is improved.
        auto bestParams = std::vector<double>{};
        const auto sendBestFit = [&](
                const std::vector<double> & v_, double cost )
        {
            bestParams = v_;
            params.receiveBestFit(bestParams,cost,nSamples,nIter,f);
        };

        // perform the optimization.
        swarm = cu::differentialEvolution(
            std::move(swarm), params.crossOverProb, params.diffWeight,
            cost, shallTerminate, sendBestFit, rng );

        const auto bestImf = calculateImfFromPairsOfReals(
                    getSamplesFromParams( bestParams, nSamples ) );
        if ( oldImf >= imfs.size() )
        {
            CU_ASSERT_THROW( oldImf == imfs.size(),
                "Invalid IMF index. Have some indexes been skipped?" );
            imfs.push_back({});
            bestParamSets.push_back({});
        }
        imfs[oldImf] = bestImf;
        bestParamSets[oldImf] = bestParams;
    } // while loop
    return imfs;
}


std::vector<double> getPreprocessedSamples(
        const OptimizationParams & params )
{
    const auto processingFunctions = createProcessingFunctions();
    return processSamples(
        params.samples,
        params.preprocessing,
        processingFunctions );
}


} // namespace dimf
