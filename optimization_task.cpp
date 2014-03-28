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

struct OptimizationTask::Impl
{
    /////////////////////////////////
    // Data shared between threads //
    /////////////////////////////////
    struct SharedData
    {
        // Holds whether the currently running optimization task (if any)
        // has been cancelled.
        bool cancelled = true;
        // Holds whether the currently running optimization task (if any)
        // should proceed with calculating the imf of the current residue
        // function.
        bool shall_calculate_next_imf = false;
    };
    cu::Monitor<SharedData> shared;

    // single-threaded for background ops.
    cu::ParallelExecutor worker;
};


OptimizationTask::OptimizationTask()
    : m( std::make_unique<Impl>() )
{
}


OptimizationTask::~OptimizationTask()
{
    cancel();
}


void OptimizationTask::restart( OptimizationParams params )
{
    // cancel the currently running operation (if any)
    cancel();

    // concurrently launch optimization on the worker thread
    m->worker.addTask( [=]()
    {
    QU_HANDLE_ALL_EXCEPTIONS_FROM
    {
        auto f = params.samples;
        using cu::pi;
        m->shared( []( Impl::SharedData & shared )
        {
            shared.cancelled = false;
            shared.shall_calculate_next_imf = false;
        } );

        auto done = false;

        const auto processingFunctions = createProcessingFunctions();
        auto residue = processSamples( f, params.preprocessing, processingFunctions );

        while ( !done )
        {
            f = processSamples( residue, params.interprocessing, processingFunctions );
            const auto nSamples = f.size();

            // calculate an initial approximation and swarm
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

            auto swarm = std::vector<std::vector<double>>( params.swarmSize );
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
            }

            // cost function for optimization
            const auto cost = [&f, nSamples]( std::vector<double> v ) -> double
            {
                return costFunction( f,
                    getSamplesFromParams( std::move(v), nSamples ) );
            };

            // This variable is shared between 'shallTerminate' and 'sendBestFit'.
            auto nIter = 0;

            // function which returns whether the
            // optimization algorithm shall terminate.
            const auto shallTerminate = [&]( const decltype(swarm) & )-> bool
            {
                ++nIter;
                return m->shared( [&]( Impl::SharedData & shared ) -> bool
                {
                    if ( shared.shall_calculate_next_imf )
                    {
                        shared.shall_calculate_next_imf = false;
                        nIter = 0;
                        return true;
                    }
                    if ( shared.cancelled )
                    {
                        done = true;
                        return true;
                    }
                    return false;
                } );
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
            cu::subAssign( begin(residue), end(residue), 1.,
                           begin(bestImf), end(bestImf) );
        } // while loop
    }; // QU_HANDLE_ALL_EXCEPTIONS_FROM
    } ); // task for worker
}


void OptimizationTask::cancel()
{
    m->shared( []( Impl::SharedData & shared )
    {
        shared.cancelled = true;
    } );
}


void OptimizationTask::continueWithNextImf()
{
    m->shared( []( Impl::SharedData & shared )
    {
        shared.shall_calculate_next_imf = true;
    } );
}

} // namespace dimf
