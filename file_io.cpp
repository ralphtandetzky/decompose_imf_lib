#include "file_io.h"

#include "../cpp_utils/exception.h"
#include <fstream>
#include <iterator>

namespace dimf
{

std::vector<double> readSamplesFromFile( const std::string & fileName )
{
    std::ifstream file{ fileName };
    if ( !file )
        CU_THROW( "Could not open the file \"" + fileName + "\"." );
    auto result = std::vector<double>( std::istream_iterator<double>(file),
                                       std::istream_iterator<double>() );
    if ( file.bad() )
        CU_THROW( "The file \"" + fileName +
                  "\" could not be read." );
    if ( result.empty() )
        CU_THROW( "The file \"" + fileName +
                  "\" does not contain samples." );
    if ( !file.eof() )
        CU_THROW( "The end of the file \"" + fileName +
                  "\" has not been reached." );

    return result;
}

} // namespace dimf
