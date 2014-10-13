// ==========================================================================
//                               BAM Realigner
// ==========================================================================
// Copyright (c) 2014, Manuel Holtgrewe
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef BAM_REALIGNER_SRC_BAM_REALIGNER_OPTIONS_H_
#define BAM_REALIGNER_SRC_BAM_REALIGNER_OPTIONS_H_

#include <iosfwd>
#include <string>
#include <stdexcept>

// ----------------------------------------------------------------------------
// Class InvalidArgumentsException
// ----------------------------------------------------------------------------

class InvalidCommandLineArgumentsException : public std::runtime_error
{
public:
    InvalidCommandLineArgumentsException() : std::runtime_error("")
    {}
};

// ----------------------------------------------------------------------------
// Class BamRealignerOptions
// ----------------------------------------------------------------------------

class BamRealignerOptions
{
public:
    // Verbosity: 0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Indexed input alignment file (.bam).
    std::string inAlignmentPath;
    // Input reference (.fasta), will build FAI index for it.
    std::string inReferencePath;
    // Input intervals file (.bed)
    std::string inIntervalsPath;
    // Output BAM file.
    std::string outAlignmentPath;
    // Output text file with MSAs.
    std::string outMsasPath;

    // Additional radius around target intervals to extract reads from.
    int windowRadius;

    BamRealignerOptions() : verbosity(1), windowRadius(100)
    {}

    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

BamRealignerOptions parseCommandLine(int argc, char ** argv);

#endif  // #ifndef BAM_REALIGNER_SRC_BAM_REALIGNER_OPTIONS_H_
