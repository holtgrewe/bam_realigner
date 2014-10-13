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

#include "bam_realigner_options.h"

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/intervals_io.h>

// ----------------------------------------------------------------------------
// Class BamRealignerOptions
// ----------------------------------------------------------------------------

void BamRealignerOptions::print(std::ostream & out) const
{
    out << "__OPTIONS________________________________________________________\n"
        << "\n"
        << "VERBOSITY      \t" << verbosity << "\n"
        << "\n"
        << "INPUT REFERENCE\t" << inReferencePath << "\n"
        << "INPUT ALIGNMENT\t" << inAlignmentPath << "\n"
        << "INPUT INTERVALS\t" << inIntervalsPath << "\n";
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

BamRealignerOptions parseCommandLine(int argc, char ** argv)
{
    BamRealignerOptions result;

    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bam_realigner");

    // Set short description, version, and date.
    setShortDescription(parser, "BAM realigner");
    setVersion(parser, "0.1");
    setDate(parser, "October 2014");

    // Define usage line and long description.
    addUsageLine(parser, "--in-alignment ALI.bam --in-reference REF.fa --in-intervals INT.bed");
    addDescription(parser, "Read realignments from BAM files.");

    addOption(parser, seqan::ArgParseOption("q",  "quiet",        "Quiet output"));
    addOption(parser, seqan::ArgParseOption("v",  "verbose",      "Verbose output"));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose output"));

    // Define Options -- Section Input / Output Optiosn
    addSection(parser, "Input / Output Options");

    addOption(parser, seqan::ArgParseOption("", "in-alignment", "Input alignment file.",
                                            seqan::ArgParseArgument::INPUT_FILE, "BAM"));
    setRequired(parser, "in-alignment", true);
    setValidValues(parser, "in-alignment", seqan::BamFileIn::getFileFormatExtensions());

    addOption(parser, seqan::ArgParseOption("", "in-reference", "Input reference file.",
                                            seqan::ArgParseArgument::INPUT_FILE, "FASTA"));
    setRequired(parser, "in-reference", true);
    setValidValues(parser, "in-reference", "fa fasta");

    addOption(parser, seqan::ArgParseOption("", "in-intervals", "Input bed file.",
                                            seqan::ArgParseArgument::INPUT_FILE, "FASTA"));
    setRequired(parser, "in-intervals", true);
    setValidValues(parser, "in-intervals", seqan::IntervalsFileIn::getFileFormatExtensions());

    // Define Options -- Algorithm Parameters
    addSection(parser, "Algorithm Parameters");

    addOption(parser, seqan::ArgParseOption("", "window-radius", "Window radius to use for extension.",
                                            seqan::ArgParseArgument::INTEGER, "LEN"));
    setDefaultValue(parser, "window-radius", 100);

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        throw InvalidCommandLineArgumentsException();

    // Extract option values.
    result.verbosity = isSet(parser, "quiet") ? 0 : result.verbosity;
    result.verbosity = isSet(parser, "verbose") ? 2 : result.verbosity;
    result.verbosity = isSet(parser, "very-verbose") ? 3 : result.verbosity;

    getOptionValue(result.inAlignmentPath, parser, "in-alignment");
    getOptionValue(result.inReferencePath, parser, "in-reference");
    getOptionValue(result.inIntervalsPath, parser, "in-intervals");

    getOptionValue(result.windowRadius, parser, "window-radius");

    return result;
}
