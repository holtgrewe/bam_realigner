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

#include "realigner_step.h"

#include <iostream>
#include <vector>
#include <string>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/intervals_io.h>

#include "bam_realigner_options.h"

namespace {  // anonymous namespace

}  // anonymous namespace


// ---------------------------------------------------------------------------
// Class RealignerStepImpl
// ---------------------------------------------------------------------------

class RealignerStepImpl
{
public:
    RealignerStepImpl(seqan::BamFileIn & bamFileIn,
                      seqan::BamIndex<seqan::Bai> & baiIndex,
                      seqan::FaiIndex & faiIndex,
                      seqan::GenomicRegion const & region,
                      BamRealignerOptions const & options) :
            bamFileIn(bamFileIn), baiIndex(baiIndex), faiIndex(faiIndex), region(region), options(options)
    {
        extendRegion();
    }

    void run();

private:

    // Extend region by options.windowRadius.
    void extendRegion()
    {
        region.beginPos -= options.windowRadius;
        region.beginPos = std::max(0u, region.beginPos);
        region.endPos += options.windowRadius;
    }
    // Extend region by extents of alignment.
    void extendRegion(seqan::BamAlignmentRecord const & record)
    {
        if (record.rID != (int)region.rID)
            return;  // do not update if on difference contig
        region.beginPos = std::min((int)region.beginPos, record.beginPos);
        region.endPos = std::max((int)region.endPos, (int)(record.beginPos + getAlignmentLengthInRef(record)));
    }

    // Load reference sequence.
    void loadReference();
    // Load alignments;
    void loadAlignments();

    // The reference sequence window.
    seqan::Dna5String ref;
    // The alignment records overlapping with the window.
    std::vector<seqan::BamAlignmentRecord> records;

    // Input BAM and FAI index.
    seqan::BamFileIn & bamFileIn;
    seqan::BamIndex<seqan::Bai> & baiIndex;
    seqan::FaiIndex & faiIndex;
    // The region to realign.
    seqan::GenomicRegion region;

    // Options.
    BamRealignerOptions const & options;
};

void RealignerStepImpl::loadReference()
{
    if (options.verbosity >= 2)
        std::cerr << "Loading reference...\n";
    readRegion(ref, faiIndex, region);
    if (options.verbosity >= 2)
        std::cerr << "  => DONE\n";
}

void RealignerStepImpl::loadAlignments()
{
    if (options.verbosity >= 2)
        std::cerr << "Loading alignments...\n";

    // Translate region reference name to reference ID in BAM file.
    if (!getIdByName(region.rID, nameStoreCache(context(bamFileIn)), region.seqName))
    {
        std::string msg = std::string("Unknown reference ")  + toCString(region.seqName);
        throw seqan::IOError(msg.c_str());
    }

    // Jump to region using BAI file.
    bool hasAlignments = false;
    if (!jumpToRegion(bamFileIn, hasAlignments, region.rID, region.beginPos, region.endPos, baiIndex))
        throw seqan::IOError("Problem jumping in file.\n");
    if (!hasAlignments)
    {
        // Handle the case of no alignments in region.
        seqan::CharString buffer;
        region.toString(buffer);
        if (options.verbosity >= 1)
            std::cerr << "\nWARNING: No alignments in region " << buffer << "\n";
        return;
    }

    // Load alignments.
    seqan::BamAlignmentRecord record;
    while (true)
    {
        readRecord(record, bamFileIn);
        if (record.rID == seqan::BamAlignmentRecord::INVALID_REFID)
            break;  // done, no more aligned records
        if (std::make_pair(record.rID, record.beginPos) > std::make_pair((int)region.rID, (int)region.endPos))
            break;  // done, no more records
        extendRegion(record);
        records.push_back(record);
    }

    if (options.verbosity >= 2)
        std::cerr << "  => DONE\n";
}

void RealignerStepImpl::run()
{
    // Load alignments, updates positions in region.
    loadAlignments();
    // Load reference sequence in regions.
    loadReference();
}

// ---------------------------------------------------------------------------
// Class RealignerStep
// ---------------------------------------------------------------------------

RealignerStep::RealignerStep(seqan::BamFileIn & bamFileIn,
                             seqan::BamIndex<seqan::Bai> & baiIndex,
                             seqan::FaiIndex & faiIndex,
                             seqan::GenomicRegion const & region,
                             BamRealignerOptions const & options) :
        impl(new RealignerStepImpl(bamFileIn, baiIndex, faiIndex, region, options))
{}

RealignerStep::~RealignerStep()
{}

void RealignerStep::run()
{
    impl->run();
}

