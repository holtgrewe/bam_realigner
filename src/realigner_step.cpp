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
#include <seqan/realign.h>
#include <seqan/store.h>
#include <seqan/misc/misc_interval_tree.h>

#include "bam_realigner_options.h"

namespace {  // anonymous namespace

}  // anonymous namespace


// ---------------------------------------------------------------------------
// Class RealignerStepImpl
// ---------------------------------------------------------------------------

class RealignerStepImpl
{
public:
    RealignerStepImpl(seqan::BamFileOut & bamFileOut,
                      seqan::VirtualStream<char, seqan::Output> & msasTxtOut,
                      seqan::BamFileIn & bamFileIn,
                      seqan::BamIndex<seqan::Bai> & baiIndex,
                      seqan::FaiIndex & faiIndex,
                      seqan::GenomicRegion const & region,
                      BamRealignerOptions const & options) :
            bamFileOut(bamFileOut), msasTxtOut(msasTxtOut), bamFileIn(bamFileIn), baiIndex(baiIndex),
            faiIndex(faiIndex), region(region), options(options)
    {
        extendRegion();
    }

    void run();

private:

    // Typedefs for using the store a big more comfortably.
    typedef seqan::FragmentStore<> TFragmentStore;
    typedef TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef seqan::Value<TAlignedReadStore>::Type TAlignedRead;
    typedef TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef seqan::Value<TReadSeqStore>::Type TReadSeq;
    typedef seqan::Gaps<TReadSeq, seqan::AnchorGaps<TAlignedRead::TGapAnchors> > TReadGaps;

    typedef TFragmentStore::TContigStore TContigStore;
    typedef seqan:: Value<TContigStore>::Type TContig;
    typedef TFragmentStore::TContigSeq TContigSeq;
    typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<TContig::TGapAnchors> > TContigGaps;

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
    // Build FragmentStore from aligned records.
    void buildFragmentStore();
    // Perform realignment on store.
    void performRealignment();
    // Update the BAM records from MSA stored in store.
    void updateBamRecords();
    // Write out BAM records.
    void writeBamRecords();

    // The reference sequence window.
    seqan::Dna5String ref;
    // The alignment records overlapping with the window.
    std::vector<seqan::BamAlignmentRecord> records;

    // Output files.
    seqan::BamFileOut & bamFileOut;
    seqan::VirtualStream<char, seqan::Output> & msasTxtOut;
    // Input BAM and FAI index.
    seqan::BamFileIn & bamFileIn;
    seqan::BamIndex<seqan::Bai> & baiIndex;
    seqan::FaiIndex & faiIndex;
    // The region to realign.
    seqan::GenomicRegion region;
    // The used FragmentStore.
    seqan::FragmentStore<> store;

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
    // Build FragmentStore from the aligned alignment records.
    buildFragmentStore();
    // Perform realignment.
    performRealignment();
    // Update the BAM records before writing out.
    updateBamRecords();
    // Write out BAM records.
    writeBamRecords();
}

// TODO(holtgrew): This function is much too big, split into smaller ones!

void RealignerStepImpl::buildFragmentStore()
{
    // Set sequence into store.
    resize(store.contigStore, 1);
    store.contigStore[0].seq = ref;
    resize(store.contigNameStore, 1);
    region.toString(store.contigNameStore[0]);

    // Stores (refPos, numInsertions) for each read, used for distributing gaps to other reads below.
    std::vector<std::map<int, int> > readInsertions(records.size());
    // Stores (refPos, numGaps) gaps to insert into the reference.
    std::map<int, int> refGaps;

    // TODO(holtgrew): The code below does NOT handle soft- and hard-clipping.

    // We append the reads ignoring pairing and forward/reverse information.
    for (auto const & record : records)
    {
        // -------------------------------------------------------------------
        // Append read's sequence and id information.
        // -------------------------------------------------------------------

        auto readID = appendRead(store, record.seq, record.qName);

        // -------------------------------------------------------------------
        // Append alignment for read if it is aligned in the BAM file.
        // -------------------------------------------------------------------

        if (hasFlagUnmapped(record))
            continue;
        int beginPos = record.beginPos - region.beginPos;
        int clippedLength = 0;
        _getClippedLength(record.cigar, clippedLength);
        int endPos = beginPos + clippedLength;
        auto alignmentID = appendAlignedRead(store, readID, 0, beginPos, endPos);

        // -------------------------------------------------------------------
        // Update read gaps and begin offset in case of leading gaps.
        // -------------------------------------------------------------------

        TReadGaps readGaps(store.readSeqStore[readID],
                           store.alignedReadStore[alignmentID].gaps);
        unsigned leadingGaps = cigarToGapAnchorRead(record.cigar, readGaps);
        store.alignedReadStore[alignmentID].beginPos += leadingGaps;
        store.alignedReadStore[alignmentID].beginPos += leadingGaps;

        // Update readGapPositions and refGapPositions.
        int refPos = store.alignedReadStore[alignmentID].beginPos;
        int readPos = 0;
        auto readGapsIt = begin(readGaps, seqan::Standard());
        if (options.verbosity >= 3)
            std::cerr << "READ\t" << store.readNameStore[readID] << "\n";
        for (auto cigar : record.cigar)
        {
            switch (cigar.operation)
            {
                case 'D':  // deletion from read => gap in read
                    // no need to insert gap, already registered in cigarToGapAnchorRead()
                    readGapsIt += cigar.count;
                    readPos += cigar.count;
                    if (options.verbosity >= 3)
                        std::cerr << "\t" << cigar.operation << "\tcigar.count=" << cigar.count
                                  << "\treadPos=" << readPos
                                  << "\trefPos=" << refPos
                                  << "\n";
                    break;

                case 'I':  // insertion into reference => gap in ref
                    refGaps[refPos] = std::max(refGaps[refPos], (int)cigar.count);
                    readInsertions[readID][refPos] = cigar.count;
                    readGapsIt += cigar.count;
                    readPos += cigar.count;
                    if (options.verbosity >= 3)
                        std::cerr << "\t" << cigar.operation << "\tcigar.count=" << cigar.count
                                  << "\treadPos=" << readPos
                                  << "\trefPos=" << refPos
                                  << "\n";
                    break;

                case 'X':  // aligned characters
                case '=':
                case 'M':
                    readGapsIt += cigar.count;
                    refPos += cigar.count;
                    readPos += cigar.count;
                    if (options.verbosity >= 3)
                        std::cerr << "\t" << cigar.operation << "\tcigar.count=" << cigar.count
                                  << "\treadPos=" << readPos
                                  << "\trefPos=" << refPos
                                  << "\n";
                    break;

                case 'P':  // ignore paddings in read
                default:
                    if (options.verbosity >= 3)
                        std::cerr << "\t" << cigar.operation << "\tcigar.count=" << cigar.count
                                  << "\treadPos=" << readPos
                                  << "\trefPos=" << refPos
                                  << "\n";
                    break;
            }
        }

        if (options.verbosity >= 3)
            std::cerr << "\t\t" << readGaps << "\n";
    }

    // -----------------------------------------------------------------------
    // Project individual insertions to MSA
    // -----------------------------------------------------------------------

    // Sort aligned reads, so we can lowerBound() below.
    sortAlignedReads(store.alignedReadStore, seqan::SortEndPos());
    sortAlignedReads(store.alignedReadStore, seqan::SortBeginPos());

    // Obtain contig gaps.
    TContigGaps contigGaps(store.contigStore[0].seq, store.contigStore[0].gaps);

    // Build interval tree for overlapping reads lookup (cargo is alignment id, not read id!)
    typedef seqan::IntervalAndCargo<int, int> TInterval;
    seqan::String<TInterval> intervals;
    for (auto const & el : store.alignedReadStore)
        appendValue(intervals, TInterval(toSourcePosition(contigGaps, el.beginPos),
                                         toSourcePosition(contigGaps, el.endPos),
                                         el.id));
    seqan::IntervalTree<int, int> tree(intervals);

    // Insert gaps into overlapping reads.
    seqan::String<int> results;
    for (auto it = refGaps.rbegin(); it != refGaps.rend(); ++it)
    {
        // Update overlapping reads.
        clear(results);
        findIntervals(tree, it->first, it->first + 1, results);
        for (auto alignmentID : results)
        {
            auto & el = store.alignedReadStore[alignmentID];
            int viewPos = it->first - el.beginPos;

            TReadGaps readGaps(store.readSeqStore[el.readId],
                               store.alignedReadStore[alignmentID].gaps);

            // Leading gaps are handled below in shifting step and trailing gaps are ignored.
            if (viewPos == 0 || viewPos == (int)(el.endPos - el.beginPos))
                continue;

            int delta = readInsertions[el.readId].count(it->first) ? readInsertions[el.readId][it->first] : 0;
            insertGaps(readGaps, viewPos, it->second - delta);
            el.endPos += it->second - delta;
            if (options.verbosity >= 3)
                std::cerr << "INSERTING READ GAPS\t" << el.readId << "\t" << readGaps << "\t" << viewPos
                          << "\t" << (it->second - delta) << "\n";
        }
        // Update reads left of where we are inserting gaps.
        auto itA = lowerBoundAlignedReads(store.alignedReadStore, toViewPosition(contigGaps, it->first), seqan::SortBeginPos());
        for (; itA != end(store.alignedReadStore, seqan::Standard()); ++itA)
        {
            if (options.verbosity >= 3)
                std::cerr << "SHIFTING LEFT\t" << itA->readId << "\t" << store.readSeqStore[itA->readId]
                          << " by " << it->second << "\n";
            itA->beginPos += it->second;
            itA->endPos += it->second;
        }
    }
    
    // Insert gaps into contig.
    for (auto it = refGaps.rbegin(); it != refGaps.rend(); ++it)
    {
        // Insert gaps into contigs.
        if (options.verbosity >= 3)
            std::cerr << "INSERTING CONTIG GAPS\t" << it->first << "\t" << it->second << "\n";
        insertGaps(contigGaps, it->first, it->second);
    }

    // Print store after loading.
    if (options.verbosity >= 2 || msasTxtOut.good())
    {
        if (msasTxtOut.good())
            msasTxtOut << ">" << store.contigNameStore[0] << " before realignment\n";
        if (options.verbosity >= 2)
            std::cerr << ">" << store.contigNameStore[0] << " before realignment\n";

        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        if (msasTxtOut.good())
            printAlignment(msasTxtOut, layout, store, 0, 0, (int)(region.endPos - region.beginPos), 0, 10000);
        if (options.verbosity >= 2)
            printAlignment(std::cerr, layout, store, 0, 0, (int)(region.endPos - region.beginPos), 0, 10000);
    }
}

void RealignerStepImpl::performRealignment()
{
    if (options.verbosity >= 2)
        std::cerr << "Performing realignment\n";
    reAlignment(store, 0, 1, 10, 1, 0, 0, /*debug=*/(options.verbosity >= 3),
                /*printTiming=*/(options.verbosity >= 2));
    if (options.verbosity >= 2)
        std::cerr << "  => DONE\n";

    // Print store after realignment.
    if (options.verbosity >= 2 || msasTxtOut.good())
    {
        if (msasTxtOut.good())
            msasTxtOut << ">" << store.contigNameStore[0] << " before realignment\n";
        if (options.verbosity >= 2)
            std::cerr << ">" << store.contigNameStore[0] << " after realignment\n";

        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        int minPos = seqan::maxValue<int>(), maxPos = seqan::minValue<int>();
        for (auto const & el : store.alignedReadStore)
        {
            minPos = std::min(minPos, (int)el.beginPos);
            maxPos = std::max(maxPos, (int)el.endPos);
        }
        if (minPos == seqan::maxValue<int>())
            minPos = maxPos = 0;
        if (msasTxtOut.good())
            printAlignment(msasTxtOut, layout, store, 0, minPos, maxPos, 0, 10000);
        if (options.verbosity >= 2)
            printAlignment(std::cerr, layout, store, 0, minPos, maxPos, 0, 10000);
    }
}

void RealignerStepImpl::updateBamRecords()
{
    // Obtain contig gaps.
    TContigGaps contigGaps(store.contigStore[0].seq, store.contigStore[0].gaps);

    for (auto const & el : prefix(store.alignedReadStore, length(store.alignedReadStore) - 1))
    {
        auto & record = records[el.readId];

        // Obtain read gaps and clipped contig gaps.
        TReadGaps readGaps(store.readSeqStore[el.readId], el.gaps);
        TContigGaps clippedContigGaps(store.contigStore[0].seq,
                                      store.contigStore[0].gaps);
        setClippedEndPosition(contigGaps, el.endPos);
        setClippedBeginPosition(contigGaps, el.beginPos);

        // Update alignment position and alignment info.
        record.beginPos = region.beginPos + toSourcePosition(contigGaps, el.beginPos);
        getCigarString(record.cigar, clippedContigGaps, readGaps);
    }
}

void RealignerStepImpl::writeBamRecords()
{
    for (auto const & record : records)
        writeRecord(bamFileOut, record);
}

// ---------------------------------------------------------------------------
// Class RealignerStep
// ---------------------------------------------------------------------------

RealignerStep::RealignerStep(seqan::BamFileOut & bamFileOut,
                             seqan::VirtualStream<char, seqan::Output> & msaTxtOut,
                             seqan::BamFileIn & bamFileIn,
                             seqan::BamIndex<seqan::Bai> & baiIndex,
                             seqan::FaiIndex & faiIndex,
                             seqan::GenomicRegion const & region,
                             BamRealignerOptions const & options) :
        impl(new RealignerStepImpl(bamFileOut, msaTxtOut, bamFileIn, baiIndex, faiIndex,
                                   region, options))
{}

RealignerStep::~RealignerStep()
{}

void RealignerStep::run()
{
    impl->run();
}

