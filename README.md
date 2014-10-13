BAM Realigner
=============

SeqAn-Based Realignment of BAM files.

Building
--------

Using
-----

    # bam_realigner [-v] --in-alignment ALI.bam --in-reference REF.fa \
                         --in-intervals REGIONS.intervals

Caveats
-------

* Soft clippings are not interpreted (yet).
* The program only writes out records overlapping with the target regions at
  the moment.
* Target regions are given with a window radius of 100bp (user parameter) and
  if two such regions overlap then the records are written out twice.
* When the position of a read changes then its mate's PNEXT field is not updated.
* Unaligned reads are not written out.