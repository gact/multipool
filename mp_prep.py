#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
mp_prep.py - Prepare a SNP allele depth file for Multipool inference.

This script takes as input one or more VCF files, and outputs a SNP allele
depth file that can be used as input in the Multipool inference script. For
a given pool sample, SNP alleles for that sample are matched to those of the
specified founder samples, and the depths of those alleles are output in the
same order as the founders are specified. All input variants are assumed to
have been quality-checked and appropriately filtered.

For more information on the Multipool package, see the original paper by
Edwards and Gifford (2012) [PubMed PMID: 22537047] and the GitHub repository
at <https://github.com/matted/multipool>.

Released under the MIT license as part of a version of the Multipool package.
Multipool prep script by Thomas Walsh. Multipool package by Matt Edwards.
"""

# TODO: handle regions that are part of a chromosome (e.g. 'chr04:1-200000').
# TODO: handle multiple regions?

import argparse
from collections import defaultdict
from collections import OrderedDict
from contextlib import contextmanager
import os
import platform
import pysam
import sys

from about import about

################################################################################

known_allele_depth_tagsets = (
    ("RO", "AO"), # used by FreeBayes
    ("AD",)       # used by GATK and version 4.3+ of the VCF specification
)

known_allele_depth_tags = {

    "AD": OrderedDict([ # GATK
        ("ID", "AD"),
        ("Number", "."),
        ("Type", "Integer"),
        ("Description", "Allelic depths for the ref and alt "
            "alleles in the order listed")
    ]),

    "AO": OrderedDict([ # FreeBayes
        ("ID", "AO"),
        ("Number", "A"),
        ("Type", "Integer"),
        ("Description", "Alternate allele observations, "
            "with partial observations recorded fractionally")
    ]),

    "RO": OrderedDict([ # FreeBayes
        ("ID", "RO"),
        ("Number", "1"),
        ("Type", "Integer"),
        ("Description", "Reference allele observation count, "
            "with partial observations recorded fractionally")
    ])
}

################################################################################

@contextmanager
def bust_open(name=None, mode='r'):
    """Robustly open file for reading or writing. If hyphen ('-') is specified
    or no filename is given, open standard input or output, respectively.

    Based on the solution at <stackoverflow.com/questions/17602878>
    by Stack Overflow user Wolph <stackoverflow.com/users/54017>.
    """

    if mode == 'r':
        default_handle = sys.stdin
    elif mode == 'w':
        default_handle = sys.stdout
    else:
        raise ValueError("unsupported file open mode: {!r}".format(mode))

    if name not in (None, '-'):
        fh = open(name=name, mode=mode)
    else:
        fh = default_handle

    try:
        yield fh
    finally:
        if fh is not default_handle:
            fh.close()

def cdl(argument, dtype=unicode):
    """Validate command-line argument as a comma-delimited list."""
    try:
        values = map( dtype, argument.split(',') )
    except AttributeError:
        raise ValueError("argument is not a comma-delimited list")
    except ValueError:
        raise ValueError("argument cannot be converted to a list "
            "of type {!r}".format(dtype.__name__))
    return values

def doOutput(allele_depth_info, output_file):
    """Output allele depth file."""
    with bust_open(output_file, 'w') as fh:
        for contig in allele_depth_info:
            for pos in sorted( allele_depth_info[contig] ):
                data = [pos] + allele_depth_info[contig][pos]
                line = ' '.join( str(x) for x in data )
                fh.write('{}\n'.format(line))

def get_allele_depth_info(pool, region, sample_file_info, founder_allele_info):
    """Get dictionary of allele depth info from input variant files."""

    allele_depth_info = OrderedDict()
    expected_tagset = None

    # Get allele depth of each observed biallelic
    # pool SNP in the given region(s) of interest.
    for input_file in sample_file_info[pool]:

        with pysam.VariantFile(input_file) as reader:

            d = OrderedDict()

            for record in reader:

                alleles = record.alleles
                contig = record.contig
                pos = record.pos

                # Skip if variant is not in a region of interest.
                if contig != region:
                    continue

                # Check for multiple records describing the same variant.
                if contig in d and pos in d[contig]:
                    raise RuntimeError("multiple records describe the same "
                        "variant at '{}' in pool sample '{}'".format(
                        '{}:{}'.join(contig, pos), pool))

                # Skip if variant is not a SNP.
                if any( a not in ('A', 'C', 'G', 'T') for a in alleles ):
                    continue

                # Get full set of genotype subfield tags for this record.
                record_tagset = set(record.format)

                # Identify allele depth tagset for this record,
                # then get the allele depths for the pool sample.
                if 'AD' in record_tagset:

                    ad_tagset = ('AD',)

                    allele_depths = record.samples[pool]['AD']

                elif all( tag in record_tagset for tag in ('RO', 'AO') ):

                    ad_tagset = ('RO', 'AO')

                    ref_obs = record.samples[pool]['RO']
                    alt_obs = record.samples[pool]['AO']
                    allele_depths = (ref_obs,) + alt_obs

                else:
                    raise RuntimeError("cannot find a known allele depth "
                        "tagset in pool variant data")

                # Check that allele depth tagsets are consistent.
                if expected_tagset is None:
                    expected_tagset = ad_tagset
                elif ad_tagset != expected_tagset:
                    raise RuntimeError("inconsistent allele depth tagsets "
                        "found in pool variant data")

                # Set non-zero allele depths for this variant.
                for allele, depth in zip(alleles, allele_depths):
                    if depth > 0:
                        d.setdefault(contig, defaultdict(lambda: defaultdict(int)))
                        d[contig][pos][allele] += depth

    # Get the allele depth of each pool variant
    # in terms of the matching founder allele.
    for contig in d:

        for pos in d[contig]:

            # Skip if no founder variant at this locus.
            if contig not in founder_allele_info or pos not in founder_allele_info[contig]:
                continue

            # Get founder alleles at this variant locus.
            founder_alleles = founder_allele_info[contig][pos]

            # Get variant allele depths with respect to each founder allele.
            variant_depths = list()
            for founder_allele in founder_alleles:
                try:
                    allele_depth = d[contig][pos].pop(founder_allele)
                except KeyError: # no depth available for founder allele
                    continue
                variant_depths.append(allele_depth)

            # Skip if not all founder alleles are represented.
            if len(variant_depths) != len(founder_alleles):
                continue

            # Skip if non-founder alleles are represented.
            if sum( d[contig][pos].values() ) > 0:
                continue

            allele_depth_info.setdefault(contig, dict())
            allele_depth_info[contig][pos] = variant_depths

    # Check that pool allele depths were found for region(s) of interest.
    for region in regions:
        if region not in allele_depth_info:
            raise RuntimeError("pool allele depths not found "
                "for region '{}'".format(region))

    return allele_depth_info

def get_founder_allele_info(founders, region, sample_file_info):
    """Get dictionary of founder allele info from input variant files."""

    founder_allele_info = defaultdict( lambda: defaultdict(set) )

    # If two founder samples are specified, get homozygous genotypes from
    # the given regions of interest that are segregating with respect to
    # the specified founder samples..
    if len(founders) == 2:

        d = defaultdict( lambda: defaultdict( lambda: defaultdict(set)) )

        for founder in founders:

            for input_file in sample_file_info[founder]:

                with pysam.VariantFile(input_file) as reader:

                    for record in reader:

                        alleles = record.alleles
                        contig = record.contig
                        pos = record.pos

                        # Skip if variant is not in a region of interest.
                        if contig != region:
                            continue

                        # Check for multiple records describing the same variant.
                        if contig in d and pos in d[contig] and founder in d[contig][pos]:
                            raise RuntimeError("multiple records describe the "
                                "same variant at '{}' in founder '{}'".format(
                                '{}:{}'.join(contig, pos), founder))

                        # Skip if variant is not a SNP.
                        if any( a not in ('A', 'C', 'G', 'T') for a in alleles ):
                            continue

                        try: # Get founder allele indices.
                            allele_indices = [ i for i in set(
                                record.samples[founder]['GT'] ) ]
                        except KeyError: # Skip if no genotype subfield.
                            continue

                        try: # Get founder allele symbols.
                            founder_alleles = [ alleles[i]
                                for i in allele_indices ]
                        except TypeError: # Skip if no genotype call.
                            continue

                        # Skip if founder genotype is not homozygous.
                        if len(founder_alleles) != 1:
                            continue

                        # Set founder allele symbol for this variant.
                        d[contig][pos][founder] = founder_alleles[0]

        # Keep variants for which each of the founder
        # samples has a different homozygous genotype.
        for contig in d:

            for pos in d[contig]:

                if all( f in d[contig][pos] for f in founders ):

                    homolog_alleles = [ list(d[contig][pos][f])[0]
                        for f in founders ]

                    if len(set(homolog_alleles)) == len(homolog_alleles):
                        founder_allele_info[contig][pos] = homolog_alleles

    else:
        raise RuntimeError("cannot get founder alleles for {} "
            "founder samples".format( len(founders) ))

    # Check that founder alleles were found for region(s) of interest.
    for region in regions:
        if region not in founder_allele_info:
            raise RuntimeError("founder alleles not found "
                "for region '{}'".format(region))

    return founder_allele_info

def get_sample_file_info(samples, input_files):
    """Get mapping of sample names to variant files."""

    sample_file_info = defaultdict(list)

    for input_file in input_files:
        with pysam.VariantFile(input_file) as reader:
            for sample in samples:
                if sample in reader.header.samples:
                    sample_file_info[sample].append(input_file)

    for sample in samples:
        if sample not in sample_file_info:
            raise RuntimeError("sample '{}' not found "
                "in input variant data".format(sample))

    return sample_file_info

def parseArgs():
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("files", nargs="+", metavar="FILE",
        help="Input variant file[s] in VCF format.")

    parser.add_argument("-p", "--pool", metavar="STR",
        required=True, help="Pool sample name.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--founders", metavar="LIST", type=cdl,
        help="Comma-separated list of founder sample names. "
        "Two founders must be specified.")
    group.add_argument("-F", "--founders-file", metavar="FILE",
        help="Text file listing founder sample names, one per line. "
        "Two founders must be specified.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-r", "--regions", metavar="LIST", type=cdl,
        help="Comma-separated list of regions for which to create a "
        "Multipool allele depth file. One region must be specified, "
        "and this must refer to an entire chromosome/contig.")
    group.add_argument("-R", "--regions-file", metavar="FILE",
        help="Text file listing regions for which to create a Multipool "
        "allele depth file, one per line. One region must be specified, "
        "and this must refer to an entire chromosome/contig.")

    parser.add_argument("-o", "--output", dest="outfile", metavar="FILE",
        type=str, help="Output file of allele counts.")

    parser.add_argument("-v", "--version", action="version",
        version="%(prog)s {}".format(about['version']))

    return parser.parse_args()

################################################################################

if __name__ == "__main__":

    args = parseArgs()

    input_files = sorted( set([ os.path.realpath(x) for x in args.files ]) )

    pool = args.pool

    if args.founders_file is not None:
        with open(args.founders_file) as fh:
            founders = [ x for x in fh.read().splitlines() if x != '' ]
    else:
        founders = [ x for x in args.founders ]

    if args.regions_file is not None:
        with open(args.regions_file) as fh:
            regions = [ x for x in fh.read().splitlines() if x != '' ]
    else:
        regions = [ x for x in args.regions ]

    output_file = args.outfile

    print >> sys.stderr, "Python version: {}".format(platform.python_version())
    print >> sys.stderr, "Multipool version: {}".format(about['version'])
    print >> sys.stderr, "PySAM version: {}".format(pysam.__version__)
    print >> sys.stderr, "Input variant file(s): {}".format( os.pathsep.join(input_files) )
    print >> sys.stderr, "Pool sample: {}".format(pool)
    print >> sys.stderr, "Founder sample(s): {}".format( ','.join(founders) )
    print >> sys.stderr, "Region(s): {}".format( ','.join(regions) )
    print >> sys.stderr, "Output allele depth file: {}".format(output_file)

    for input_file in input_files:
        if not os.path.exists(input_file):
            raise RuntimeError("input file not found: {}".format(input_file))

    if pool in founders:
        raise RuntimeError("cannot set sample '{}' as both pool and founder".format(pool))

    if len(founders) != 2:
        raise ValueError("please specify two founder samples")

    if len(regions) != 1:
        raise ValueError("please specify one region for which to prepare Multipool input")

    region = regions[0]

    print >> sys.stderr, "Checking for pool and founder samples in input variant data..."
    sample_file_info = get_sample_file_info(founders + [pool], input_files)

    print >> sys.stderr, "Getting founder allele info..."
    founder_allele_info = get_founder_allele_info(founders, region, sample_file_info)

    print >> sys.stderr, "Getting pool allele depths..."
    allele_depth_info = get_allele_depth_info(pool, region,
        sample_file_info, founder_allele_info)

    print >> sys.stderr, "Writing allele depth file..."
    doOutput(allele_depth_info, output_file)

    print >> sys.stderr, "Done."

################################################################################
