# -*- coding: utf-8 -*-
"""
chanjo.annotator.core
~~~~~~~~~~~~~~~~~~~~~~

Central pipeline for the Chanjo annotator module.
"""
from __future__ import absolute_import, unicode_literals

from toolz import concat
from toolz.curried import do

from .._compat import text_type
from ..depth_reader import BamFile
from .stages import (
  calculate_metrics,
  comment_sniffer,
  extend_interval,
  group_intervals,
  prefix,
  process_interval_group
)
from ..utils import bed_to_interval, validate_bed_format


def annotate_bed_stream(bed_stream, bam_path, cutoff=10, extension=0,
                        contig_prefix='', bp_threshold=17000):
  """Annotate all intervals from a BED-file stream.

  Yields tuple data for each interval with calculated coverage and
  completeness.

  Args:
    bed_stream (sequence): usually a BED-file handle to read from
    bam_path (str): path to BAM-file
    cutoff (int, optional): threshold for completeness calculation,
      defaults to 10
    extension (int, optional): number of bases to extend each interval
      with (+/-), defaults to 0
    contig_prefix (str, optional): rename contigs by prefixing,
      defaults to empty string
    bp_threshold (int, optional): optimization threshold for reading
      BAM-file in chunks, default to 17000

  Yields:
    tuple: :class:`chanjo.BaseInterval`, coverage (float), and
      completeness (float)
  """
  # setup: connect to BAM-file
  bam = BamFile(bam_path)

  # filter out comments
  no_comments = (line for line in bed_stream if not comment_sniffer(line))

  # strip invisble characters
  cleaned_up = (text_type.rstrip(line) for line in no_comments)

  # prefix to contig
  prefixed = (prefix(contig_prefix, line) for line in cleaned_up)

  # split lines
  splitted = (text_type.split(line) for line in prefixed)

  # check correct format
  validated = (do(validate_bed_format, row) for row in splitted)

  # convert to objects
  objectified = (bed_to_interval(*row) for row in validated)

  # extend intervals
  extended = (extend_interval(row, extension=extension) for row in objectified)

  # group rows by threshold
  grouped = group_intervals(extended, bp_threshold=bp_threshold)

  # read coverage
  read_depths = (process_interval_group(bam, row) for row in grouped)

  # flatten list of lists
  flattened = concat(read_depths)

  # calculate coverage/completeness
  return (calculate_metrics(depths, threshold=cutoff) for depths in flattened)
