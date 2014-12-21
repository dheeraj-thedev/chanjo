# -*- coding: utf-8 -*-
"""
chanjo.builder.core
~~~~~~~~~~~~~~~~~~~~

Central pipeline for the Chanjo builder module.
"""
from __future__ import absolute_import
import errno
import os

from path import path
from toolz import pipe, reduce, concat
from toolz.curried import map

from .._compat import text_type
from .consumers import commit_per_contig
from .stages import aggregate, build_block, build_interval, build_superblock
from ..utils import bed_to_interval, split


def init_db(chanjo_db, bed_stream, overwrite=False):
  """Build a new database instance from the Chanjo BED stream.

  Args:
    chanjo_db (Store): initialized Store class instance
    bed_stream (sequence): Chanjo-style BED-stream
    overwrite (bool, optional): whether to automatically overwrite an
      existing database, defaults to False
  """
  # check if the database already exists (expect 'mysql' to exist)
  # 'dialect' is in the form of '<db_type>+<connector>'
  if chanjo_db.dialect == 'mysql' or path(chanjo_db.uri).exists():
    if overwrite:
      # wipe the database clean with a warning
      chanjo_db.tear_down()
    elif chanjo_db.dialect == 'sqlite':
      # prevent from wiping existing database to easily
      raise OSError(errno.EEXIST, os.strerror(errno.EEXIST), chanjo_db.uri)

  # set up new tables
  chanjo_db.set_up()

  # strip invisble characters
  cleaned_up = (line.rstrip() for line in bed_stream)

  # split lines
  splitted = (text_type.split(line) for line in cleaned_up)

  # convert to objects
  objectified = (bed_to_interval(*row) for row in splitted)

  # build interval objects and load into database
  nested_intervals = (build_interval(chanjo_db, interval)
                      for interval in objectified)

  # flatten the nested list of lists of intervals
  flattened = concat(nested_intervals)

  # aggregate intervals based on common block IDs
  block_groups = aggregate(flattened)

  # build interval objects and load into database
  blocks = (build_block(chanjo_db, group) for group in block_groups)

  # aggregate blocks based on common superblock IDs
  superblock_groups = aggregate(blocks)

  # build interval objects and load into database
  superblocks = (build_superblock(chanjo_db, group)
                 for group in superblock_groups)

  # reduce the superblocks and commit every contig
  reduce(commit_per_contig(chanjo_db), superblocks, 'chr0')

  # commit also the last contig
  chanjo_db.save()
