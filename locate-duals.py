#!/usr/bin/env python3
import argparse
import collections
import datetime
import logging
import pathlib
import sys
from bfx import getreads
from bfx import samreader
import duallib

PAD_CHAR = '*'
REVCOMP_TABLE = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
DESCRIPTION = """Correlate ambiguous bases in reads with their positions in alignments."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('bam', type=pathlib.Path)
  options.add_argument('fq1', type=pathlib.Path)
  options.add_argument('fq2', type=pathlib.Path)
  options.add_argument('-h', '--help', action='help',
    help='Print this argument help text and exit.')
  logs = parser.add_argument_group('Logging')
  logs.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = logs.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  
  dual_readses = [get_dual_reads_dict(fq_path) for fq_path in (args.fq1, args.fq2)]
  for mate, dual_reads in enumerate(dual_readses, 1):
    logging.info(f'Info: Mate {mate}: {len(dual_reads)} reads with dual bases.')

  missing = found = 0
  for read, align, dual in get_duals(args.bam, dual_readses):
    if dual.ref_coord is None:
      missing += 1
    else:
      found += 1
      print(dual.format())
  logging.info(f'Info: {found} dual bases with a reference coordinate, {missing} without.')


def get_duals(bam_path, dual_readses):
  for read, align, duals in get_duals_per_read(bam_path, dual_readses):
    for dual in duals:
      yield read, align, dual


def get_duals_per_read(bam_path, dual_readses):
  for read, align in get_matched_reads_and_alignments(bam_path, dual_readses):
    try:
      duals = list(duallib.Dual.from_raws(filter_for_duals(*align_sites(read, align)), read, align))
    except RuntimeError:
      logging.critical(
        f'Error comparing read {read.name} with alignment {align.qname}/{align.mate} '
        f'({align.cigar})'
      )
      raise
    if duals:
      yield read, align, duals


def get_dual_reads_dict(fq_path):
  return {read.id: read for read in duallib.get_dual_reads(getreads.getparser(fq_path, 'fastq'))}


def get_matched_reads_and_alignments(bam_path, dual_readses):
  stats = collections.Counter()
  for align in samreader.read_bam(bam_path):
    if align.secondary:
      stats['secondary'] += 1
      continue
    elif align.duplicate:
      stats['duplicate'] += 1
      continue
    elif align.supplemental:
      stats['supplemental'] += 1
      continue
    try:
      read = dual_readses[align.mate-1][align.qname]
    except KeyError:
      continue
    yield read, align
  for key, count in stats.items():
    if count > 0:
      logging.info(f'Info: Skipped {count} {key} alignments.')


def align_sites(read, align):
  read_coords = list(range(1, len(read.seq)+1))
  ref_coords = [align.to_ref_coord(coord) for coord in read_coords]
  read_seq = read.seq.upper()
  align_seq_padded = get_padded_seq(align).upper()
  if align.reverse:
    read_bases = get_complement(read_seq)
    align_bases = ''.join(reversed(align_seq_padded))
  else:
    read_bases = read_seq
    align_bases = align_seq_padded
  if not len(read_coords) == len(ref_coords) == len(read_bases) == len(align_bases):
    raise RuntimeError(
      f'Read has a different computed length than alignment:\n'
      f'  {read_seq}\n'
      f'  {align_seq_padded}'
    )
  return read_coords, ref_coords, read_bases, align_bases


def filter_for_duals(read_coords, ref_coords, read_bases, align_bases):
  for read_coord, ref_coord, read_base, align_base in zip(
    read_coords, ref_coords, read_bases, align_bases
  ):
    if read_base != align_base and align_base != PAD_CHAR:
      if read_base in duallib.DUAL_BASES and align_base == 'N':
        yield read_coord, ref_coord, read_base, align_base
      else:
        raise RuntimeError(
          f'Read has an unexpected discrepancy with alignment: '
          f'{read_base} ({read_coord}) != {align_base} ({ref_coord})'
        )


def get_padded_seq(align):
  """Replace hard-clipped bases with PAD_CHAR and return a sequence the length of the original read.
  """
  align_seq = align.seq
  if align._cigar_list:
    oplen, op = align._cigar_list[0]
    if op == 'H':
      align_seq = PAD_CHAR*oplen + align_seq
    if len(align._cigar_list) > 1:
      oplen, op = align._cigar_list[-1]
      if op == 'H':
        align_seq += PAD_CHAR*oplen
  return align_seq


def print_read_by_alignment(read, align, read_seq, align_seq):
  for read_coord, (read_base, align_base) in enumerate(zip(read_seq, align_seq), 1):
    ref_coord = align.to_ref_coord(read_coord)
    if ref_coord is None:
      ref_coord_str = 'None'
    else:
      ref_coord_str = f'{ref_coord:4d}'
    print(f'{ref_coord_str} ← {read_coord:3d}: {align_base} ← {read_base}')


def get_complement(seq):
  return seq.translate(REVCOMP_TABLE)


def get_revcomp(seq):
  return get_complement(seq)[::-1]


def fail(message):
  logging.critical(f'Error: {message}')
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
