#!/usr/bin/env python3
import argparse
import collections
import datetime
import logging
import pathlib
import sys
from bfx import getreads
from bfx import samreader

DUAL_BASES = {'R':'AG', 'Y':'CT', 'S':'GC', 'W':'AT', 'K':'GT', 'M':'AC'}
REVCOMP_TABLE = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
DESCRIPTION = """Correlate ambiguous bases in reads with their positions in alignments."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('bam', type=pathlib.Path)
  options.add_argument('fq1', type=pathlib.Path)
  options.add_argument('fq2', type=pathlib.Path)
  options.add_argument('-o', '--outformat', choices=('vcf', 'tsv'), default='tsv')
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

  if args.outformat == 'vcf':
    print(create_header(args.refname))

  missing = found = 0
  for dual_data, read, align in collate_dual_data(args.bam, dual_readses):
    if dual_data.ref_coord is None:
      missing += 1
    else:
      found += 1
      print(dual_data.format(args.outformat))
  logging.info(f'Info: {found} dual bases with a reference coordinate, {missing} without.')


DUAL_DATA_FIELDS = ('ref_name', 'read_id', 'mate', 'ref_coord', 'read_coord', 'alt1', 'alt2')
class DualData(collections.namedtuple('DualData', DUAL_DATA_FIELDS)):
  __slots__ = ()
  @classmethod
  def from_raw(cls, raw_dual, align, read):
    read_coord, ref_coord, read_base, align_base = raw_dual
    alt1, alt2 = DUAL_BASES[read_base]
    return cls(align.rname, read.id, align.mate, ref_coord, read_coord, alt1, alt2)
  def format(self, format_):
    if format_ == 'vcf':
      strs = (self.ref_name, str(self.ref_coord), '.', self.alt1, self.alt2, '.', '.', '.', '.', '.')
    elif format_ == 'tsv':
      strs = (
        self.ref_name, self.read_id, str(self.mate), str(self.ref_coord), str(self.read_coord),
        self.alt1, self.alt2
      )
    return '\t'.join(strs)


def collate_dual_data(bam_path, dual_readses):
  for read, align, raw_duals in get_raw_dual_data(bam_path, dual_readses):
    for raw_dual in raw_duals:
      dual_data = DualData.from_raw(raw_dual, align, read)
      yield dual_data, read, align


def get_raw_dual_data(bam_path, dual_readses):
  for read, align in get_matched_reads_and_alignments(bam_path, dual_readses):
    duals = match_duals_with_ns(read, align)
    if duals:
      yield read, align, duals


def get_dual_reads_dict(fq_path):
  return {read.id: read for read in get_dual_reads(getreads.getparser(fq_path, 'fastq'))}


def get_dual_reads(reads):
  for read in reads:
    if seq_has_duals(read.seq):
      yield read


def seq_has_duals(seq):
  for base in seq:
    if base in DUAL_BASES:
      return True
  return False


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


def match_duals_with_ns(read, align):
  read_seq = read.seq.upper()
  align_seq_raw = get_padded_seq(align).upper()
  if align.reverse:
    align_seq = get_revcomp(align_seq_raw)
  else:
    align_seq = align_seq_raw
  if len(read_seq) != len(align_seq):
    raise RuntimeError(
      f'Read {read.name} has a different computed length than alignment {align.qname} '
      f'({align.cigar}):\n  {read_seq}\n  {align_seq}'
    )
  duals = []
  for i, (read_base, align_base) in enumerate(zip(read_seq, align_seq)):
    if align.reverse:
      read_coord = len(read_seq) - i
    else:
      read_coord = i + 1
    ref_coord = align.to_ref_coord(read_coord)
    if read_base != align_base and align_base != '*':
      if read_base in DUAL_BASES and align_base == 'N':
        duals.append((read_coord, ref_coord, read_base, align_base))
      else:
        raise RuntimeError(
          f'Read {read.name} has an unexpected discrepancy with alignment {align.qname}: '
          f'{read_base} ({read_coord}) != {align_base} ({ref_coord})'
        )
  return duals


def get_padded_seq(align):
  """Replace hard-clipped bases with '*' and return a sequence the length of the original read."""
  align_seq = align.seq
  if align._cigar_list:
    oplen, op = align._cigar_list[0]
    if op == 'H':
      align_seq = '*'*oplen + align_seq
    if len(align._cigar_list) > 1:
      oplen, op = align._cigar_list[-1]
      if op == 'H':
        align_seq += '*'*oplen
  return align_seq


def create_header(refname=None):
  dt = datetime.datetime.now()
  lines = [
    '##fileformat=VCFv4.0',
    f'##fileDate={dt.year}{dt.month:02d}{dt.day:02d}',
    '##source=ambig-to-vcf.py'
  ]
  if refname:
    lines.append(f'##reference={refname}')
  lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample')
  return '\n'.join(lines)


def print_read_by_alignment(read, align, read_seq, align_seq):
  for read_coord, (read_base, align_base) in enumerate(zip(read_seq, align_seq), 1):
    ref_coord = align.to_ref_coord(read_coord)
    if ref_coord is None:
      ref_coord_str = 'None'
    else:
      ref_coord_str = f'{ref_coord:4d}'
    print(f'{ref_coord_str} ← {read_coord:3d}: {align_base} ← {read_base}')


def get_revcomp(seq):
  return seq.translate(REVCOMP_TABLE)[::-1]


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
