#!/usr/bin/env python3
import argparse
import datetime
import logging
import sys
import duallib

DESCRIPTION = """"""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('duals', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
    help='The duals. In tsv format (the output of add-ref.py). Omit to read from stdin.')
  options.add_argument('--ref-name',
    help='The name of the reference sequence (for the VCF header).')
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

  print(create_header(args.ref_name))

  mismatches = 0
  for dual in duallib.parse_tsv(args.duals):
    try:
      ref_base, alt = get_ref_alt(dual.ref_base, dual.alt1, dual.alt2)
    except ValueError:
      mismatches += 1
      logging.info(
        f'Info: Dual {dual.alt1}/{dual.alt2} from {dual.read_id}/{dual.mate} at '
        f'{dual.ref_name}:{dual.ref_coord} does not match the reference base {dual.ref_base!r}.'
      )
      continue
    strs = (dual.ref_name, str(dual.ref_coord), '.', ref_base, alt, '.', '.', '.', '.', '.')
    print(*strs, sep='\t')
  
  logging.warning(f'Saw {mismatches} sites that did not match the reference.')


def create_header(refname=None):
  dt = datetime.datetime.now()
  lines = [
    '##fileformat=VCFv4.0',
    f'##fileDate={dt.year}{dt.month:02d}{dt.day:02d}',
    '##source=g4/to-vcf.py'
  ]
  if refname:
    lines.append(f'##reference={refname}')
  lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample')
  return '\n'.join(lines)


def get_ref_alt(ref_base, alt1, alt2):
  if ref_base == alt1:
    return ref_base, alt2
  elif ref_base == alt2:
    return ref_base, alt1
  else:
    raise ValueError(f'Neither alt ({alt1}/{alt2}) match the reference ({ref_base})')


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
