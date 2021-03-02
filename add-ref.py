#!/usr/bin/env python3
import argparse
import dataclasses
import importlib
import logging
import sys
from bfx import getreads
import duallib

DESCRIPTION = """Add the reference base to dual sites produced by locate-duals.py."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('ref', metavar='ref.fa', type=argparse.FileType('r'),
    help='The reference sequence.')
  options.add_argument('duals', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
    help='The duals. In tsv format (the output of locate-duals.py). Omit to read from stdin.')
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

  ref = parse_ref(args.ref)

  for dual in add_ref_bases(duallib.parse_tsv(args.duals), ref):
    print(dual.format())


def add_ref_bases(duals, ref):
  for dual in duals:
    chrom = ref[dual.ref_name]
    ref_base = chrom.seq[dual.ref_coord-1]
    yield dataclasses.replace(dual, ref_base=ref_base)


def parse_ref(ref_file):
  ref = {}
  for chrom in getreads.getparser(ref_file, 'fasta'):
    ref[chrom.name] = chrom
  return ref


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
