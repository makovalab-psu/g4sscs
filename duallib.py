import dataclasses

DUAL_BASES = {'R':'AG', 'Y':'CT', 'S':'GC', 'W':'AT', 'K':'GT', 'M':'AC'}
NULL_STR = '.'

@dataclasses.dataclass(frozen=True)
class Dual:
  ref_name: str
  read_id: str
  mate: int
  ref_coord: int
  read_coord: int
  alt1: str
  alt2: str
  ref_base: str = None
  @classmethod
  def from_strs(cls, raw_values, null_str=NULL_STR):
    values = []
    for (field, field_type), raw_value in zip(cls.__annotations__.items(), raw_values):
      if raw_value == null_str:
        value = None
      else:
        value = field_type(raw_value)
      values.append(value)
    return cls(*values)
  @classmethod
  def from_raw(cls, raw_dual, read, align):
    read_coord, ref_coord, read_base, align_base = raw_dual
    alt1, alt2 = DUAL_BASES[read_base]
    return cls(align.rname, read.id, align.mate, ref_coord, read_coord, alt1, alt2)
  @classmethod
  def from_raws(cls, raw_duals, read, align):
    for raw_dual in raw_duals:
      yield cls.from_raw(raw_dual, read, align)
  def format(self, sep='\t', null_str=NULL_STR):
    values = dataclasses.astuple(self)
    value_strs = [null_str if value is None else str(value) for value in values]
    return sep.join(value_strs)


def parse_tsv(lines):
  for line_raw in lines:
    raw_values = line_raw.rstrip('\r\n').split('\t')
    yield Dual.from_strs(raw_values)


def get_dual_reads(reads):
  for read in reads:
    if seq_has_duals(read.seq):
      yield read


def seq_has_duals(seq):
  for base in seq:
    if base in DUAL_BASES:
      return True
  return False
