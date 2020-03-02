from functools import reduce
from typing import List, Optional, Dict, Union

from tools.resources import CODON_MAP, AA_WEIGHTS, START, STOP_CODONS, CODON_SIZE, COMPLEMENT_DNA, \
    COMPLEMENT_RNA


class FastaTypeError(Exception):
    pass


class UnknownCodonError(Exception):
    pass


DNA: str = 'DNA'
RNA: str = 'RNA'
PROTEIN: str = 'PROTEIN'


class Fasta:
    def __init__(self, id_: str = 'noid', sequence: str = ''):
        if ' ' in id_:
            raise ValueError('No spaces allowed in id')
        if sequence == '':
            raise ValueError('Sequence cannot be empty')
        self.id = id_
        self.sequence = sequence.upper()
        self.type: str = self._infer_type()

    def __str__(self) -> str:
        return f'>{self.id}\n{self.sequence}'

    def __repr__(self) -> str:
        return f'Fasta({self.id}, {self.sequence})'

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Fasta):
            return NotImplemented()
        return self.sequence == other.sequence

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, item):
        return self.sequence[item]

    @property
    def gc(self) -> float:
        c = self.sequence.count('C')
        g = self.sequence.count('G')
        total = len(self.sequence)
        return ((c + g) / total) * 100

    @classmethod
    def from_str(cls, s: str) -> List["Fasta"]:
        raw_seqs = s.split('>')[1:]
        fastas = []
        for rs in raw_seqs:
            id_, seq = rs.split('\n', 1)
            seq = ''.join(seq.split())  # Remove all whitespace/newlines
            fastas.append(cls(id_, seq))
        return fastas

    def count(self, s: str) -> int:
        return self.sequence.count(s)

    @staticmethod
    def _complement_sequence(seq: str, type_: str) -> str:
        if type_ not in [DNA, RNA]:
            raise FastaTypeError('Only Fasta with type == DNA or RNA can be complemented')
        comp_dict: Dict[str, str] = COMPLEMENT_DNA if type == DNA else COMPLEMENT_RNA
        return ''.join([comp_dict[base] for base in seq])

    def complement(self) -> "Fasta":
        comp = self._complement_sequence(self.sequence, type_=self.type)
        return Fasta(self.id + '_complement', comp)

    def reverse_complement(self) -> "Fasta":
        rev_comp = self._complement_sequence(''.join(reversed(self.sequence)), type_=self.type)
        return Fasta(self.id + '_reverse_complement', rev_comp)

    def transcribe(self) -> "Fasta":
        if self.type != DNA:
            raise FastaTypeError('Only Fasta with type == DNA can be transcribed')
        return Fasta(self.id + '_transcribed', self.sequence.replace('T', 'U'))

    def reverse_transcribe(self) -> "Fasta":
        if self.type != RNA:
            raise FastaTypeError('Only Fasta with type == RNA can be reverse transcribed')
        return Fasta(self.id + '_reverse_transcribed', self.sequence.replace('U', 'T'))

    def _translatable_sequences(self) -> Dict[str, List[str]]:
        orf1: List[str] = self._seq_to_codons(self.sequence)
        orf2: List[str] = self._seq_to_codons(self.sequence[1:])
        orf3: List[str] = self._seq_to_codons(self.sequence[2:])
        rc: 'Fasta' = self.reverse_complement()
        rc_orf1: List[str] = self._seq_to_codons(rc.sequence)
        rc_orf2: List[str] = self._seq_to_codons(rc.sequence[1:])
        rc_orf3: List[str] = self._seq_to_codons(rc.sequence[2:])
        return {'orf1': orf1, 'orf2': orf2, 'orf3': orf3, 'rc_orf1': rc_orf1, 'rc_orf2': rc_orf2, 'rc_orf3': rc_orf3}

    @staticmethod
    def _seq_to_codons(s) -> List[str]:
        r: int = len(s) % 3
        if r != 0:
            s = s[:-r]
        return [s[i:i + CODON_SIZE] for i in range(0, len(s), CODON_SIZE)]

    @staticmethod
    def _translate(codons: Optional[List[str]]) -> List[str]:
        if codons is None:
            return []

        stop_indices: List[int] = [i for i, x in enumerate(codons) if x in STOP_CODONS]
        if not stop_indices:
            return []
        start_indices: List[int] = [i for i, x in enumerate(codons) if x == START and i < max(stop_indices)]

        peptides: List[str] = []
        for start in start_indices:
            try:
                stop: int = [x for x in stop_indices if x > start][0]  # Find first stop codon after start codon
            except IndexError:
                continue
            peptides.append(''.join(CODON_MAP[c] for c in codons[start:stop]))

        return peptides

    def translate_all_orfs(self) -> Dict[str, List["Fasta"]]:
        f: Fasta = self._transcribe_if_DNA()
        trans_seq: Dict[str, List[str]] = f._translatable_sequences()
        peptides: Dict[str, List[Fasta]] = {}
        for frame, seq in trans_seq.items():
            if not seq:
                peptides[frame] = []
            else:
                # Create Fasta Objects from possible proteins in seq
                proteins = Fasta._translate(seq)
                peptides[frame] = [Fasta(id_=f'{f.id}_{frame}_{i+1}', sequence=x) for i, x in enumerate(proteins)]
        return peptides

    def translate_first_orf(self) -> List["Fasta"]:
        return self.translate_all_orfs()['orf1']

    def _transcribe_if_DNA(self) -> 'Fasta':
        if self.type == PROTEIN:
            raise FastaTypeError('Only Fasta with type == DNA or RNA can be translated')
        if self.type != RNA:
            return self.transcribe()
        else:
            return self

    def molecular_weight(self) -> float:
        if self.type != PROTEIN:
            raise FastaTypeError('Molecular weight can only be found for amino acid sequences')
        return reduce(lambda x, y: x+y, [AA_WEIGHTS[aa] for aa in self.sequence])

    def _infer_type(self) -> str:
        if set(self.sequence) == set('ATGC'):
            return DNA
        elif set(self.sequence) == set('AUGC'):
            return RNA
        else:
            # TODO: Do this more robustly, raise error on non-amino acid letters
            return PROTEIN

    def simple_splice(self, s: Union[str, 'Fasta'], force: bool = False) -> "Fasta":
        """Remove sequences s from self.sequence"""
        if isinstance(s, Fasta):
            s = s.sequence
        if s in self.sequence:
            new_seq = self.sequence.replace(s, '')
        elif not force:
            raise ValueError('Sequence cannot be spliced because it is not in Fasta')
        else:
            new_seq = s
        return Fasta(id_=f'{self.id}_spliced', sequence=new_seq)


if __name__ == "__main__":
    pass
