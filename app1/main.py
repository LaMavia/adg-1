from collections import defaultdict
from mmh3 import hash
from tqdm import tqdm
from Bio import SeqIO
from sys import argv
from array import array
from contextlib import contextmanager
from time import perf_counter


@contextmanager
def catchtime(label: str):
    t1 = t2 = perf_counter()
    yield
    t2 = perf_counter()
    print(f'{label} took {t2 - t1:.3f} s')

def kmer_hash(kmer: str) -> int:
    return hash(kmer)

def get_minimizers(seq: str, k: int = 15, w: int = 10):
    """
    Stream (yield) minimizers instead of building a list.
    Yields: (window_position, minimizer_hash)
    """
    for i in range(len(seq) - w + 1):
        window = seq[i:i + w]
        minim = min(kmer_hash(window[j:j+k]) for j in range(w - k + 1))
        yield i, minim

class MashMapIndex:
    """Hierarchical minimizer index."""
    def __init__(self, ref_seq: str, k: int = 15, windows=[10,20,40]):
        self.k = k
        self.windows = windows
        self.index = defaultdict(lambda: array('I'))
        self.ref_seq = ref_seq
        self.build_index(ref_seq)

    def build_index(self, ref_seq: str):
        for w in self.windows:
            for pos, m in get_minimizers(ref_seq, k=self.k, w=w):
                self.index[(w, m)].append(pos)

    def query(self, read_seq: str):
        """
        Now uses streaming minimizers.
        """
        hits = defaultdict(int)

        for w in self.windows:
            for read_pos, m in get_minimizers(read_seq, k=self.k, w=w):
                ref_positions = self.index.get((w, m))
                if not ref_positions:
                    continue
                for ref_pos in ref_positions:
                    offset = ref_pos - read_pos
                    hits[offset] += 1

        return hits


def mashmap_map_read(read_seq: str, index: MashMapIndex, min_hits_ratio=0.2):
    hits = index.query(read_seq)

    first_w = index.windows[0]
    minimizer_count = sum(1 for _ in get_minimizers(read_seq, k=index.k, w=first_w))

    min_hits = max(1, int(minimizer_count * min_hits_ratio))

    try:
        candidate_positions = [
            next(pos for pos, count in hits.items() if count >= min_hits)
        ]
    except StopIteration:
        candidate_positions = []

    return candidate_positions


def main():
    seq_rec = next(SeqIO.parse(argv[1], "fasta"))
    genome = str(seq_rec.seq)
    print("Building MashMap index...")
    with catchtime("mashmap index"):
        index = MashMapIndex(genome, k=15, windows=[60])

    reads = list(SeqIO.parse(argv[2], "fasta"))
    fout = open(argv[3], "w")

    for read in tqdm(reads):
        read_seq = str(read.seq)
        positions = mashmap_map_read(read_seq, index, min_hits_ratio=0.2)
        for pos in positions:
            fout.write(f"{read.id}\t{pos}\t{pos + len(read_seq)}\n")
            break  # report only first candidate for simplicity
    fout.close()


if __name__ == "__main__":
    main()