from itertools import islice
from collections import Counter
import gzip


class QualityControl:
    def __init__(self, path):
        self.path = path
        self.li = []
        self.reads = []
        self.repeats = 0
        self.n_rates = []
        self.avg = 0
        self.gc_avg = 0
        self.n_per_seq = 0
        self.reads = self.get_avg_read()

    def count_gc(self, read):
        gc = (read.count('G') + read.count('C')) / len(read)
        return round(gc, 4)

    def get_avg_read(self):
        with gzip.open(self.path, 'r') as f:
            contents = f.read()  # one string
            data = contents.decode('utf-8').splitlines()
            for line in islice(data, 1, None, 4):
                length = len(line.strip())
                self.li.append(length)
                self.reads.append(line.strip())
        print(self.li)
        if len(self.li) == 0:
            exit()
        self.avg = round(sum(self.li) / len(self.li))
        print('self.avg', self.avg)
        self.li.sort()
        c = Counter(self.li)
        total_gc = 0
        for x in self.reads:
            total_gc += self.count_gc(x)
        temp = total_gc / len(self.reads)
        self.gc_avg = round(temp * 100, 2)
        return self.reads

    def get_repeats(self, reads):
        uniques = set(reads)
        self.repeats = len(reads) - len(uniques)

        for r in reads:
            if 'N' in r:
                x = r.count('N')
                n_rate = x / len(r)
                self.n_rates.append(n_rate)
        self.n_per_seq = round(sum(self.n_rates) / len(reads), 4)
        return self.repeats, len(self.n_rates)

    def print_info(self):
        print(f'Reads in the file = {len(self.li)}:')
        print(f'Reads sequence average length = {self.avg}\n')
        print('Repeats = ', self.repeats)
        print('Reads with Ns = ', len(self.n_rates), '\n')
        print(f'GC content average = {self.gc_avg}%')
        print(f'Ns per read sequence = {self.n_per_seq * 100: .2f}%')

    def main(self):
        self.get_repeats(self.reads)
        return self.repeats


if __name__ == '__main__':
    read1 = QualityControl(input())
    read2 = QualityControl(input())
    read3 = QualityControl(input())
    info1 = sum(read1.get_repeats(read1.reads))
    info2 = sum(read2.get_repeats(read2.reads))
    info3 = sum(read3.get_repeats(read3.reads))
    infos = [info1, info2, info3]
    best = min(infos)
    if best == info1:
        read1.print_info()
    elif best == info2:
        read2.print_info()
    else:
        read3.print_info()
