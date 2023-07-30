"""Map the official Hg38 contig names onto the chr1,chr2 etc that we use in our standard Hg38 reference on nanopore boxes onto any given file as passed by the command line

Example usage
--------------
python rename_contigs.py <file_in> <file_out> <direction>

<direction> can be any character. If present, we convert from official HG38 contig names - NC_000001.11 etc, to human readable - chr1 etc. Otherwise the conversion happens the other direction.
"""
import sys

hg38 = """NC_000001.11
NC_000002.12
NC_000003.12
NC_000004.12
NC_000005.10
NC_000006.12
NC_000007.14
NC_000008.11
NC_000009.12
NC_000010.11
NC_000011.10
NC_000012.12
NC_000013.11
NC_000014.9
NC_000015.10
NC_000016.10
NC_000017.11
NC_000018.10
NC_000019.10
NC_000020.11
NC_000021.9
NC_000022.11
NC_000023.11
NC_000024.10
NC_012920.1"""

human_readable = """chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY
chrM"""

if __name__ == "__main__":
    print(sys.argv[1])
    with open(sys.argv[1]) as fh, open(sys.argv[2], "w") as out:
        convert = (
            dict(zip(hg38.splitlines(), human_readable.splitlines()))
            if len(sys.argv) > 3
            else dict(zip(human_readable.splitlines(), hg38.splitlines()))
        )
        print(convert)
        for line in fh:
            print(line)
            delim = "\t" if "\t" in line else ","
            splitter = None if "\t" in line else ","
            line = list(map(lambda x: convert.get(x, x), line.strip().split(splitter)))
            out.write(delim.join(line) + "\n")
