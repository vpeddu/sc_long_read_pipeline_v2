import sys, argparse, csv
from collections import defaultdict

def parse_commandline():
  default_fn = 'sqanti_corrected.genePred'
  parser=argparse.ArgumentParser()
  parser.add_argument('--genePred', '-g', help='gtf file in genePred format', type=str, default=default_fn, required=False)
  args=parser.parse_args()
  print(args, file=sys.stderr)
  return args

args = parse_commandline()
genes_fn = args.genePred

groups = defaultdict(list)

with open(genes_fn) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        key = tuple(fields[i] for i in [1,2,3,4,5,6,7,8,9,10,12,13,14])
        groups[key].append(fields)

with open("tx_dups_xref.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    for rows in sorted(groups.values(), key=len, reverse=True):
        if len(rows) > 1:
            # Sort rows by Ensembl transcript id (descending order)
            rows = sorted(rows, key=lambda r: r[0])
            last_tx = rows[-1]  # row with largest ENST transcript number

            for row in rows[:-1]:
                writer.writerow([row[0], row[11], last_tx[0], last_tx[11]])
