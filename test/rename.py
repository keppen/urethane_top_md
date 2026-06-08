import sys
from collections import defaultdict

if len(sys.argv) != 3:
    print("Usage: script.py input.pdb output.pdb")
    sys.exit(1)

inp = sys.argv[1]
out = sys.argv[2]

counters = defaultdict(int)

with open(inp) as fin, open(out, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16].strip()

            # element = first letter (works for your case)
            element = atom_name[0]

            counters[element] += 1
            new_name = f"{element}{counters[element]}"

            # PDB atom name field = columns 13–16 (width 4, right-aligned usually)
            new_name_fmt = new_name.rjust(4)

            line = line[:12] + new_name_fmt + line[16:]

        fout.write(line)
