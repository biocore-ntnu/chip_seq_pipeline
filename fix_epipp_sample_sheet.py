handle = open("sample_sheet.txt")
lines = []
lines.append(next(handle).strip())

mods = ["H3K27me3", "H3K4me3", "PolII"]

for l in handle:
    if "Input" in l:
        f, n, g, c, m = l.split()
        for mod in mods:
            nt = n.replace("Input", mod + "_Input")
            gt = g.replace("Input", mod)
            lines.append(" ".join([f, nt, gt, c, m]))
    else:
        lines.append(l.strip())

from operator import itemgetter
get_name_group = itemgetter(1, 2)
sorted_lines = [lines[0]] + sorted(lines[1:], key=lambda l: get_name_group(l.split()))
print("\n".join(sorted_lines))
