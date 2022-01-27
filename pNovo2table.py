import sys

pnovo = sys.argv[1]

out_title = ["File_Name","Sequence","Score","Rank"]

psm_list = []
with open(pnovo) as f:
    while True:
        line = f.readline()
        if line == "": break
        if line.startswith("S"):
            filename = line.split("\t")[1]
        if line.startswith("P\t"):
            items = line.strip().split("\t")
            sequence = items[1]
            score = items[2]
            rank = items[0][1:]
            psm_list.append((filename, sequence, score, rank))
with open(pnovo[:-4] + "-table.txt","w") as f:
    f.write("\t".join(out_title) + "\n")
    for psm in psm_list:
        f.write("\t".join(psm) + "\n")
    