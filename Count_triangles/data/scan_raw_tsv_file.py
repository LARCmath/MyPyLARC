
import sys
import math

def doit(filename):
    count = 0
    v1min = 10000000
    v1max = -1
    v2min = v1min
    v2max = v1max
    eup = 0
    edown = 0
    with open(filename, "r") as fp:
        for x in fp:
            y = x.split()
            assert len(y) == 3, f"y = {y}"
            assert int(y[2]) == 1, f"y = {y}"
            v1min = min(v1min, int(y[0]))
            v1max = max(v1max, int(y[0]))
            v2min = min(v2min, int(y[1]))
            v2max = max(v2max, int(y[1]))
            if y[0] < y[1]:
                eup += 1
            else:
                edown += 1
            count += 1
    assert (v1min == 1) and (v2min == 1)
    assert v1max == v2max
    assert eup == edown
    assert filename[-4:] == ".tsv"
    level = math.ceil(math.log(v1max) / math.log(2))
    outfilename = filename[:-4] + "_header" + str(level) + ".txt"
    with open(outfilename, "w") as outfp:
        print("%%MatrixMarket matrix coordinate integer general", file=outfp)
        print(f"{v1max} {v2max} {count}", file=outfp)


# read_matrixMarketExchange_file_return_matPTR("filename")

doit(sys.argv[1])

