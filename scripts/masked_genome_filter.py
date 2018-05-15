#!/usr/bin/env python3

import sys

def scaffold_filter(record, min_len):
    length  = len(record[1])
    N_num   = record[1].count("N")
    N_ratio = N_num / (length + 1)
    if record[1] == "":
        return None
    elif length < min_len:
        return None
    elif N_ratio > 0.3:
        return None
    else:
        print(">{}\n{}".format(record[0], record[1]))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit()

    filepath = sys.argv[1]
    record = ["",""]
    min_len = 2000

    for line in open(filepath, 'r'):
        if line.startswith(">"):
            scaffold_filter(record, min_len)
            record[0] = line.strip()[1:]
            record[1] = ""
        else:
            record[1] = record[1] + line.strip()
