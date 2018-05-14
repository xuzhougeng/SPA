import sys

if len(sys.argv) < 2:
    sys.exit()

fp = open(sys.argv[1], 'r')

last_state = ""
pre_seqid  = ""
pre_start  = "0"
pre_end    = "0"

for line in  fp:
    items = line.strip().split("\t")
    now_seqid = items[0]
    now_start = items[1]
    now_end   = items[2]
    # reset the position when seqid change
    if not now_seqid == pre_seqid:
        pre_start = "0"    
        pre_end   = "0"
    # if current start is large then the end of last line
    if  int(now_start) > int(pre_end): 
        print(line.strip())
        pre_seqid = items[0]
        pre_start = items[1]
        pre_end   = items[2]
fp.close()
