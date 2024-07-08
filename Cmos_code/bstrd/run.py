#!/usr/bin/env python

from time import sleep
from bstrd import BS, bsstream


v0 = BS("pid")
v1 = BS("SATFE10-PEPG046:FCUP-INTENSITY-CAL")
#v2 = BS("SATES21-GES1:A1_VALUES")
#v3 = BS("SIN-CVME-TIFGUN-EVR0:RX-PULSEID")


#from collections import defaultdict
#h = defaultdict(int)


for i, data in enumerate(bsstream):
    sleep(0.005)

#    h[mc.queue.qsize()] += 1

#    print(v0.value)

#    if v0.value is not None and v0.value % 100 == 0:
#        for val, count in sorted(h.items()):
#            print(val, "x" * count if count <100 else "=" + str(count))
#        print()

    if i == 200:
        v2 = BS("SATES21-GES1:A1_VALUES")

    if i < 200:
        print(i, v0, v1)
    else:
        print(i, v0, v1, v2)


    if i == 400:
        break


for data in bsstream:
    print(data)


