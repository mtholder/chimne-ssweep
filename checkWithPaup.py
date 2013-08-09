#!/usr/bin/env python
import sys
import subprocess
import re
lnLPat = re.compile(r"^-ln\sL\s+([.0-9]+)")
f = open("paupStdOut", "w")
subprocess.call(["paup", "-n", "masterPaup.nex"], stdout=f)
f = open("paup.log", "rU")
first = None
for line in f:
    m = lnLPat.match(line)
    if m:
        print 
        if first is None:
            first = float(m.groups(1)[0])
        else:
            second = float(m.groups(1)[0])
            if abs(first-second) > 10e-4:
                sys.exit(1)
            print "It worked!"    
            sys.exit(0)
sys.exit(1)
