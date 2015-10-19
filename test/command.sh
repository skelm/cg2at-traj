#!/bin/bash
#
# Example command for cg2at.py:
#
# For information on all the options, use ../cg2at.py --help
#

../cg2at.py --select_protein --center_protein --skip 50 --processes -1 --out_prefix trajectory md.xtc md.tpr

#
# The minimal command is:
#
# ../cg2at.py md.xtc md.tpr
#
# but that would take a long time to run. You can add --skip 50 to speed it up and only consider every 50th frame.
#

