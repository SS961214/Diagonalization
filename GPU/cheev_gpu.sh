#!/bin/bash -x

date
time for n in $(seq 76 100);do ./testing_cheevd $((100*($n+50)));done
