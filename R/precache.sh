#!/bin/bash
Rscript precache.R --plot rd &
Rscript precache.R --plot rl &
Rscript precache.R --plot fs &
Rscript precache.R --plot na12878 &
wait
