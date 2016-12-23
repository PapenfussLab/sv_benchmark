#!/bin/bash
Rscript precache.R rd &
Rscript precache.R rl &
Rscript precache.R fs &
Rscript precache.R na12878 &
wait