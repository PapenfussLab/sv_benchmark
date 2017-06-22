#!/bin/bash
# Invokes all indel callers
#
. common.sh

for CALLER in $(ls call_*.sh); do
	if [ "$CALLER" != "call_all.sh" ] ; then
		C=${CALLER%.sh}
		echo "Calling ${C:5}"
		./$CALLER $1 $2 $3 $4 $5
	fi
done



