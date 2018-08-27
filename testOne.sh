#!/bin/bash
for f in ../combine/*.combine;
	do
timeout 120 /usr/bin/time -f"wtime=%e" ./plf2 $f &>$f.plf2
	done
