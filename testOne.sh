#!/bin/bash
for f in ../combine/*.combine;
	do
timeout 200s /usr/bin/time -f"wtime=%e" ./pl $f &>$f.pl
	done
