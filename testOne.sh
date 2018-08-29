#!/bin/bash
for f in ../combine/*.combine;
	do
timeout 400s /usr/bin/time -f"wtime=%e" ./pl $f &>$f.pl1
	done
