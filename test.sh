timeout 2m /usr/bin/time -f "wtime=%e" ./pl1 $1 &> $1.testres1
timeout 2m /usr/bin/time -f "wtime=%e" ./pl2 $1 &> $1.testres2

