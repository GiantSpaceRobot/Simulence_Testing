#!/bin/bash

#Create directories
mkdir -p AllTests

# Loop 10 times
for i in {1..10}
do
	# your-unix-command-here
	echo "
	Loop #$i
	"
	# Human
	./Test_Simulence.sh human
	awk -F '\t' '{print $2}' human/All-human-ncRNAs_spiked_Simulated-Reads_Stats.txt > AllTests/human_$i.txt
	# Mouse
	./Test_Simulence.sh mouse
	awk -F '\t' '{print $2}' mouse/All-mouse-ncRNAs_spiked_Simulated-Reads_Stats.txt > AllTests/mouse_$i.txt
	# Mouse
	./Test_Simulence.sh rat
	awk -F '\t' '{print $2}' rat/All-rat-ncRNAs_spiked_Simulated-Reads_Stats.txt > AllTests/rat_$i.txt
done

paste AllTests/human*.txt > AllTests/human-results.txt
paste AllTests/mouse*.txt > AllTests/mouse-results.txt
paste AllTests/rat*.txt > AllTests/rat-results.txt

