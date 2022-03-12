#!/bin/sh

for i in $(cat pdb.lst)
	do
	mkdir $i
	cd $i
	cp ../2ONC-cac-$i.pdb .
	~/programs/pdb2plif-0.0.2/pdb2plif.sh 2ONC-cac-$i.pdb ligand_CAC1_0.mol2
	rm 2ONC-cac-$i.pdb
	cd ..
	done
