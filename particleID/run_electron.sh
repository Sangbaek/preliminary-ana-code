#!/usr/bin/bash

export pdir=`pwd`
export groovy=$COATJAVA"/bin/run-groovy"
data_path="/volatile/clas12/rg-a/production/recon/pass0/v15/"
dir=`ls $data_path`
for run in $dir
do
	run=${run:2}
	$groovy electron_pID.groovy $run `find $data_path"00"$run -name "*.hipo"`
done
