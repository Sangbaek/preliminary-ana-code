#!/usr/bin/bash

export pdir=`pwd`
export groovy=$COATJAVA"/bin/run-groovy"
data_path=$pdir"/../"
dir=`ls $data_path`
$groovy elec_vz.groovy `find dir -name "*.hipo"`
