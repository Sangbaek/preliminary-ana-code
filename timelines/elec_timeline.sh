#!/usr/bin/bash

export pdir=`pwd`
export groovy=$COATJAVA"/bin/run-groovy"
data_path=$pdir"/../"
$groovy elec_vz.groovy `find $data_path -name "*.hipo"`
