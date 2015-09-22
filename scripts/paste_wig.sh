#! /usr/bin/sh

dir="./inst/wig_files"

set=$1
files=`ls $dir | grep $set`

out_file=$set.wig

fl=200

echo -n "" > $dir/$out_file

for file in $files; do
    ff=${file//$set/}
    ff=${ff//.wig/}
    ff=${ff//_/}
    echo 'variableStep chrom='$ff 'span='$fl >> $dir/$out_file
    cat $dir/$file >> $dir/$out_file
done
