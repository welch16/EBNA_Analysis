#! /usr/bin/sh

dir="./inst/wig_files"

set=$1
fl=$2
files=`ls $dir | grep $set`

out_file=$set.wig

printf "" > $dir/$out_file
printf "track type=wiggle_0\n" >> $dir/$out_file

for file in $files; do
    ff=${file//$set/}
    ff=${ff//.wig/}
    ff=${ff//_/}
    printf 'variableStep chrom='$ff' span='$fl'\n'>> $dir/$out_file
    cat $dir/$file >> $dir/$out_file
done
