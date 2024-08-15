#!/bin/bash
freyja update
freyja update --outdir . #to ensure files end up in this dir
my_func() {
    fn=$1
    depthfolder=$2
    output=$3

    fn_out=${fn##*/}
    baseName=${fn##*/}
    baseName=${baseName%.*}
    depthfile0="${depthfolder}${fn_out%.tsv}.depths"
    output0="${output}${fn_out%.*}.demix.tsv"
    echo $fn
    echo $depthfile0
    echo $output0
    freyja demix $fn $depthfile0 --output $output0 --eps 0.0000001
    # if [ -f "$output0" ]; then
    #     echo "file already exists. skipping. "
    # else
    #     # echo "file doesn't exist. running. "
    #     freyja demix $fn $depthfile0 --output $output0 --eps 0.0000001
    # fi

}

export -f my_func
parallel -j 24 my_func ::: ../variants/* ::: ../depths/ ::: ../outputs/
