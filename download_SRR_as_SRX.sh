#download SRR


try_down () {
    target=$1
    name=$2
    echo "downloading $target as $name"
    num123=`echo $target | cut -c 4-6`
    if [ -z $name ]; then
        wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR$num123/$target/$target'.sra'
    else
        wget -O $name'.sra' ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR$num123/$target/$target'.sra'
    fi
    echo "downloading  $target as $name done"
}


# list handle
date

if [ ! -z $3 ]; then
    echo "downloading list from $3, the first 2 parameters would be ignored"
    while read p; do
        try_down $p
    done < $3
else
    echo "downloading only 1 file"
    try_down $1 $2
fi


