diarizeme_bin=/idiap/temp/msrikanth/repo/idiap-aib-diarization-dimseng/src/diarization/cmake/diarizeme
beta_val=10

numargs=3
if [ $# -lt $numargs ]; then
    echo "Usage: bash run.diarizeme.sh iplist scpfolder opfolder"
    exit
fi

iplist=$1
scpfolder=$2
opdir=$3

export nthreads=1

if [ -e $opdir -a ! -d $opdir ]; then
    echo "Warning: $opdir exists but is not a directory"
    exit
elif [ ! -d $opdir ]; then
    mkdir $opdir
fi    

for fname in `cat $iplist`; do
    fileid=`basename $fname .fea`
    scpfile=$scpfolder/$fileid.scp
echo "$diarizeme_bin \
    --mfcc $fname 1.0 \
    --recid $fileid \
    --outdir $opdir \
    --tmpdir $opdir \
    -s $scpfile \
    --beta $beta_val \
    --nthread 1 > $opdir/$fileid.out" > temp.$fileid.sh
qsub -l q1d temp.$fileid.sh     
done
