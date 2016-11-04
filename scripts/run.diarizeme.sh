trap "exit" SIGINT


diarizeme_bin=src/diarization/cmake/diarizeme
beta_val=10

numargs=4
if [ $# -lt $numargs ]; then
    echo "Usage: bash run.diarizeme.sh ipfile scpfile opfolder fileid [betaval]"
    echo ""
    echo "ipfile :  Input file with MFCC features in HTK format"
    echo "scpfile : file with speech boundaries in the following format"
    echo "          segment_name=file_name[start_frame,end_frame]"
    echo "outdir :  directory to store temporary files and output rttm in"
    echo "fileid :  ID of the file (to be used in rttm output)"
    echo "betaval   :  (optional) beta value for IB clustering"
    exit
fi

fname=$1
scpfile=$2
opdir=$3
fileid=$4

if [ $# -eq 5 ]; then
    beta_val=$5
fi

if [ -z $nthreads ]; then
    export nthreads=1
fi

if [ -e $opdir -a ! -d $opdir ]; then
    echo "Warning: $opdir exists but is not a directory"
    exit
elif [ ! -d $opdir ]; then
    mkdir -p $opdir
fi    

$diarizeme_bin \
    --mfcc $fname 1.0 \
    --recid $fileid \
    --outdir $opdir \
    --tmpdir $opdir \
    -s $scpfile \
    --beta $beta_val \
    --nthread $nthreads > $opdir/$fileid.out
