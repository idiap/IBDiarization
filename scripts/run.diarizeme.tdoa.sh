trap "exit" SIGINT

diarizeme_bin=src/diarization/cmake/diarizeme
beta_val=10

numargs=7
if [ $# -lt $numargs ]; then
    echo -ne "Usage: bash run.diarizeme.sh mfccfile mfccwt"
    echo  " tdoafile tdoawt scpfile opfolder fileid [betaval]"
    echo ""
    echo "mfccfile :  Input file with MFCC features in HTK format"
    echo "mfccwt   :  a value between 0.0 and 1.0 as a weight on MFCC features"
    echo "tdoafile : Input file with TDOA features in TDOA format"
    echo "tdoawt   : a value between 0.0 and 1.0 as a weight on TDOA features"
    echo "scpfile : file with speech boundaries in the following format"
    echo "          segment_name=file_name[start_frame,end_frame]"
    echo "outdir :  directory to store temporary files and output rttm in"
    echo "fileid :  ID of the file (to be used in rttm output)"
    echo "betaval   :  (optional) beta value for IB clustering"
    exit
fi

fname=$1
fwt=$2
tdoafname=$3
tdoawt=$4
scpfile=$5
opdir=$6
fileid=$7

if [ $# -eq 8 ]; then
    beta_val=$8
fi

export nthreads=1

if [ -e $opdir -a ! -d $opdir ]; then
    echo "Warning: $opdir exists but is not a directory"
    exit
elif [ ! -d $opdir ]; then
    mkdir $opdir
fi    

echo $diarizeme_bin \
    --mfcc $fname $fwt \
    --tdoa $tdoafname $tdoawt \
    --recid $fname \
    --outdir $opdir \
    --tmpdir $opdir \
    -s $scpfile \
    --beta $beta_val \
    --nthread 1 


$diarizeme_bin \
    --mfcc $fname $fwt \
    --tdoa $tdoafname $tdoawt \
    --recid $fileid \
    --outdir $opdir \
    --tmpdir $opdir \
    -s $scpfile \
    --beta $beta_val \
    --nthread 1  > $opdir/$fileid.out
