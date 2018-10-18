This README file explains, in brief, the commands to 
compile and run the IB diarization toolkit.

# PREREQUISITES

The toolkit has been tested to run on linux environments.
There are three prerequisites to the package

1. Cmake -- required to build the package
2. libboost -- boost library for C++
3. openmp library -- for multithreaded processing

# CLEAN COMPILATION

```
$> cd src/diarization/cmake/
$> cmake .
$> # make clean if cmake has already been run once
$> make
```

# RUNNING THE TOOLKIT

```
$> # cd to $IB_DIARIZATION_HOME
$> bash scripts/run.diarizeme.sh ipfile scpfile opfolder fileid [betaval]
```


Example command

```
$> bash scripts/run.diarizeme.sh data/mfcc/AMI_20050204-1206.fea data/scp/AMI_20050204-1206.scp result.dir/ AMI_20050204-1206
```

The files required to run the script are as follows:
  - AMI_20050204-1206.fea: Feature file (e.g. MFCCs) in HTK format.
  - AMI_20050204-1206.scp: Segmentation file with each line marking a segment in the feature file.
                           Each line follows the format (without quotes) "segment_name=filename[start_frame,end_frame]".
  - result.dir: The directory in which the result file will be stored. All intermediate files will also be stored here.
  - AMI_20050204-1206: The name of the audio file to be used while storing result files. In the result folder, which is
                       the folder result.dir in the above command, will contain a result file called AMI_20050204-1206.rttm.

To test the result use md-perl-eval tool available on the NIST website

```
$> perl md-eval-v21.pl -m -afc -c 0.25 -r data/rttm/AMI_20050204-1206.rttm -s result.dir/AMI_20050204-1206.rttm
```

The expected DER is 8.79%

To add TDOA features (or other complementary features), use the following script

```
$>bash scripts/run.diarizeme.tdoa.sh data/mfcc/AMI_20050204-1206.fea 0.8 data/tdoa/AMI_20050204-1206.fea 0.2 data/scp/AMI_20050204-1206.scp result.dir.tdoa AMI_20050204-1206
```

The expected DER is 7.12%

# Using Diarization toolkit as a command

The diarization engine can also be accessed as a linux command. The relevant binary is src/diarization/cmake/diarizeme.
To run the command instead of the script use the following command

```
$> src/diarization/cmake/diarizeme \
    --mfcc data/mfcc/AMI_20050204-1206.fea 1.0 \
    --recid AMI_20050204-1206 \
    --outdir result.dir \
    --tmpdir result.dir \
    -s data/scp/AMI_20050204-1206.scp \
    --beta 10 \
    --nthread 1 
```

In the above command, the MFCC features are given a weight of 1.0. The tmpdir option
is also pointed to result.dir directory. However, it can be different from the outdir
option as well. 

The beta value, which is the Lagrangian parameter used during IB clustering, is optimized for the AMI corpus. 
This may have to be optimized for other datasets.

The nthread option sets the number of threads to be used for clustering.

To add other features along with MFCCs simply use the --other option. A weight along with the file name will
have to be supplied. The sum of weights for all features should sum to 1.0.
For example, to use TDOA features along with MFCC features run

```
$> src/diarization/cmake/diarizeme \
    --mfcc data/mfcc/AMI_20050204-1206.fea 0.8 \
    --tdoa data/tdoa/AMI_20050204-1206.fea 0.2 \
    --recid AMI_20050204-1206 \
    --outdir result.dir \
    --tmpdir result.dir \
    -s data/scp/AMI_20050204-1206.scp \
    --beta 10 \
    --nthread 1 
```
#  Pydiarization - a python wrapper

The [pydiarization](https://github.com/idiap/IBDiarization/tree/master/src/pydiarization) subproject provides high-level API around the Diarization toolkit.
