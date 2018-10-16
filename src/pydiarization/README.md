# pydiarization

pydiarization is a wrapper around the [IBDiarization](https://github.com/idiap/IBDiarization) toolkit.

## Requirements

You must have the follow binaries in your path:
*  **ffmpeg**
*  **compute-mfcc-feats** and **copy-feats-to-htk** from Kaldi
*  **diarizeme** from the IBDiarization toolkit.

## Installation

pydiarization can be installed by either conda
```
conda install -c wdroz pydiarization
```
or pip
```
pip install pydiarization
```

## Usage

Before using the pywrapper, you have to create a folder that will contains the results of the IBDiarization toolkit. 
```bash
mkdir result.dir
```

## Test the installation

To check if all binaries works and are recognized by pydiarization, you can run the tests by typing: ```python3 -m pydiarization.test_diarization_wrapper```

### from code

Example that convert a single video to a string (.rttm content)
```python
from pydiarization.diarization_wrapper import video_to_rttm_string

rttm_content = video_to_rttm_string('MY_VIDEO_OR_URL.avi')
```

Here the list of all the high-levels API:
```python
def video_to_rttm_string(video_path):
    """ High-level function that return the rttm as string from a video
    Arguments:
    video_path -- where the video is

    Return: the rttm content as string
    """
def wav_to_rttm_string(wav_path):
    """ High-level function that return the rttm as string from a wav
    Arguments:
    wav_path -- where the wav is

    Return: the rttm content as string
    """
def rttm_to_string(rttm_path):
    """ transform rttm file to string
    Arguments:
    rttm_path -- where is the rttm file

    Return: rttm content as string
    """
def rttm_from_video(video_path, rttm_path):
    """ create a .rttm file from a video
    Arguments:
    video_path -- path to the video
    rttm_path -- path where the .rrtm file will be saved
    """
def rttm_from_wav(wav_path, rttm_path):
    """ create a .rttm file from a wav
    Arguments:
    wav_path -- path to the wav file
    rttm_path -- path where the .rrtm file will be saved
    """

```

### from CLI

The usage is the follow:
<pre>
$ python3 -m pydiarization.run --help
usage: run.py [-h] [--video VIDEO] [--wav WAV] [--output OUTPUT]

Diarization Toolkit Wrapper CLI

optional arguments:
  -h, --help       show this help message and exit
  --video VIDEO    take a video file or URL as input
  --wav WAV        take a wav file as input
  --output OUTPUT  specify to write as .rttm file
</pre>

#### some examples

Ask to get the .rttm content from a url:
`python3 -m pydiarization.run --video http://data.cstr.inf.ed.ac.uk/summa/data/test.mp4`
<pre>SPEAKER tmp86hrhwsd 1 0.01 292.17 <NA> <NA> tmp86hrhwsd_spkr_9 <NA>
SPEAKER tmp86hrhwsd 1 292.18 2.98 <NA> <NA> tmp86hrhwsd_spkr_2 <NA></pre>

Ask to write to .rttm file instead:
`python3 -m pydiarization.run --video http://data.cstr.inf.ed.ac.uk/summa/data/test.mp4 --output toto.rttm`