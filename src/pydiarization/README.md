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

TODO