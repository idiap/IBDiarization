"""
module that test diarization wrapper (not unit testing)

Author: Droz William <william.droz@idiap.ch>
"""
import os
import subprocess
from pydiarization.diarization_wrapper import video_to_rttm_string, wav_to_rttm_string
import pydiarization.config as config

def _run_cmd(command):
    """ run command using subprocess
    Arguments:
    command -- command to run
    """
    process = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.communicate(command.encode('utf-8'))
    code = process.returncode
    return code == 0

def test_binaries():
    """ test all the binaries """
    binary_parameters = {config.p_binary_compute_mfcc_feats: '--help',
                         config.p_binary_copy_feats_to_htk: '--help',
                         config.p_binary_diarizeme: '--help',
                         config.p_binary_ffmpeg: '--help'}
    is_ok = True
    for binary, parameters in binary_parameters.items():
        if not test_binary(binary, parameters):
            is_ok = False
            print('something wrong with {}'.format(binary))
    return is_ok

def test_binary(binary, parameters):
    """ test a specific binary with specific parameters """
    command = '{} {}'.format(binary, parameters)
    return _run_cmd(command)
    
def test_end_to_end_audio_with_scp_file(audio_file, external_scp, reference_rttm):
    """ test"""
    res = wav_to_rttm_string(audio_file, external_scp)
    print(res)
    return True

def test_end_to_end(video_file):
    """ test the Diarization Toolkit with the video_file
    Arguments:
    video_file -- video file from either in path or URL
    """
    res = video_to_rttm_string(video_file)
    if len(res) > 10:
        return True
    else:
        print('end-to-end system seem broken')
        return False

if __name__ == '__main__':
    audio_file = 'AMI_20041210-1052_seg.wav'
    external_scp = 'AMI_20041210-1052.scp'
    reference_rttm = 'AMI_20041210-1052.rttm'
    is_ok = True
    is_ok &= test_binaries()
    is_ok &= test_end_to_end_audio_with_scp_file(audio_file, external_scp, reference_rttm)
    if is_ok:
        print('diarization_wrapper -> OK')
    else:
        print('diarization_wrapper -> failed')
    exit(int(not is_ok))
    
