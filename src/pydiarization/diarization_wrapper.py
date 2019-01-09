"""
Module that wrappe the diarization toolkit

Author: Droz William <william.droz@idiap.ch>
"""
import tempfile
import os
import re
import subprocess
import pydiarization.config as config

def video_to_rttm_string(video_path, external_scp=None):
    """ High-level function that return the rttm as string from a video
    Arguments:
    video_path -- where the video is
    external_scp -- external scp file used to create sequences

    Return: the rttm content as string
    """
    rttm_path = tempfile.NamedTemporaryFile(suffix='.rttm', dir=config.p_tmp_folder).name
    rttm_from_video(video_path, rttm_path, external_scp)
    rttm_str = rttm_to_string(rttm_path)
    os.remove(rttm_path)
    return rttm_str

def wav_to_rttm_string(wav_path, external_scp=None):
    """ High-level function that return the rttm as string from a wav
    Arguments:
    wav_path -- where the wav is
    external_scp -- external scp file used to create sequences

    Return: the rttm content as string
    """
    rttm_path = tempfile.NamedTemporaryFile(suffix='.rttm', dir=config.p_tmp_folder).name
    rttm_from_wav(wav_path, rttm_path, external_scp)
    rttm_str = rttm_to_string(rttm_path)
    os.remove(rttm_path)
    return rttm_str

def rttm_to_string(rttm_path):
    """ transform rttm file to string
    Arguments:
    rttm_path -- where is the rttm file

    Return: rttm content as string
    """
    with open(rttm_path, 'r') as f:
        return f.read()

def rttm_from_video(video_path, rttm_path, external_scp=None):
    """ create a .rttm file from a video
    Arguments:
    video_path -- path to the video
    rttm_path -- path where the .rrtm file will be saved
    external_scp -- external scp file used to create sequences
    """
    wav_path = tempfile.NamedTemporaryFile(suffix='.wav', dir=config.p_tmp_folder).name
    _video_to_wav(video_path, wav_path)
    rttm_from_wav(wav_path, rttm_path, external_scp)
    os.remove(wav_path)

def rttm_from_wav(wav_path, rttm_path, external_scp=None):
    """ create a .rttm file from a wav
    Arguments:
    wav_path -- path to the wav file
    rttm_path -- path where the .rrtm file will be saved
    external_scp -- external scp file used to create sequences
    """
    fea_path = tempfile.NamedTemporaryFile(suffix='.fea', dir=config.p_tmp_folder).name
    kaldi_wav_scp_path = tempfile.NamedTemporaryFile(suffix='.scp', dir=config.p_tmp_folder).name
    kaldi_ark_scp_path = tempfile.NamedTemporaryFile(suffix='.scp', dir=config.p_tmp_folder).name
    ark_path = tempfile.NamedTemporaryFile(suffix='.ark', dir=config.p_tmp_folder).name
    scp_path = tempfile.NamedTemporaryFile(suffix='.scp', dir=config.p_tmp_folder).name
    tmp_rttm_file = tempfile.NamedTemporaryFile(suffix='.rttm', dir=config.p_result_folder).name

    fileid = _extract_name_without_ext(fea_path)

    _create_kaldi_wav_scp_file(wav_path, kaldi_wav_scp_path, fileid)
    _create_kaldi_ark(kaldi_wav_scp_path, ark_path)
    _create_kaldi_ark_scp_file(ark_path, kaldi_ark_scp_path, fileid)
    _convert_kaldi_mfcc_ark_to_htk_fea(config.p_tmp_folder, ark_path)
    if external_scp is not None:
        _create_scp_file_without_htk_from_external_scp(fea_path, scp_path, fileid, external_scp)
    else:
        _create_scp_file_without_htk(fea_path, scp_path, fileid)
    _create_rttm_file(fea_path, scp_path, rttm_path, tmp_rttm_file)

    os.remove(fea_path)
    os.remove(scp_path)
    os.remove(kaldi_wav_scp_path)
    os.remove(kaldi_ark_scp_path)
    os.remove(ark_path)

def _get_sequences(external_scp):
    """ retreive all the sequences from a scp file
    Arguments:
    external_scp -- external scp file used to create sequences
    """
    with open(external_scp, 'r') as f:
        return re.findall(r'\.fea\[(\d+),\s*(\d+)\]', f.read())

def _create_scp_file_without_htk_from_external_scp(fea_path, scp_path, fileid, external_scp):
    """ create a scp file based on a existing one (external)
    Arguments:
    fea_path -- path to the fea file
    scp_path -- path where the scp file will be saved
    fileid -- identifiant of file rttm
    external_scp -- external scp file used to create sequences
    """
    fea_nb_bytes = _get_nb_bytes_from_file(fea_path)
    nb_ceps = _get_num_ceps()
    nb_frames = int(((fea_nb_bytes - 14) / 4) / nb_ceps)
    with open(scp_path, 'w') as f:
        for seq in _get_sequences(external_scp):
            f.write('{}_{}_{}={}[{},{}]\n'.format(fileid, *seq, fea_path, *seq))


def _create_kaldi_wav_scp_file(wav_path, scp_path, fileid):
    """ create a scp file for getting the .wav for kaldi
    Arguments:
    wav_path -- path to the wav file
    scp_path -- path where the .scp file will be created
    fileid -- the unique name for identify the file
    """

    scp_content = '{} {}'.format(fileid, wav_path)
    with open(scp_path, 'w') as f:
        f.write(scp_content)


def _create_kaldi_ark(scp_path, ark_path):
    """ create a kaldi ark file from the scp wav file
        and using a config file. (that compute mfcc)
    Arguments:
    scp_path -- path to the scp wav file
    ark_path -- path where the ark file will be written
    """
    command = '{} --config={} scp:{} ark:{}'.format(config.p_binary_compute_mfcc_feats,
                                                    config.p_kaldi_mfcc_conf_file, scp_path, ark_path)
    _run_cmd(command, 'compute-mfcc-feats')

def _video_to_wav(video_path, wav_path):
    """ extract a wav from a video using ffmpeg
    Arguments:
    video_path -- path to the video
    wav_path -- path where the .wav file will be saved
    """
    command = '{} -i {} -y -f wav -flags bitexact -ar {} -ac 1 -acodec pcm_s16le {}'.format(
        config.p_binary_ffmpeg,
        video_path,
        config.p_default_sampling,
        wav_path
    )
    _run_cmd(command, 'ffmpeg')

def _create_kaldi_ark_scp_file(ark_path, scp_path, fileid):
    """ create a scp file to access the ark file from kaldi
    Arguments:
    ark_path -- the path to the ark file
    scp_path -- where the scp file will be written
    fileid -- identifiant of the file
    """
    scp_content = '{} {}'.format(fileid, ark_path)
    with open(scp_path, 'w') as f:
        f.write(scp_content)

def _convert_kaldi_mfcc_ark_to_htk_fea(fea_folder_path, ark_path):
    """ convert a kaldi mfcc file in ark format to htk fea file
    Arguments:
    fea_folder_path -- the folder where the .fea file will be created, 
                       the named inside the .scp file will by used.
    ark_path -- the path of the ark file
    """
    command = '{} --output-dir={} --output-ext=fea ark:{}'.format(config.p_binary_copy_feats_to_htk, fea_folder_path, ark_path)
    _run_cmd(command, 'copy-feats-to-htk')

def _get_parameters_from_mfcc_conf(parameter, param_value_regexp):
    """ Retreive a parameter from the kaldi mfcc config file
    Arguments:
    parameter -- the parameter to look at
    param_value_regexp -- the sub part of the regex to get the value of the parameter
    """
    expression = '{}\\s*=\\s*({})'.format(parameter, param_value_regexp)
    with open(config.p_kaldi_mfcc_conf_file, 'r') as f:
        res = re.search(expression, f.read())
        return res.group(1)

def _get_num_ceps():
    """ get the num_ceps parameter from the mfcc conf file """
    return int(_get_parameters_from_mfcc_conf('--num-ceps', "\\d+"))

def _get_nb_bytes_from_file(file_path):
    """ return the number of bytes of a file """
    with open(file_path, 'rb') as f:
        return len(f.read())

def _create_scp_file_without_htk(fea_path, scp_path, fileid):
    """ create a dummy scp file from a .fea file using pure python
    Arguments:
    fea_path -- path to the fea file
    scp_path -- path where the scp file will be saved
    fileid -- identifiant of file rttm    
    """
    fea_nb_bytes = _get_nb_bytes_from_file(fea_path)
    nb_ceps = _get_num_ceps()
    nb_frames = int(((fea_nb_bytes - 14) / 4) / nb_ceps)
    with open(scp_path, 'w') as f:
        f.write('{}_0_{}={}[0,{}]\n'.format(fileid, nb_frames, fea_path, nb_frames))

def _create_rttm_file(fea_path, scp_path, rttm_path, tmp_rttm_file):
    """ Create a .rttm file by calling the diarization toolkit
    Arguments:
    fea_path -- the path to the .fea file
    scp_path -- the path to the .scp file
    rttm_path -- the path where the .rttm will be saved
    tmp_rttm_file -- temporally .rttm file
    """
    fileid = _extract_name_without_ext(tmp_rttm_file)
    commands = '\n'.join([
        '{} --mfcc {} 1.0 --recid {} --outdir {} --tmpdir {} -s {} --beta {} --nthread {}'.format(
            config.p_binary_diarizeme,
            fea_path,
            fileid,
            config.p_result_folder,
            config.p_tmp_folder,
            scp_path,
            config.p_beta_val,
            config.p_nb_threads
        ),
        'mv {} {}'.format(tmp_rttm_file, rttm_path)
    ])
    _run_cmd(commands, 'create rttm file')
    return fileid


def _run_cmd(command, name, debug=False):
    """ run command using subprocess
    Arguments:
    command -- command to run
    name -- name of the command for printing error in case of return code != 0
    debug -- if true, print stdout and stderr
    """
    process = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate(command.encode('utf-8'))
    code = process.returncode
    if code != 0:
        print('error with {}: {}'.format(name, err.decode('utf-8')))
    if debug:
        print(out.decode('utf-8'))
        print(err.decode('utf-8'))
    return out, err

def _extract_name_without_ext(file_path):
    """ extract the name of a file without the extension
    Arguments:
    file_path -- file to extract the name
    Return -- the name without the extension
    """
    basename = os.path.basename(file_path)
    file_without_ext, _ = os.path.splitext(basename)
    return file_without_ext

def _debug_cat_file(file_path):
    """ put the content of file_path in stdout
    Arguments:
    file_path -- file to read
    """
    with open(file_path, 'r') as f:
        print(f.read())
