""" Config file for the diarization toolkit wrapper

Author: Droz William <william.droz@idiap.ch>
"""
from os import environ, path

__module_folder=path.dirname(path.abspath(__file__))

# default sampling for creating wav file (in Hz)
p_default_sampling = 16000

# where to find the runnable binary of the diarization toolkit
p_binary_diarizeme = 'diarizeme'

# where to find ffmpeg
p_binary_ffmpeg = 'ffmpeg'

# where to find compute-mfcc-feats
p_binary_compute_mfcc_feats = 'compute-mfcc-feats'

# where to find copy-feats-to-htk
p_binary_copy_feats_to_htk = 'copy-feats-to-htk'

# where to find the kaldi config file to compute mfcc
p_kaldi_mfcc_conf_file = path.join(__module_folder, 'kaldi_mfcc.conf')

# tmp folder for diarization toolkit (must already exist)
p_tmp_folder = 'result.dir'

# folder that will contains results of the diarization toolkit (must already exist)
p_result_folder = 'result.dir'

# value of beta for diarization toolkit
p_beta_val = 10

# number of thread for the diarization toolkit
p_nb_threads = 4

def _update_parameters_from_env():
    """ update the config values from env
    """
    for k, v in globals().items():
        if k.startswith('p_'):
            if k in environ:
                new_v = type(v)(environ[k])
                print('Parameters {} was updated from {} to {} by environment override'.format(k, v, new_v))
                globals()[k] = new_v

# update the config by using env
_update_parameters_from_env()