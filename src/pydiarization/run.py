"""
module for CLI with the diarization toolkit wrapper
"""
import argparse
import sys
import pydiarization.diarization_wrapper as diarization_wrapper

def handle_video(video, output):
    if not output:
        return diarization_wrapper.video_to_rttm_string(video)
    else:
        diarization_wrapper.rttm_from_video(video, output)

def handle_audio(audio, output):
    if not output:
        return diarization_wrapper.wav_to_rttm_string(audio)
    else:
        diarization_wrapper.rttm_from_wav(audio, output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Diarization Toolkit Wrapper CLI')
    parser.add_argument('--video', help='take a video file or URL as input')
    parser.add_argument('--wav', help='take a wav file as input')
    parser.add_argument('--output', help='specify to write as .rttm file')
    args = parser.parse_args()
    output = args.output
    if args.video:
        res = handle_video(args.video, output)
    elif args.wav:
        res = handle_audio(args.wav, output)
    else:
        print('error, you have to specify either --video or --wav', file=sys.stderr)
        exit(-1)
    if res is not None:
        print(res)