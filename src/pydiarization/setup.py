from setuptools import setup, find_packages

long_description = ''
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='pydiarization',
      version='0.4',
      packages=find_packages(where='../pydiarization'),
      package_dir={'pydiarization': '.'},
      package_data={'pydiarization': ['kaldi_mfcc.conf']},
      author="William Droz",
      author_email="william.droz@idiap.ch",
      description="Diarization toolkit wrapper",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="https://github.com/wdroz/IBDiarization/tree/master/src/pydiarization",
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Topic :: Multimedia :: Sound/Audio :: Speech"
    ])
