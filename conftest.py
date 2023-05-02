from os import environ
from slitlessutils.config import download_config_files


if 'slitlessutils_config' not in environ:
    download_config_files()