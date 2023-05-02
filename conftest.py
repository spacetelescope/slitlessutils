from os import environ
from pathlib import Path
from tarfile import open as tar_open
from tempfile import gettempdir
from astropy.utils.data import download_file


def download_config_files(extract_path=gettempdir()):
    '''Method to automatically download the config files.'''

    # Cast to pathlib path object (in case we got a string)
    extract_path = Path(extract_path)

    # Download the tar file
    config_url = 'https://stsci.box.com/shared/static/fzlb7y36ofi18ziy6mkbyg710stmygjf.gz'
    config_archive = Path(download_file(config_url, cache=True))

    with tar_open(config_archive) as tar:
        # Only extract the files that are missing
        for file in tar:
            if not (extract_path / Path(file.name)).exists():
                tar.extract(file, path=extract_path)

    # Prepare environment required for slitlessutils to be imported
    environ['slitlessutils_config'] = str(extract_path / 'slitlessutils_config')


if 'slitlessutils_config' not in environ:
    download_config_files()