import pytest
from slitlessutils import config

cfg = config.Config()

# download the latest reference files
reffile = cfg.retrieve_reffiles(update=True)

# Copied over from https://github.com/spacetelescope/ci_watson
@pytest.fixture(scope='function')
def _jail(tmp_path):
    """Perform test in a pristine temporary working directory."""
    old_dir = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield str(tmp_path)
    finally:
        os.chdir(old_dir)