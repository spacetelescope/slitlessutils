from pathlib import Path
from tempfile import TemporaryDirectory

from astroquery.mast import Observations

from slitlessutils.core.preprocess.crrej.drizzle import drizzle
from slitlessutils.examples.wr96 import WR96_DRIZZLE_PARAMS


def test_wr96_drizzle():
    '''
    Runs the drizzle step embedded in the wr96 example
    '''
    uris = {
        'jdql01jnq_flc.fits': 'mast:HST/product/jdql01jnq_flc.fits',
        'jdql01jpq_flc.fits': 'mast:HST/product/jdql01jpq_flc.fits',
        'jdql01jvq_flc.fits': 'mast:HST/product/jdql01jvq_flc.fits',
        'jdql01jxq_flc.fits': 'mast:HST/product/jdql01jxq_flc.fits'
    }
    # Setup temporary directory
    with TemporaryDirectory() as tempdir:
        rawdata_dir = Path(tempdir) / 'raw_data'
        rawdata_dir.mkdir(parents=True)
        mosaic_dir = Path(tempdir) / 'mosaic_files'
        mosaic_dir.mkdir(parents=True)
        # Download wr96 visit files
        for filename in uris:
            Observations.download_file(uris[filename], local_path=str(rawdata_dir / filename))
        rawdata_filepaths = [str(filepath) for filepath in rawdata_dir.iterdir()]

        # Check that our temp folder is indeed empty
        assert len(list(mosaic_dir.iterdir())) == 0
        # Actually perform drizzle
        drizzle(rawdata_filepaths, outdir=mosaic_dir, **WR96_DRIZZLE_PARAMS)
        # Confirm we have our output mosaics
        assert len(list(mosaic_dir.iterdir())) > 0

