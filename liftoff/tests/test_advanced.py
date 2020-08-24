from glob import glob
import liftoff
import liftoff.run_liftoff
import os
import pytest

def test_yeast(tmp_path):
    # cleanup
    for tfile in glob('liftoff/tests/*.mmi') + glob('liftoff/tests/*.fai') + glob('liftoff/tests/*_db'):
        try:
            os.unlink(tfile)
        except OSError:
            pass

    # inputs
    asmbl = 'liftoff/tests/GCA_000146045.2_R64_genomic.fna.gz'
    annot = 'liftoff/tests/GCA_000146045.2_R64_genomic.gff.gz'
    target = 'liftoff/tests/GCA_000146045.2_R64_genomic.fna.gz'
    chroms = 'liftoff/tests/chroms.txt'
    unplaced = 'liftoff/tests/unplaced.txt'

    # outputs
    output = str(tmp_path / 'GCA_000146045.2_R64_to_GCA_000146045.2_advanced.gff3')
    unmapped = str(tmp_path / 'unmapped_features_GCA_000146045.2_R64_advanced.txt')
    tempdir = str(tmp_path / 'sandbox')
    expout = 'liftoff/tests/GCA_000146045.2_R64_to_GCA_000146045.2_R64_expected_advanced.gff'
    expunmapped = 'liftoff/tests/unmapped_features_GCA_000146045.2_R64_expected_advanced.txt'

    # run the program
    args = ['-g', annot, '-o', output, '-u', unmapped, '-dir', tempdir, '-chroms',chroms, '-unplaced', unplaced, '-copies',  target, asmbl]
    liftoff.run_liftoff.main(args)

    # verify the output
    with open(output, 'r') as fh1, open(expout, 'r') as fh2:
        observed_output = fh1.read().strip()
        expected_output = fh2.read().strip()
        assert observed_output == expected_output
    with open(unmapped, 'r') as fh1, open(expunmapped, 'r') as fh2:
        observed_output = fh1.read().strip()
        expected_output = fh2.read().strip()
        assert observed_output == expected_output
