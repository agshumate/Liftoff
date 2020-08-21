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
    asmbl = 'liftoff/tests/GCF_000146045.2_R64_genomic.fna.gz'
    annot = 'liftoff/tests/GCF_000146045.2_R64_genomic_subset.gff.gz'
    target = 'liftoff/tests/GCF_001298625.1_SEUB3.0_genomic.fna.gz'

    # outputs
    output = str(tmp_path / 'lifted.gff3')
    unmapped = str(tmp_path / 'unmapped.txt')
    tempdir = str(tmp_path / 'sandbox')
    expout = 'liftoff/tests/Scer-to-Seub_liftover.gff3'
    expunmapped = 'liftoff/tests/Scer-to-Seub_unmapped.txt'

    # run the program
    args = ['-g', annot, '-o', output, '-u', unmapped, '-dir', tempdir, target, asmbl]
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
