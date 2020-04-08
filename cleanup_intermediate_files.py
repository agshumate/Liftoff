import os
import glob


def cleanup_intermediate_files():
    for file in glob.glob("*_genes.fa"):
        os.remove(file)


    for file in glob.glob("*_blastdb*"):
        os.remove(file)


    for file in glob.glob("*.xml"):
        os.remove(file)