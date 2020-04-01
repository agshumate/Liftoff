class exon_node:
    def __init__(self, chrm, lifted_start, lifted_end, strand, original_start, original_end, alignment):
        self.chrm = chrm
        self.lifted_start = lifted_start
        self.lifted_end = lifted_end
        self.strand = strand
        self.original_start = original_start
        self.original_end = original_end
        self.alignment = alignment
