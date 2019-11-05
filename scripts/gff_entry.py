class CDS:
    def __init__(self, start, end, frame, status):
        self.start = start
        self.end = end
        self.frame = frame
        self.status = status


class Exon:
    def __init__(self, parent, id, start, end, status, strand, source, CDS):
        self.start = start
        self.end = end
        self.id = id
        self.source = source
        self.parent = parent
        self.status = status
        self.strand = strand
        self.CDS = CDS


class Transcript:
    def __init__(self, id, start, end, exons, strand, status, source):
        self.id = id
        self.start = start
        self.end = end
        self.Exons = exons
        self.strand = strand
        self.status = status
        self.source = source


class Gene:
    def __init__(self, chrm, name, start, end, transcripts, strand, status, source):
        self.chrm = chrm
        self.name = name
        self.start = start
        self.end = end
        self.Transcripts = transcripts
        self.strand = strand
        self.status = status
        self.source = source
