class aligned_seg:
    def __init__(self, aln_id, query_name, reference_name, query_block_start, query_block_end, reference_block_start,
                 reference_block_end, is_reverse, mismatches):
        self.aln_id = aln_id
        self.query_name = query_name
        self.reference_name = reference_name
        self.query_block_start = query_block_start
        self.query_block_end = query_block_end
        self.is_reverse = is_reverse
        self.reference_block_start = reference_block_start
        self.reference_block_end = reference_block_end
        self.mismatches = mismatches
