import process_blast_alignments as pba
import networkx as nx
import lifted_exon
import numpy as np



def find_best_mapping(coordinate_map, alignment_scores, gene_db, gene, exon_alignments, coords_to_exclude, weight_threshold):
    unique_exons = pba.find_unique_exons(gene_db, gene)
    node_dict, exon_graph = intialize_graph()
    mapped_exons = {}
    for i in range(len(unique_exons)):
        exon = unique_exons[i]
        exon_coords, exon_alignment_scores = convert_coords(coordinate_map, alignment_scores,  gene, exon)
        add_alignments_to_graph(exon_coords, exon_alignment_scores, i,  node_dict, exon_graph, gene, exon, exon_alignments, unique_exons, coords_to_exclude)
    prune_graph(exon_graph, node_dict)
    target_node, num_end_nodes = add_target_node(node_dict, exon_graph, len(unique_exons), gene, unique_exons)
    shortest_path = nx.shortest_path(exon_graph, source="0.0", target = target_node, weight='cost')
    shortest_path_weight = nx.shortest_path_length(exon_graph, source="0.0", target = target_node, weight='cost')
    if weight_threshold is None:
        add_mapping = True
    elif shortest_path_weight <= weight_threshold:
        add_mapping=True
    else:
        add_mapping = False
    if add_mapping:
        for node_name in shortest_path:
            node = node_dict[node_name]
            mapped_exons[str(node.original_start) + ":" + str(node.original_end)] = node
    return mapped_exons, shortest_path_weight


def same_alignment(node_successors, aln_num):
    for successor in node_successors:
        if successor.split(".")[1] == aln_num:
            return successor
    return None


def prune_graph(exon_graph, node_dict):
    for node in node_dict:
        if node != "0.0":
            aln_num = node.split(".")[1]
            node_successors = list(exon_graph.successors(node))
            same_aln_node = same_alignment(node_successors, aln_num)
            if same_aln_node is not None:
                for successor in node_successors:
                    if successor != same_aln_node:
                        exon_graph.remove_edge(node, successor)


def intialize_graph():
    exon_graph = nx.DiGraph()
    node_dict = {}
    node_dict["0.0"] = lifted_exon.exon_node("start", -1, -1, None, -1, -1, [])
    exon_graph.add_node("0.0")
    return node_dict, exon_graph


def convert_coords(coordinate_map, alignment_scores,  gene, exon):
    gene_start, gene_end = gene.start, gene.end
    relative_start = exon.start - gene_start
    relative_end = exon.end - gene_start
    if gene.strand == '+':
        return coordinate_map[relative_start: relative_end + 1], alignment_scores[relative_start: relative_end + 1]
    else:
        return coordinate_map[(gene_end-gene_start) - relative_end: (gene_end-gene_start) - relative_start + 1], \
               alignment_scores[(gene_end-gene_start) - relative_end: (gene_end-gene_start) - relative_start + 1]


def find_score_coords(exon, gene_start, gene_db):
    transcript = list(gene_db.parents(exon, level=1))
    all_cds = list(gene_db.children(transcript[0], featuretype = "CDS"))
    cds_list = pba.find_cds(exon, all_cds)
    if pba.exon_has_cds(cds_list):
        score_start = cds_list[0].start - gene_start
        score_end = cds_list[-1].end - gene_start
    else:
        score_start = exon.start - gene_start
        score_end = exon.end - gene_start
    return score_start, score_end


def add_alignments_to_graph(all_exon_coords, all_exon_alignment_scores, exon_num,  node_dict, exon_graph, gene, exon, exon_alignments, unique_exons, coords_to_exclude):
    num_alignments = np.shape(all_exon_coords)[1]
    for i in range(num_alignments):
        exon_coords = all_exon_coords[:,i]
        exon_scores = all_exon_alignment_scores[:,i]
        hsp_and_chrm = exon_alignments[i][:]
        node_key, is_valid_mapping, exon_scores = add_exon_node(exon_num +1, i, node_dict, exon_graph, exon_coords, hsp_and_chrm, exon,  coords_to_exclude, exon_scores)
        if is_valid_mapping:
            add_edges(gene, exon_scores, exon_graph, node_dict, exon_num +1, node_key,unique_exons, coords_to_exclude)



def add_exon_node(exon_num, alignment_num, node_dict, exon_graph, exon_coords, hsp_and_chrm, exon,  coords_to_exclude, exon_scores):
    strand = get_hsp_strand(hsp_and_chrm[0])
    target_chrm = hsp_and_chrm[1]
    node_key = str(exon_num) + "." + str(alignment_num)
    for criteria in coords_to_exclude:
        if strand == criteria[3] and target_chrm == criteria[2]:
            exon_coords, exon_scores = exclude_coords(exon_coords, criteria, exon_scores)
    lifted_exon_start, lifted_exon_end = find_boundary_coords(exon_coords)
    is_valid_mapping = True
    if lifted_exon_start !=0 and lifted_exon_end != 0:
        for criteria in coords_to_exclude:
            if strand == criteria[3] and target_chrm == criteria[2]:
                if pba.count_overlap(lifted_exon_start, lifted_exon_end, criteria[0], criteria[1]) > 0:
                    is_valid_mapping = False
    else:
        is_valid_mapping = False
    if is_valid_mapping:
        new_exon_node = lifted_exon.exon_node(target_chrm, lifted_exon_start, lifted_exon_end, strand, exon.start, exon.end, exon_coords)
        node_dict[node_key] = new_exon_node
        exon_graph.add_node(node_key)
    return node_key, is_valid_mapping, exon_scores


def exclude_coords(exon_coords, coords_to_exclude,  exon_alignment_scores):
    exclude_start, exclude_end = coords_to_exclude[0] -1  , coords_to_exclude[1] -1
    for i in range(len(exon_coords)):
        if exon_coords[i] >= exclude_start and exon_coords[i] <= exclude_end:
            exon_coords[i] = 0
            exon_alignment_scores[i] = 2
    return exon_coords, exon_alignment_scores


def get_hsp_strand(hsp):
    if hsp.sbjct_start > hsp.sbjct_end:
        strand = '-'
    else:
        strand = '+'
    return strand


def find_boundary_coords(exon_coords):
    if len(exon_coords) != 0:
        aligned_bases = exon_coords[exon_coords != 0]
        if len(aligned_bases) !=0:
            start, end = min(aligned_bases[0], aligned_bases[-1]), max(aligned_bases[0], aligned_bases[-1])
        else:
            start, end = 0,0
    else:
        start, end = 0,0
    return min(start, end), max(start, end)


def add_edges(gene, exon_scores,  exon_graph, node_dict, exon_num, node_key,  unique_exons, coords_to_exclude):
    mismatch_penalty = count_mismatches(exon_scores)
    adjacent_nodes = find_adjacent_nodes(node_dict, exon_num, node_key, exon_graph, gene, coords_to_exclude)
    for adjacent_node_name in adjacent_nodes:
        order_penalty = find_skipped_exons_penalty(adjacent_node_name, node_key, unique_exons)
        edge_weight = order_penalty + mismatch_penalty
        exon_graph.add_edge(adjacent_node_name, node_key, cost=edge_weight)


def count_mismatches(exon_scores):
    return np.sum(exon_scores)


def find_adjacent_nodes(node_dict, exon_num, node_key, exon_graph, gene, coords_to_exclude):
    adjacent_nodes = []
    current_node = node_dict[node_key]
    for previous_node_name in node_dict:
        previous_node = node_dict[previous_node_name]
        if is_adjacent(previous_node, current_node, exon_num, previous_node_name, node_key, exon_graph, gene, coords_to_exclude):
            adjacent_nodes.append(previous_node_name)
    if len(adjacent_nodes) == 0:
        adjacent_nodes.append("0.0")
    return adjacent_nodes



def is_adjacent(previous_node, current_node, exon_num, previous_node_name, current_node_name, exon_graph, gene, coords_to_exclude):
    if exon_num == 1 and previous_node.chrm == "start":
        return True
    elif int(previous_node_name.split(".")[0]) == exon_num:
        return False
    elif find_successors(exon_graph, current_node_name, previous_node_name) is True:
        return False
    elif current_node.chrm == "end":
        return True
    elif previous_node.strand != current_node.strand:
        return False
    elif previous_node.chrm != current_node.chrm:
        return False
    else:
        new_distance, original_distance = calculate_distance(previous_node, current_node)
    if np.abs(new_distance) < 0.5*(np.abs(original_distance)) or np.abs(new_distance) > 1.5 * (np.abs(original_distance)):
        return False
    elif pba.count_overlap(previous_node.original_start, previous_node.original_end, current_node.original_start,
                               current_node.original_end) <= 0 and pba.count_overlap(previous_node.lifted_start,
                                                                                     previous_node.lifted_end,
                                                                                     current_node.lifted_start,
                                                                                     current_node.lifted_end) >0:
        return False
    elif pba.count_overlap(previous_node.original_start, previous_node.original_end, current_node.original_start,
                               current_node.original_end) > 0 and pba.count_overlap(previous_node.lifted_start,
                                                                                     previous_node.lifted_end,
                                                                                     current_node.lifted_start,
                                                                                     current_node.lifted_end) <= 0:
        return False
    elif intron_exclusion(coords_to_exclude, previous_node, current_node):
        return False
    elif gene.strand == previous_node.strand and previous_node.lifted_start <= current_node.lifted_start:
        return True
    elif gene.strand != previous_node.strand and previous_node.lifted_end >= current_node.lifted_end:
        return True
    else:
        return False


def intron_exclusion(coords_to_exclude, previous_node, current_node):
    for exclusion_coords in coords_to_exclude:
        min_coord = min(previous_node.lifted_start, previous_node.lifted_end, current_node.lifted_start, current_node.lifted_end)
        max_coord = max(previous_node.lifted_start, previous_node.lifted_end, current_node.lifted_start, current_node.lifted_end)
        if pba.count_overlap(min_coord, max_coord, exclusion_coords[0], exclusion_coords[1]) >0:
            return True

    return False



def find_successors(exon_graph, current_node_name, previous_node_name):
    exon_key = current_node_name.split(".")[0]
    successors = list(exon_graph.successors(previous_node_name))
    if len(successors) == 0:
        return False
    else:
        for successor in successors:
            if successor.split(".")[0] != exon_key:
                return True
    return False

def calculate_distance(previous_node, current_node):
    max_original = max(previous_node.original_start, previous_node.original_end, current_node.original_start, current_node.original_end)
    min_original = min(previous_node.original_start, previous_node.original_end, current_node.original_start, current_node.original_end)
    original_distance = max_original - min_original
    max_lifted = max(previous_node.lifted_start, previous_node.lifted_end, current_node.lifted_start, current_node.lifted_end)
    min_lifted = min(previous_node.lifted_start, previous_node.lifted_end, current_node.lifted_start, current_node.lifted_end)
    new_distance = max_lifted - min_lifted
    return original_distance, new_distance


def find_skipped_exons_penalty(previous_node_name, current_node_name, unique_exons):
    exon_number1 = int(previous_node_name.split(".")[0])
    exon_number2 = int(current_node_name.split(".")[0])
    penalty = 0
    if exon_number2 - exon_number1 == 1:
        return 0
    else:
        for i in range(int(exon_number1), int(exon_number2)-1):
            penalty += (unique_exons[i].end - unique_exons[i].start) * 5
    return penalty


def add_target_node(node_dict, exon_graph, num_exons, gene, unique_exons):
    node_key = str(num_exons +1) + ".0"
    node_dict[node_key] = lifted_exon.exon_node("end", -1, -1, None, -1, -1, [])
    exon_graph.add_node(node_key)
    adjacent_nodes = find_adjacent_nodes(node_dict, num_exons + 1 , node_key, exon_graph, gene, [])
    for node in adjacent_nodes:
        cost = find_skipped_exons_penalty(node, node_key, unique_exons)
        exon_graph.add_edge(node, node_key, cost=cost)
    return node_key, len(adjacent_nodes)
