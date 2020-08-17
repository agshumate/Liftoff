import networkx as nx
from liftoff import aligned_seg, liftoff_utils, new_feature
import numpy as np


def find_best_mapping(alignments, query_length, parent, feature_heirarchy, previous_gene_start, previous_gene_seq,
                      inter,
                      lifted_features_list, distance_factor, allow_chrom_split):
    children = feature_heirarchy.children[parent.id]
    children_coords = liftoff_utils.merge_children_intervals(children)
    node_dict, aln_graph = intialize_graph()
    head_nodes, has_full_length = add_single_alignments(node_dict, aln_graph, alignments, children_coords, parent,
                                       previous_gene_start, previous_gene_seq, inter, feature_heirarchy.parents,
                                                        lifted_features_list)

    chain_alignments(head_nodes, node_dict, aln_graph, parent, children_coords, inter,
                     feature_heirarchy.parents, lifted_features_list, distance_factor, allow_chrom_split,
                     has_full_length)
    add_target_node(aln_graph, node_dict, query_length, children_coords, parent)
    shortest_path_nodes = find_shortest_path(node_dict, aln_graph)



    shortest_path_nodes_by_chrm = split_shortest_path_nodes(shortest_path_nodes)
    all_converted_results = []
    if len(shortest_path_nodes) == 0:
        return [[{}, 0,0]]
    else:
        for node_group in shortest_path_nodes_by_chrm:
            mapped_children, alignment_coverage, seq_id = convert_all_children_coords(node_group, children, parent)
            all_converted_results.append([mapped_children, alignment_coverage, seq_id])
    return all_converted_results


def intialize_graph():
    aln_graph = nx.DiGraph()
    node_dict = {}
    node_dict[0] = aligned_seg.aligned_seg("start", "start", "start", -1, -1, -1, -1, -1, [], False)
    aln_graph.add_node(0)
    return node_dict, aln_graph


def add_single_alignments(node_dict, aln_graph, alignments, children_coords, parent,
                          previous_gene_start, previous_gene_seq, inter, parent_dict, lifted_features_list):

    alignments = sort_alignments(previous_gene_start, previous_gene_seq, alignments)
    node_num = 1
    previous_node_id = -1
    head_nodes = []
    previous_node = 0
    has_full_length = False
    for original_aln in alignments:
        if original_aln.full_length:
            has_full_length = True
        if inter is not None:
            aln = trim_overlap_coords(original_aln, parent, original_aln.query_name, inter,
                                      parent_dict, lifted_features_list)
        else:
            aln = original_aln
        if aln.aln_id != previous_node_id:
            previous_node = 0
            previous_node_id = aln.aln_id
        if is_valid_alignment(previous_node, aln, node_dict, parent, inter, parent_dict,
                              lifted_features_list):
            if previous_node == 0:
                head_nodes.append(node_num)
            current_node = add_to_graph(aln_graph, node_num, children_coords, parent, node_dict, previous_node, aln)
            previous_node = current_node
            node_num += 1
    return head_nodes, has_full_length


def sort_alignments(previous_gene_start, previous_gene_seq, alignments):
    order_dict = {}
    order = 0
    alignments = sorted(alignments, key=lambda x: [is_different_chrom(previous_gene_seq, x.reference_name), np.abs(
        previous_gene_start - x.reference_block_start)])
    for aln in alignments:
        if aln.aln_id not in order_dict:
            order_dict[aln.aln_id] = order
            order += 1
    alignments = sorted(alignments, key=lambda x: (order_dict[x.aln_id], x.query_block_start))
    return alignments


def is_different_chrom(previous_gene_seq, current_seq):
    return previous_gene_seq != current_seq

def trim_overlap_coords(aln, parent, feature_name, intervals, parent_dict, lifted_features_list):
    ref_start, ref_end = aln.reference_block_start, aln.reference_block_end
    target_chrm = aln.reference_name
    target_strand = get_strand(aln, parent)
    overlaps = liftoff_utils.find_overlaps(ref_start, ref_end, target_chrm, target_strand, feature_name, intervals,
                                           parent_dict, lifted_features_list)
    for overlap in overlaps:
        start, end, chrm, strand = overlap[0], overlap[1], overlap[2][1].seqid, overlap[2][
            1].strand
        if target_strand == strand and target_chrm == chrm:
            if start <= ref_start and end >= ref_end:
                ref_start, ref_end = -1, -1
                break
            elif end < ref_start or start > ref_end:
                ref_start, ref_end = ref_start, ref_end
            elif start >= ref_start and end <= ref_end:
                if start - ref_start > ref_end - end:
                    ref_start, ref_end = ref_start, start - 1
                else:
                    ref_start, ref_end = end + 1, ref_end
            elif start > ref_start:
                ref_start, ref_end = ref_start, start - 1
            else:
                ref_start, ref_end = end + 1, ref_end
    updated_aln = update_coords(ref_start, ref_end, aln)
    return updated_aln


def get_strand(aln, parent):
    if aln.is_reverse:
        if parent.strand == "-":
            strand = '+'
        else:
            strand = "-"
    else:
        strand = parent.strand
    return strand


def update_coords(ref_start, ref_end, aln):
    start_diff = ref_start - aln.reference_block_start
    end_diff = aln.reference_block_end - ref_end
    if ref_start != aln.reference_block_start or ref_end != aln.reference_block_end:
        new_query_start = aln.query_block_start + start_diff
        new_query_end = aln.query_block_end - end_diff
        new_aln = aligned_seg.aligned_seg(aln.aln_id, aln.query_name, aln.reference_name, new_query_start,
                                          new_query_end, ref_start, ref_end, aln.is_reverse, aln.mismatches,
                                          aln.full_length)
    else:
        new_aln = aln
    return new_aln


def is_valid_alignment(previous_node, aln, node_dict, parent, inter, parent_dict,
                       lifted_features_list):
    if aln.reference_block_start == -1:
        return False
    elif previous_node != 0:
        if spans_overlap_region(node_dict[previous_node], aln, parent, inter, parent_dict, lifted_features_list):
            return False
    return True


def spans_overlap_region(from_node, to_node, parent, intervals, parent_dict, lifted_features_list):
    if from_node.reference_name != to_node.reference_name:
        return False
    if intervals == None:
        return False
    target_chrm = from_node.reference_name
    node_overlap = get_node_overlap(from_node, to_node)
    strand = get_strand(from_node, parent)
    overlaps = liftoff_utils.find_overlaps(from_node.reference_block_end, to_node.reference_block_start + node_overlap,
                                           target_chrm, strand, from_node.query_name, intervals,
                                           parent_dict, lifted_features_list)
    return len(overlaps) > 0


def get_node_overlap(from_node, to_node):
    if from_node.reference_name == "start" or to_node.reference_name == "start":
        return 0
    else:
        return max(0, from_node.query_block_end - to_node.query_block_start +1)

def add_to_graph(aln_graph, node_num, children_coords, parent, node_dict, previous_node, aln):
    node_dict[node_num] = aln
    current_node = node_num
    node_weight = get_node_weight(aln, children_coords, parent)
    aln_graph.add_node(node_num, weight=node_weight)
    edge_weight = get_edge_weight(node_dict[previous_node], node_dict[current_node], children_coords, parent)
    aln_graph.add_edge(previous_node, current_node, cost=edge_weight)
    return current_node


def get_node_weight(aln, children_coords, parent):
    weight = 0
    for child_interval in children_coords:
        relative_start = liftoff_utils.get_relative_child_coord(parent, child_interval[0], aln.is_reverse)
        relative_end = liftoff_utils.get_relative_child_coord(parent, child_interval[1], aln.is_reverse)
        child_start, child_end = min(relative_start, relative_end), max(relative_start, relative_end)
        weight += len(aln.mismatches[(aln.mismatches >= child_start) & (aln.mismatches <= child_end)])
    return weight


def get_edge_weight(from_node, to_node, children_coords, parent):
    node_overlap = get_node_overlap(from_node, to_node)
    unaligned_range = [from_node.query_block_end + 1, to_node.query_block_start + node_overlap - 1]
    unaligned_exon_bases = 0
    for child_interval in children_coords:
        if from_node.reference_name == "start":
            is_reverse = to_node.is_reverse
        else:
            is_reverse = from_node.is_reverse
        relative_start = liftoff_utils.get_relative_child_coord(parent, child_interval[0], is_reverse)
        relative_end = liftoff_utils.get_relative_child_coord(parent, child_interval[1], is_reverse)
        child_start, child_end = min(relative_start, relative_end), max(relative_start, relative_end)
        overlap = liftoff_utils.count_overlap(child_start, child_end, min(unaligned_range[0], unaligned_range[1]), max(unaligned_range[0], unaligned_range[1]))
        if overlap==1 and unaligned_range[0] == unaligned_range[1]+1 and from_node.reference_name == \
                to_node.reference_name:
            unaligned_exon_bases += (to_node.reference_block_start + node_overlap) - from_node.reference_block_end - 1
        else:
            unaligned_exon_bases += max(0, overlap)
    return unaligned_exon_bases


def chain_alignments(head_nodes, node_dict, aln_graph, parent, children_coords, intervals,
                     parent_dict, lifted_features_list, distance_factor, allow_chrom_split, has_full_length):
    for node_name in head_nodes:
        add_edges(node_name, node_dict, aln_graph, parent, children_coords, intervals,
                  parent_dict, lifted_features_list, distance_factor, allow_chrom_split, has_full_length)


def add_edges(head_node_name, node_dict, aln_graph, parent, children_coords, intervals,
              parent_dict, lifted_features_list, distance_factor, allow_chrom_split, has_full_length):
    for node_name in node_dict:
        from_node = node_dict[node_name]
        if is_valid_edge(node_name, head_node_name, parent, intervals,
                         parent_dict, lifted_features_list, distance_factor, allow_chrom_split, has_full_length
                         ,node_dict, aln_graph):
            edge_weight = get_edge_weight(from_node, node_dict[head_node_name], children_coords, parent)
            aln_graph.add_edge(node_name, head_node_name, cost=edge_weight)


def is_valid_edge(from_node_name, to_node_name, parent, intervals, parent_dict, lifted_features_list, distance_factor,
                  allow_chrom_split, has_full_length,  node_dict, aln_graph):

    from_node, to_node = node_dict[from_node_name], node_dict[to_node_name]
    if has_full_length:
        allow_chrom_split = False
    if from_node.aln_id == to_node.aln_id:
        return False
    if   from_node.query_block_end >= to_node.query_block_end:
        return False
    if from_node.is_reverse != to_node.is_reverse:
        return False
    if from_node.reference_name != to_node.reference_name:
        if allow_chrom_split is False:
            return False
        if chrom_in_ancestors(aln_graph, from_node_name, to_node_name, node_dict):
            return False
    else:
        expected_distance = to_node.query_block_end - from_node.query_block_start
        actual_distance = to_node.reference_block_end - from_node.reference_block_start
        if to_node.reference_block_start < from_node.reference_block_end:
            return False
        if (actual_distance > distance_factor * expected_distance):
            return False
        if spans_overlap_region(from_node, to_node, parent, intervals, parent_dict,
                            lifted_features_list):
            return False
    return True



def chrom_in_ancestors(aln_graph, from_node_name, to_node_name, node_dict):
    for ancestor in nx.ancestors(aln_graph, from_node_name):
        if node_dict[ancestor].reference_name == node_dict[to_node_name].reference_name:
            return True
    return False





def add_target_node(aln_graph, node_dict, query_length, children_coords, parent):
    num_nodes = max(node_dict.keys())
    node_dict[num_nodes + 1] = aligned_seg.aligned_seg("end", "end", "end", query_length, query_length ,
                                                       query_length,
                                                       query_length, 0, [], False)
    aln_graph.add_node(num_nodes + 1)
    for node_name in node_dict:
        if is_terminal_node(node_name, aln_graph) and node_name != num_nodes+1:
            edge_weight = get_edge_weight(node_dict[node_name], node_dict[num_nodes + 1], children_coords, parent)
            aln_graph.add_edge(node_name, num_nodes + 1, cost=edge_weight)


def is_terminal_node(node, aln_graph):
    successors = aln_graph.successors(node)
    for successor in successors:
        return False
    return True


def find_shortest_path(node_dict, aln_graph):
    shortest_path = nx.shortest_path(aln_graph, source=0, target=len(node_dict) - 1,
                                     weight=lambda u, v, d: get_weight(u, v, d, aln_graph))
    shortest_path_nodes = []
    for i in range(1, len(shortest_path) - 1):
        node_name = shortest_path[i]
        shortest_path_nodes.append(node_dict[node_name])
    trim_path_boundaries(shortest_path_nodes)
    return shortest_path_nodes


def get_weight(u, v, d, G):
    node_u_wt = G.nodes[u].get('weight', 0)
    node_v_wt = G.nodes[v].get('weight', 0)
    edge_wt = d.get('cost', 1)
    return node_u_wt / 2 + node_v_wt / 2 + edge_wt


def trim_path_boundaries(shortest_path_nodes):

    for i in range(1, len(shortest_path_nodes)):
        from_node = shortest_path_nodes[i-1]
        to_node = shortest_path_nodes[i]
        node_overlap = get_node_overlap(from_node, to_node)
        if node_overlap !=0:
            print(node_overlap)
        to_node.query_block_start += node_overlap
        to_node.reference_block_start += node_overlap




def split_shortest_path_nodes(shortest_path_nodes):
    current_chrm = ""
    nodes_by_chrm = []
    current_chrm_group = []
    for node in shortest_path_nodes:
        if node.reference_name != current_chrm:
            if current_chrm_group != []:
                nodes_by_chrm.append(current_chrm_group)
            current_chrm_group = []
            current_chrm = node.reference_name
        current_chrm_group.append(node)
    if current_chrm_group != []:
        nodes_by_chrm.append(current_chrm_group)
    return nodes_by_chrm


def convert_all_children_coords(shortest_path_nodes, children, parent):
    shortest_path_nodes.sort(key=lambda x: x.query_block_start)
    mapped_children = {}
    total_bases, mismatches, insertions, deletions = 0, 0, 0, 0
    for child in children:
        total_bases += (child.end - child.start + 1)
        nearest_start_coord, nearest_end_coord = find_nearest_aligned_start_and_end(child.start, child.end,
                                                                                    shortest_path_nodes, parent)
        if nearest_start_coord != -1 and nearest_end_coord != -1:
            lifted_start, start_node = convert_coord(nearest_start_coord, shortest_path_nodes)
            lifted_end, end_node = convert_coord(nearest_end_coord, shortest_path_nodes)
            deletions += find_deletions(start_node, end_node, shortest_path_nodes)
            deletions += (nearest_start_coord - child.start) + (child.end - nearest_end_coord)
            mismatches += find_mismatched_bases(child.start, child.end, shortest_path_nodes, parent)
            insertions += find_insertions(start_node, end_node, shortest_path_nodes)
            strand = get_strand(shortest_path_nodes[0], parent)
            new_child = new_feature.new_feature(child.id, child.featuretype, shortest_path_nodes[0].reference_name,
                                                'Liftoff',
                                                strand, min(lifted_start, lifted_end) + 1,
                                                max(lifted_start, lifted_end) + 1, dict(child.attributes))
            mapped_children[new_child.id] = new_child
        else:
            deletions += (child.end - child.start + 1)
    return mapped_children, (total_bases - deletions) / total_bases, (total_bases - mismatches - (insertions +
                                                                                                  deletions)) / \
           total_bases


def find_nearest_aligned_start_and_end(child_start, child_end, shortest_path_nodes, parent):
    relative_coord1 = liftoff_utils.get_relative_child_coord(parent, child_start, shortest_path_nodes[0].is_reverse)
    relative_coord2 = liftoff_utils.get_relative_child_coord(parent, child_end, shortest_path_nodes[0].is_reverse)
    relative_start, relative_end = min(relative_coord1, relative_coord2), max(relative_coord1, relative_coord2)
    nearest_start = find_nearest_aligned_start(relative_start, relative_end, shortest_path_nodes)
    nearest_end = find_nearest_aligned_end(shortest_path_nodes, relative_end, relative_start)
    return nearest_start, nearest_end


def find_nearest_aligned_start(relative_start, relative_end, shortest_path_nodes):
    nearest_start = -1
    for node in shortest_path_nodes:
        if relative_start <= node.query_block_end:
            if relative_start >= node.query_block_start:
                nearest_start = relative_start
            else:
                if node.query_block_start < relative_end:
                    nearest_start = node.query_block_start
            return nearest_start
    return nearest_start


def find_nearest_aligned_end(shortest_path_nodes, relative_end, relative_start):
    nearest_end = -1
    for i in range(len(shortest_path_nodes)):
        node = shortest_path_nodes[i]
        if relative_end <= node.query_block_end:
            if relative_end >= node.query_block_start:
                nearest_end = relative_end
            else:
                if i > 0 and shortest_path_nodes[i - 1].query_block_end > relative_start:
                    nearest_end = shortest_path_nodes[i - 1].query_block_end
            break
    if nearest_end == -1 and node.query_block_end < relative_end and node.query_block_end > relative_start:
        nearest_end = node.query_block_end
    return nearest_end


def convert_coord(relative_coord, shortest_path_nodes):
    lifted_coord = 0
    node_index = 0
    for i in range(0, len(shortest_path_nodes)):
        if relative_coord >= shortest_path_nodes[i].query_block_start:
            if relative_coord <= shortest_path_nodes[i].query_block_end:
                lifted_coord = (shortest_path_nodes[i].reference_block_start + (
                        relative_coord - shortest_path_nodes[i].query_block_start))
                node_index = i
    return lifted_coord, node_index


def find_deletions(start, end, shortest_path_nodes):
    deletions = 0
    for i in range(start + 1, end + 1):
        to_node = shortest_path_nodes[i]
        from_node = shortest_path_nodes[i - 1]
        deletions += (to_node.query_block_start - from_node.query_block_end - 1)
    return deletions


def find_insertions(start, end, shortest_path_nodes):
    insertions = 0
    for i in range(start + 1, end + 1):
        to_node = shortest_path_nodes[i]
        from_node = shortest_path_nodes[i - 1]
        insertions += (to_node.reference_block_start  - from_node.reference_block_end - 1)
    return insertions


def find_mismatched_bases(start, end, shortest_path_nodes, parent):
    relative_coord1 = liftoff_utils.get_relative_child_coord(parent, start, shortest_path_nodes[0].is_reverse)
    relative_coord2 = liftoff_utils.get_relative_child_coord(parent, end, shortest_path_nodes[0].is_reverse)
    relative_start = min(relative_coord1, relative_coord2)
    relative_end = max(relative_coord1, relative_coord2)
    total_mismatches = 0
    for node in shortest_path_nodes:
        node_mismatches = np.array(node.mismatches)
        total_mismatches += len(
            node_mismatches[np.where((node_mismatches >= relative_start) & (node_mismatches <= relative_end))[0]])
    return total_mismatches
