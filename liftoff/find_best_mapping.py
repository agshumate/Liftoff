import networkx as nx
from liftoff import aligned_seg, liftoff_utils
import numpy as np
import copy




def find_best_mapping(alignments, query_length,  parent, coords_to_exclude, children_dict, previous_gene_start, copy_tag):
    children = children_dict[parent.id]
    children_coords = liftoff_utils.merge_children_intervals(children)
    node_dict, aln_graph = intialize_graph()
    head_nodes = add_single_alignments(node_dict, aln_graph, alignments, children_coords, parent, coords_to_exclude,
                                       previous_gene_start)
    chain_alignments(head_nodes, node_dict, aln_graph, coords_to_exclude, parent, children_coords)
    add_target_node(aln_graph, node_dict, query_length, children_coords, parent)
    shortest_path = nx.shortest_path(aln_graph, source=0, target=len(node_dict) - 1,
                                     weight=lambda u, v, d: get_weight(u, v, d, aln_graph))
    shortest_path_weight = nx.shortest_path_length(aln_graph, source=0, target=len(node_dict) - 1,
                                                   weight=lambda u, v, d: get_weight(u, v, d, aln_graph))
    shortest_path_nodes = []
    for i in range  (1,len(shortest_path)-1):
        node_name = shortest_path[i]
        shortest_path_nodes.append(node_dict[node_name])
    if len(shortest_path_nodes) == 0:
        return {}, shortest_path_weight, 0,0
    mapped_children, alignment_coverage, seq_id = convert_all_children_coords(shortest_path_nodes, children, parent, copy_tag)
    return mapped_children, shortest_path_weight, alignment_coverage, seq_id



def add_target_node(aln_graph, node_dict, query_length, children_coords, parent):
    num_nodes = max(node_dict.keys())
    node_dict[num_nodes + 1] = aligned_seg.aligned_seg("end", "end", "end", query_length-1, query_length-1, query_length-1,
                                                       query_length-1, 0, [])
    for node_name in node_dict:
        if is_terminal_node(node_name, aln_graph):
            edge_weight = get_edge_weight(node_dict[node_name],node_dict[num_nodes+1], children_coords, parent)
            aln_graph.add_edge(node_name, num_nodes+1, cost=edge_weight)



def is_terminal_node(node, aln_graph):
    successors = aln_graph.successors(node)
    for successor in successors:
        return False
    return True





def chain_alignments(head_nodes, node_dict, aln_graph, coords_to_exclude, parent, children_coords):
    for node_name in head_nodes:
        add_edges(node_name, node_dict, aln_graph, coords_to_exclude, parent, children_coords)


def add_edges(head_node_name, node_dict, aln_graph, coords_to_exclude, parent, children_coords):
    for node_name in node_dict:
        from_node = node_dict[node_name]
        if is_valid_edge(node_dict[node_name], node_dict[head_node_name], coords_to_exclude, parent):
            edge_weight = get_edge_weight(from_node, node_dict[head_node_name], children_coords, parent)
            aln_graph.add_edge(node_name, head_node_name, cost=edge_weight)



def sort_alignments(previous_gene_start, alignments):
    order_dict = {}
    order = 0
    alignments=sorted(alignments, key=lambda x: np.abs(previous_gene_start - x.reference_block_start ))
    for aln in alignments:
        if aln.aln_id not in order_dict:
            order_dict[aln.aln_id] = order
            order += 1
    alignments=sorted(alignments, key=lambda x: (order_dict[x.aln_id], x.query_block_start))
    return alignments


def is_valid_alignment(previous_node, aln, coords_to_exclude, node_dict, parent):
    if aln.reference_block_start == -1:
        return False
    elif previous_node !=0:
        if spans_overlap_region(node_dict[previous_node], aln, coords_to_exclude, parent):
            return False
    return True


def add_single_alignments(node_dict, aln_graph, alignments, children_coords, parent, coords_to_exclude, previous_gene_start):
    alignments=sort_alignments(previous_gene_start, alignments)
    node_num = 1
    previous_node_id = -1
    head_nodes = []
    previous_node = 0
    for original_aln in alignments:
        aln = trim_overlap_coords(copy.deepcopy(original_aln), coords_to_exclude, parent)
        if aln.aln_id != previous_node_id:
            previous_node = 0
            previous_node_id = aln.aln_id
        is_valid_aln = is_valid_alignment(previous_node, aln, coords_to_exclude, node_dict, parent)
        if is_valid_aln:
            if previous_node == 0:
                head_nodes.append(node_num)
            node_dict[node_num]=aln
            current_node = node_num
            node_weight = get_node_weight(aln, children_coords, parent)
            aln_graph.add_node(node_num, weight=node_weight)
            edge_weight = get_edge_weight(node_dict[previous_node], node_dict[current_node], children_coords, parent)
            aln_graph.add_edge(previous_node, current_node, cost=edge_weight)
            previous_node = current_node
            node_num +=1
    return head_nodes





def convert_all_children_coords(shortest_path_nodes, children, parent, copy_tag):
    shortest_path_nodes.sort(key=lambda x: x.query_block_start)
    mapped_children = {}
    aligned_bases, total_bases, mismatches= set([]), set([]), set([])
    for child in children:
        child_start, child_end = child.start, child.end
        lifted_feature = False
        total_bases.update(range(child_start, child_end + 1))
        while lifted_feature is False:
            lifted_start, lifted_end = convert_coord(child_start, child_end, parent, shortest_path_nodes)
            if lifted_start !=0 and lifted_end !=0:
                lifted_feature = True
            else:
                if lifted_start == 0:
                    child_start += 1
                if lifted_end == 0:
                    child_end -= 1
            if (child_start >= child_end and lifted_start ==0 and lifted_end ==0) :
                lifted_feature = True
        aligned_bases.update(range(child_start,child_end+1))
        mismatched_bases = find_mismatched_bases(child_start, child_end, shortest_path_nodes, parent)
        mismatches.update(mismatched_bases)
        if  lifted_start !=0:
            strand = get_strand(shortest_path_nodes[0], parent)
            new_child= liftoff_utils.make_new_feature(copy.copy(child), min(lifted_start, lifted_end) + 1, max(lifted_start, lifted_end) + 1, strand, shortest_path_nodes[0].reference_name)
            mapped_children[new_child.id] = new_child
    return mapped_children, len(aligned_bases)/len(total_bases), (len(aligned_bases)-len(mismatches))/len(total_bases)


def find_mismatched_bases(start,end, shortest_path_nodes, parent):
    all_mismatches = []
    relative_start=liftoff_utils.get_relative_child_coord(parent, start, shortest_path_nodes[0].is_reverse)
    relative_end=liftoff_utils.get_relative_child_coord(parent, end, shortest_path_nodes[0].is_reverse)
    for node in shortest_path_nodes:
        node_mismatches = np.array(node.mismatches)
        mismatches = node_mismatches[np.where((node_mismatches >= relative_start) & (node_mismatches <= relative_end ))[0]]
        if len(mismatches) > 0:
            all_mismatches.extend(mismatches.tolist())
    return all_mismatches



def get_strand(aln, parent):
    if aln.is_reverse:
        if parent.strand == "-":
            strand = '+'
        else:
            strand = "-"
    else:
        strand =  parent.strand
    return strand


def convert_coord(coord_start, coord_end, parent, shortest_path_nodes):
    lifted_coords = []
    for i in [coord_start, coord_end]:
        lifted_coord = 0
        relative_coord = liftoff_utils.get_relative_child_coord(parent, i, shortest_path_nodes[0].is_reverse)
        for j in range (0,len(shortest_path_nodes)):
            if relative_coord >= shortest_path_nodes[j].query_block_start:
                if relative_coord <= shortest_path_nodes[j].query_block_end:
                    lifted_coord=(shortest_path_nodes[j].reference_block_start + (
                                relative_coord - shortest_path_nodes[j].query_block_start))
        lifted_coords.append(lifted_coord)
    return lifted_coords[0], lifted_coords[1]


def intialize_graph():
    aln_graph = nx.DiGraph()
    node_dict = {}
    node_dict[0] = aligned_seg.aligned_seg("start", "start", "start", 0, 0, 0, 0, 0, [])
    aln_graph.add_node(0)
    return node_dict, aln_graph


def get_node_weight(aln, children_coords, parent):
    weight = 0
    for child_interval in children_coords:
        relative_start = liftoff_utils.get_relative_child_coord(parent, child_interval[0], aln.is_reverse)
        relative_end = liftoff_utils.get_relative_child_coord(parent, child_interval[1], aln.is_reverse)
        child_start, child_end = min(relative_start, relative_end), max(relative_start, relative_end)
        weight += len(aln.mismatches[(aln.mismatches >=child_start) & (aln.mismatches <= child_end)])
    return weight



def get_weight(u, v, d, G):
    node_u_wt = G.nodes[u].get('weight', 0)
    node_v_wt = G.nodes[v].get('weight', 0)
    edge_wt = d.get('cost', 1)
    return node_u_wt/2 + node_v_wt/2 + edge_wt



def get_edge_weight(from_node, to_node, children_coords, parent):
    unaligned_range = [from_node.query_block_end + 1, to_node.query_block_start-1]
    unaligned_exon_bases = 0
    for child_interval in children_coords:
        if from_node.reference_name == "start":
            is_reverse = to_node.is_reverse
        else:
            is_reverse=from_node.is_reverse
        relative_start= liftoff_utils.get_relative_child_coord(parent, child_interval[0], is_reverse)
        relative_end = liftoff_utils.get_relative_child_coord(parent, child_interval[1], is_reverse)
        child_start, child_end = min(relative_start, relative_end), max(relative_start, relative_end)
        overlap = liftoff_utils.count_overlap(child_start, child_end, unaligned_range[0], unaligned_range[1])
        if overlap == 1 and unaligned_range[0] == unaligned_range[1]:
            unaligned_exon_bases += to_node.reference_block_start - from_node.reference_block_end +1
        else:
            unaligned_exon_bases += max(0,overlap)
    return unaligned_exon_bases



def is_valid_edge(from_node, to_node, coords_to_exclude, parent):
    if from_node.aln_id == to_node.aln_id:
        return False
    if from_node.query_block_end > to_node.query_block_start:
        return False
    if from_node.is_reverse != to_node.is_reverse:
        return False
    if from_node.reference_name != to_node.reference_name:
        return False
    expected_distance = to_node.query_block_end - from_node.query_block_start
    actual_distance = to_node.reference_block_end - from_node.reference_block_start
    if to_node.reference_block_start < from_node.reference_block_end:
        return False
    if (actual_distance > 2 * expected_distance):
        return False
    if spans_overlap_region(from_node, to_node, coords_to_exclude, parent):
        return False
    return True


def spans_overlap_region(from_node, to_node, coords_to_exclude, parent):
    target_chrm = from_node.reference_name
    strand = get_strand(from_node,parent)
    for coords in coords_to_exclude:
        if strand == coords[3] and target_chrm == coords[2]:
            if liftoff_utils.count_overlap(coords[0], coords[1], from_node.reference_block_end, to_node.reference_block_start) > 0:
                return True
    return False


def trim_overlap_coords(aln, coords_to_exclude, parent):
    ref_start, ref_end = aln.reference_block_start, aln.reference_block_end
    target_chrm = aln.reference_name
    strand = get_strand(aln,parent)
    for coords in coords_to_exclude:
        if strand == coords[3] and target_chrm == coords[2]:
            if coords[0] <= ref_start and coords[1] >= ref_end:
                ref_start, ref_end= -1,-1
                break
            elif coords[1] < ref_start or coords[0] > ref_end:
                ref_start, ref_end = ref_start, ref_end
            elif coords[0] >= ref_start and coords[1] <= ref_end:
                if coords[0]- ref_start >  ref_end - coords[1]:
                    ref_start, ref_end= ref_start, coords[0] -1
                else:
                    ref_start, ref_end = coords[1] +1, ref_end
            elif coords[0] > ref_start:
                ref_start, ref_end =  ref_start, coords[0]-1
            else:
                ref_start, ref_end = coords[1]+1, ref_end
    start_diff = ref_start - aln.reference_block_start
    end_diff = aln.reference_block_end - ref_end
    aln.reference_block_start = ref_start
    aln.reference_block_end = ref_end
    aln.query_block_start += start_diff
    aln.query_block_end -= end_diff
    return aln