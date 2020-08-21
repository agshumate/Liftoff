input_gff = [line.rstrip() for line in open("/ccb/salz7-data/ftp.ccb/pub/data/Triticum_aestivum/Triticum_aestivum_4.0"
                           "/Triticum_aestivum_4"
                      ".0_v1"
                 ".0_ncbi_ids.gff", 'r').readlines()]
id_map = [line.rstrip() for line in open("/ccb/salz1/dpuiu/Triticum_aestivum/Assembly/Triticum_4.0/Annotation"
                                      "/NMPL03_accs",
                      'r').readlines()]

id_maps = {}
for full_line in id_map:
    line = full_line.split()
    id_maps[line[1]] = line[0]

for full_line in input_gff:
    line = full_line.split("\t")
    if line[0][0] != "#":
        new_chrom = id_maps[line[0]]

        line[0] = new_chrom
        new_line = ""
        for column in line:
            new_line += column + "\t"
        print(new_line[:-1])
    else:
        print(full_line)





