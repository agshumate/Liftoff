print(len(snakmake.input))
for i in range (len(snakemake.input)):
    f_in=open(snakemake.input[i], 'r')
    f_out = open(snakemake.output[i], 'w')
    lines=f_in.readlines()

    found_values=[]

    for line in lines:
        if(line[0] != "#"):
            values=line.split("\t")
            feature=values[2] + ":" + values[3] + "-" + values[4]
            if feature not in found_values:
                f_out.write(line)
                found_values.append(feature)
