import os

os.getcwd()
os.chdir("/Users/pramod/Desktop/Proteomics_Hippocampus/")

def parse_census(filename):
    current_protein = None
    protein_data = list()
    sequence_data = list()

    with open(filename, "r") as infile:
        for line in infile.readlines():
            line = line.strip().split('\t')
            if line[0] == "H":
                if line[1] == "PLINE":
                    header_protein = '\t'.join(line)
                elif line[1] == "SLINE":
                    header_sequence = '\t'.join(["PROTEIN"] + line[1:])
                else:
                    continue
            if line[0] == "P":
                current_protein = line[1]
                protein_data.append('\t'.join(line))
            elif line[0] == "S":
                line = [current_protein] + line
                sequence_data.append('\t'.join(line))
        protein_data = [header_protein] + protein_data
        sequence_data = [header_sequence] + sequence_data
    
    return(tuple([protein_data, sequence_data]))

if __name__ == "__main__":

    filenames = [
        "census-out-Hipp-E17.5-21488.txt",
        "census-out-Hipp-p7-21439.txt",
        "census-out-Hipp-p35-21348.txt"
    ]

    output_filenames = [
        "data/census-out-Hipp-E17.5-21488.txt",
        "data/census-out-Hipp-p7-21439.txt",
        "data/census-out-Hipp-p35-21348.txt"
    ]

    for f, o in zip(filenames, output_filenames):
        with open(o, "w") as ofile:
            ofile.write('\n'.join(parse_census(f)[1]))
