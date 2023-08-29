import os

if __name__ == '__main__':

    genome = 'mm39'
    bin_size = 100
    track_hub_folder = '/data/web/sites/PROseq_0_4_44h'
    track_hub_url = 'http://upnaesrv1.epfl.ch/PROseq_0_4_44h'

    #link tracks to track_hub folder
    link_to_data = '/data/web/sites/PROseq_0_4_44h/tracks'
    data_folder = '/bigdata/jbreda/PROseq/results/binned_norm_counts/'
    if not os.path.exists(link_to_data):
        os.system(f'ln -s {data_folder} {link_to_data}')

    # make track hub
    outfile=f'{track_hub_folder}/hub.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("hub PRO-seq\n")
        fout.write("shortLabel PRO-seq\n")
        fout.write("longLabel PRO-seq data in mouse at times 0h-44h in steps of 4h, forward and reverse stand\n")
        fout.write("genomesFile genomes.txt\n")
        fout.write("email jeremie.breda@epfl.ch\n")
        fout.write("descriptionUrl PROseq.html\n")
        fout.write("\n")

    # make genomes.txt
    outfile=f'{track_hub_folder}/genomes.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"genome {genome}\n")
        fout.write(f"trackDb {genome}/trackDb.txt\n")
        fout.write("\n")

    # make ChIP_Atlas.html
    outfile=f'{track_hub_folder}/PROseq.html'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("PRO-seq data in mouse at times 0h-44h in steps of 4h\n")

    # make url.txt
    outfile=f'{track_hub_folder}/url.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"{track_hub_url}/hub.txt\n")

    # make trackDb.txt
    Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)]
    Strands = ['forward','reverse']

    os.makedirs(f"{track_hub_folder}/{genome}", exist_ok=True)
    outfile=f"{track_hub_folder}/{genome}/trackDb.txt"
    with open(outfile,'w', encoding="utf-8") as fout:
        for sample in Samples:
            name = '_'.join(sample.split('_')[:3])
            time = int(sample.split('_')[2][2:])

            # BigWig composite tracks
            fout.write(f"track {name}\n")
            fout.write("type bigWig\n")
            fout.write("container multiWig\n")
            fout.write(f"shortLabel t={time}h\n")
            fout.write(f"longLabel {sample}\n")
            fout.write("visibility full\n")
            fout.write("aggregate solidOverlay\n")
            fout.write("showSubtrackColorOnUi on\n")
            #fout.write("maxHeightPixels 500:100:8\n")
            fout.write("autoScale on\n")
            fout.write(f"html {name}.html\n")
            fout.write("\n")
            
            for strand in Strands:

                fout.write(f"\ttrack {name}_{strand}\n")
                fout.write(f"\tshortLabel {time}_{strand[0]}\n")
                fout.write(f"\tlongLabel {sample}_{strand}\n")
                fout.write(f"\tbigDataUrl {track_hub_url}/tracks/{sample}/NormCoverage_3p_{strand}_bin{bin_size}bp.bw\n")
                fout.write(f"\tparent {name}\n")
                fout.write("\ttype bigWig\n")
                if strand == 'forward':
                    fout.write("\tcolor 0,0,255\n")
                elif strand == 'reverse':
                    fout.write("\tcolor 255,0,0\n")
                    fout.write("\tnegateValues on\n")
                fout.write(f"\n")

        
            # make description file
            outfile=f'{track_hub_folder}/{genome}/{name}.html'
            with open(outfile,'w', encoding="utf-8") as fout2:
                fout2.write(f"{sample}_{strand}\nTime point: {name.split('_')[2][2:]}h\nTop: forward stand\nBottom: reverse strand\n")


                # BAM tracks
                # Optional settings
                #
                # pairEndsByName  any value                   # presence indicates paired-end alignments
                # pairSearchRange N                           # max distance between paired alignments, default 20,000 bases
                # bamColorMode    strand|gray|tag|off         # coloring method, default is strand
                # bamGrayMode     aliQual|baseQual|unpaired   # grayscale metric, default is aliQual
                # bamColorTag     XX                          # optional tag for RGB color, default is "YC"
                # minAliQual      N                           # display only items with alignment quality at least N, default 0
                # showNames       on|off                      # if off, don't display query names, default is on
                #
                # name            track label                 # default is "User Track"
                # description     center label                # default is "User Supplied Track"
                # visibility      squish|pack|full|dense|hide # default is hide (will also take numeric values 4|3|2|1|0)
                # bigDataUrl      https://your.bam.bai.com    # default is the bigDataUrl with .bai added to the end
                # priority        N                           # default is 100
                # db              genome database             # e.g. hg18 for Human Mar. 2006
                # maxWindowToDraw N                           # don't display track when viewing more than N bases
                # chromosomes     chr1,chr2,...               # track contains data only on listed reference assembly sequences
                # doWiggle        on|off                      # if on, show data as density graph, default is off
    #
                #fout.write(f"track {name}\n")
                #fout.write("type bam\n")
                #fout.write(f"bigDataUrl http://upnaesrv1.epfl.ch/PROseq/tracks/{sample}/Aligned.sortedByCoord.out.bam\n")
                ##fout.write(f"name={name} ")
                #fout.write(f"description={sample} ")
                #fout.write("visibility=full ")
                #fout.write(f"db={genome} ")
                #fout.write("chromosomes=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,ch11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX,chrY ")
                #fout.write("autoScale on\n")
                #fout.write("\n")
