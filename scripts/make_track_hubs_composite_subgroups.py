import os

if __name__ == '__main__':

    genome = 'mm39'
    bin_size = 100
    track_hub_folder = '/data/web/sites/PROseq'
    track_hub_url = 'http://upnaesrv1.epfl.ch/PROseq'

    #link tracks to track_hub folder
    link_to_data = '/data/web/sites/PROseq/tracks'
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
        
        # Bed tracks with gene phase
        fout.write(f"track gene_phase\n")
        fout.write("type bigBed 9\n")
        fout.write("itemRgb on\n")
        fout.write(f"shortLabel Gene phase\n")
        fout.write(f"longLabel Gene phase and amplitude mapped in RGB space (red: 0, yellow: pi/3, green: 2pi/3, cyan: pi, blue: 4pi/3, magenta:5pi/3)\n")
        fout.write(f"bigDataUrl {track_hub_url}/tracks/gene_phase.bb\n")
        fout.write("visibility pack\n")
        #fout.write("maxHeightPixels 200:40:10")
        fout.write("\n")
    
        # BigWig composite tracks
        fout.write(f"track PROseq\n")
        fout.write("compositeTrack on\n")
        fout.write("subGroup1 t Time 0=0 24=24 4=4 28=28 8=8 32=32 12=12 36=36 16=16 40=40 20=20 44=44\n")
        fout.write("subGroup2 s Strand forward=forward reverse=reverse\n")
        fout.write("dimensions dimX=s dimY=t\n")
        fout.write("sortOrder s=+ t=+\n")
        fout.write("shortLabel PROseq\n")
        fout.write("longLabel PRO-seq data composite track\n")
        fout.write("type bigWig\n")
        fout.write("visibility full\n")
        fout.write("maxHeightPixels 100:50:8\n")
        fout.write("autoScale on\n")
        fout.write("descriptionUrl PROseq.html\n")
        fout.write("\n")

        for strand in Strands:    
            for sample in Samples:
                name = '_'.join(sample.split('_')[:3])
                time = int(sample.split('_')[2][2:])
                if time>24:
                    time_label = str(time-24)+"+24 " + strand[0]
                    print(time,time_label)
                else:
                    time_label = str(time)

                fout.write(f"\ttrack {name}_{strand}\n")
                fout.write(f"\tparent PROseq on\n")
                fout.write(f"\tsubGroups t={time} s={strand}\n")
                fout.write(f"\tshortLabel {time_label}\n")
                fout.write(f"\tlongLabel \n")
                fout.write(f"\tbigDataUrl {track_hub_url}/tracks/{sample}/NormCoverage_3p_{strand}_bin{bin_size}bp.bw\n")
                fout.write("\ttype bigWig\n")
                if strand == 'forward':
                    fout.write("\tcolor 0,0,255\n")
                    fout.write("\tnegateValues off\n")
                    fout.write("\tminLimit 0\n")
                elif strand == 'reverse':
                    fout.write("\tcolor 255,0,0\n")
                    fout.write("\tnegateValues on\n")
                    fout.write("\tmaxLimit 0\n")
                fout.write(f"descriptionUrl {genome}/{name}_{strand}.html\n")
                fout.write(f"\n")

                # make html files
                outfile=f'{track_hub_folder}/{genome}/{name}_{strand}.html'
                with open(outfile,'w', encoding="utf-8") as fout2:
                    fout2.write(f"<h2>{name}_{strand}</h2>\n")



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
