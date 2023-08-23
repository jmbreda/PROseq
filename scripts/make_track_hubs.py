import os

if __name__ == '__main__':

    genome = 'mm39'
    outfolder = 'track_hubs'
    os.makedirs(f"{outfolder}/{genome}", exist_ok=True)

    # make track hub
    outfile=f'{outfolder}/hub.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("hub PRO-seq\n")
        fout.write("shortLabel PRO-seq\n")
        fout.write("longLabel PRO-seq data in mouse at times 0h-44h in steps of 4h\n")
        fout.write("genomesFile genomes.txt\n")
        fout.write("email jeremie.breda@epfl.ch\n")
        fout.write("descriptionUrl PROseq.html\n")
        fout.write("\n")

    # make genomes.txt
    outfile=f'{outfolder}/genomes.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"genome {genome}\n")
        fout.write(f"trackDb {genome}/trackDb.txt\n")
        fout.write("\n")

    # make ChIP_Atlas.html
    outfile=f'{outfolder}/PROseq.html'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("PRO-seq data in mouse at times 0h-44h in steps of 4h\n")

    # make url.txt
    outfile=f'{outfolder}/url.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("http://upnaesrv1.epfl.ch/PROseq/hub.txt\n")

    # make trackDb.txt
    Samples = os.listdir('results/star/')
    Samples.sort()
    Strands = ['forward','reverse']
    outfile=f'{outfolder}/{genome}/trackDb.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        for strand in Strands:
            for sample in Samples:
                print(sample)
                name = '_'.join(sample.split('_')[:3])

                # BigWig tracks
                fout.write(f"track {name} {strand}\n")
                fout.write("type bigWig\n")
                fout.write(f"bigDataUrl http://upnaesrv1.epfl.ch/PROseq/tracks/{sample}/NormCoverage_{strand}.bw\n")
                #fout.write("type bam\n")
                #fout.write(f"bigDataUrl http://upnaesrv1.epfl.ch/PROseq/tracks/{sample}/Aligned.sortedByCoord.out.bam\n")
                fout.write(f"shortLabel {name} {strand}\n")
                fout.write(f"longLabel {sample} {strand}\n")
                fout.write("visibility full\n")
                fout.write("autoScale on\n")
                fout.write(f"\n")

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
