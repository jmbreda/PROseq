import os

if __name__ == '__main__':
    track_hub_name = 'PROseq_phase'

    genome = 'mm39'
    track_hub_folder = f'/data/web/sites/{track_hub_name}'
    if not os.path.exists(track_hub_folder):
        os.makedirs(track_hub_folder)
    track_hub_url = f'http://upnaesrv1.epfl.ch/{track_hub_name}'
    results_folder = '/bigdata/jbreda/PROseq/results'

    #link tracks to track_hub folder
    link_to_data = f'{track_hub_folder}/tracks_bw_unbinned'
    data_folder = f'{results_folder}/norm_coverage/'
    if not os.path.exists(link_to_data):
        os.system(f'ln -s {data_folder} {link_to_data}')
    
    link_to_data = f'{track_hub_folder}/tracks_bw'
    data_folder = f'{results_folder}/binned_norm_coverage/'
    if not os.path.exists(link_to_data):
        os.system(f'ln -s {data_folder} {link_to_data}')

    link_to_data = f'{track_hub_folder}/tracks_bb'
    data_folder = f'{results_folder}/phase_amp/'
    if not os.path.exists(link_to_data):
        os.system(f'ln -s {data_folder} {link_to_data}')

    # make track hub
    outfile=f'{track_hub_folder}/hub.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"hub {track_hub_name}\n")
        fout.write(f"shortLabel {track_hub_name}\n")
        fout.write("longLabel PRO-seq data in mouse at times 0h-44h in steps of 4h, forward and reverse stand, with infered phase\n")
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
    outfile=f'{track_hub_folder}/{track_hub_name}.html'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("PRO-seq data in mouse at times 0h-44h in steps of 4h\n")

    # make url.txt
    outfile=f'{track_hub_folder}/url.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"{track_hub_url}/hub.txt\n")

    # make trackDb.txt
    T = range(0,48,4)
    Samples = [f'CT{t:02d}' for t in T]
    Strands = ['forward','reverse']

    os.makedirs(f"{track_hub_folder}/{genome}", exist_ok=True)
    outfile=f"{track_hub_folder}/{genome}/trackDb.txt"
    with open(outfile,'w', encoding="utf-8") as fout:

        # Bed tracks with bin phase
        for bin_size in [100,1000,10000]:
            for strand in Strands:
                fout.write(f"track bin_phase_{bin_size}_{strand}\n")
                fout.write("type bigBed 9\n")
                fout.write("itemRgb on\n")
                fout.write(f"shortLabel Bin phase {strand} {bin_size}bp\n")
                fout.write(f"longLabel {strand} bin {bin_size}bp phase and amplitude mapped in RGB space (bordeau: 0h, yellow: 6h, green: 12h, blue: 18h)\n")
                fout.write(f"bigDataUrl {track_hub_url}/tracks_bb/bin_phase_amp_{strand}_{bin_size}bp.bb\n")
                if bin_size == 1000:
                    fout.write("visibility dense\n")
                else:
                    fout.write("visibility hide\n")
                fout.write("\n")

        # Bed tracks with gene phase
        fout.write(f"track gene_phase\n")
        fout.write("type bigBed 9\n")
        fout.write("itemRgb on\n")
        fout.write(f"shortLabel Gene phase\n")
        fout.write(f"longLabel Gene phase and amplitude mapped in RGB space (bordeau: 0h, yellow: 6h, green: 12h, blue: 18h)\n")
        fout.write(f"bigDataUrl {track_hub_url}/tracks_bb/gene_phase_amp.bb\n")
        fout.write("visibility pack\n")
        fout.write("\n")
    
        # BigWig composite tracks with bin expression
        for bin_size in [1,100,1000,10000]:
            for strand in Strands:
                fout.write(f"track PROseq_{strand}_{bin_size}bp\n")
                fout.write("compositeTrack on\n")
                fout.write("subGroup1 t Time 0=0 24=24 4=4 28=28 8=8 32=32 12=12 36=36 16=16 40=40 20=20 44=44\n")
                fout.write("dimensions dimX=t\n")
                fout.write("sortOrder t=+\n")
                #fout.write("subGroup2 s Strand forward=forward reverse=reverse\n")
                #fout.write("dimensions dimX=s dimY=t\n")
                #fout.write("sortOrder s=+ t=+\n")
                fout.write(f"shortLabel PROseq {strand} {bin_size}bp\n")
                fout.write(f"longLabel PRO-seq data composite track {strand} {bin_size}bp\n")
                fout.write("type bigWig\n")
                if bin_size == 1000:
                    fout.write("visibility full\n")
                else:
                    fout.write("visibility hide\n")
                fout.write("maxHeightPixels 100:50:8\n")
                fout.write("autoScale group\n")
                fout.write("descriptionUrl PROseq.html\n")
                fout.write("\n")
  
                for sample in Samples:
                    name = sample
                    time = int(sample[2:])
                    if time>24:
                        time_label = str(time-24)+"+24"
                    else:
                        time_label = str(time)

                    fout.write(f"\ttrack {name}_{strand}_{bin_size}\n")
                    fout.write(f"\tparent PROseq_{strand}_{bin_size}bp on\n")
                    #fout.write(f"\tsubGroups t={time} s={strand}\n")
                    fout.write(f"\tsubGroups t={time}\n")
                    fout.write(f"\tshortLabel {time_label}\n")
                    fout.write(f"\tlongLabel \n")
                    if bin_size == 1:
                        fout.write(f"\tbigDataUrl {track_hub_url}/tracks_bw_unbinned/{sample}/NormCoverage_3p_{strand}.bw\n")
                    else:
                        fout.write(f"\tbigDataUrl {track_hub_url}/tracks_bw/{sample}/NormCoverage_3p_{strand}_bin{bin_size}bp.bw\n")
                    fout.write("\ttype bigWig\n")
                    if strand == 'forward':
                        fout.write("\tcolor 0,0,255\n")
                        fout.write("\tnegateValues off\n")
                        fout.write("\tminLimit 0\n")
                    elif strand == 'reverse':
                        fout.write("\tcolor 255,0,0\n")
                        fout.write("\tnegateValues on\n")
                        fout.write("\tmaxLimit 0\n")
                    fout.write(f"\tdescriptionUrl {genome}/{name}_{strand}.html\n")
                    fout.write(f"\n")

                    # make html files
                    outfile=f'{track_hub_folder}/{genome}/{name}_{strand}.html'
                    with open(outfile,'w', encoding="utf-8") as fout2:
                        fout2.write(f"<h2>{name}_{strand}</h2>\n")