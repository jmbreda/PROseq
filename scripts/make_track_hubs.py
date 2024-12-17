import os

if __name__ == '__main__':

    # Parameters
    Genome = ['mm10','mm39']
    T = range(0,48,4)
    Samples = [f'CT{t:02d}' for t in T]
    Strands = ['forward','reverse']

    # Track hub name and url
    track_hub_name = 'PROseq'
    track_hub_url = f"https://sv-open.epfl.ch/upnae-public/sites/{track_hub_name}"

    # Create track hub folder
    track_hub_folder = f"/data/shared/sv-open/sites/{track_hub_name}"
    if not os.path.exists(track_hub_folder):
        os.makedirs(track_hub_folder)
    
    # make track hub file
    outfile=f'{track_hub_folder}/hub.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"hub {track_hub_name}\n")
        fout.write(f"shortLabel {track_hub_name}\n")
        fout.write("longLabel PRO-seq data in mouse at times 0h-44h in steps of 4h, forward and reverse stand, with infered phase\n")
        fout.write("genomesFile genomes.txt\n")
        fout.write("email jeremie.breda@epfl.ch\n")
        fout.write(f"descriptionUrl {track_hub_name}.html\n")
        fout.write("\n")

    # make genomes.txt
    outfile=f'{track_hub_folder}/genomes.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        for genome in Genome:
            fout.write(f"genome {genome}\n")
            fout.write(f"trackDb {genome}/trackDb.txt\n")
            fout.write("\n")

    # make html description
    outfile=f'{track_hub_folder}/{track_hub_name}.html'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("PRO-seq data in mouse at times 0h-44h in steps of 4h\n")

    # make url.txt
    outfile=f'{track_hub_folder}/url.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"{track_hub_url}/hub.txt\n")
    
    # make trackDb.txt
    for genome in Genome:

        if genome == 'mm10':
            track_folder = f"{track_hub_url}/results/GRCm38"
        elif genome == 'mm39':
            track_folder = f"{track_hub_url}/results/GRCm39"

        os.makedirs(f"{track_hub_folder}/{genome}", exist_ok=True)
        outfile=f"{track_hub_folder}/{genome}/trackDb.txt"
        with open(outfile,'w', encoding="utf-8") as fout:

            p = 0
            # Bed tracks with gene phase
            fout.write(f"track gene_phase\n")
            fout.write("type bigBed 9\n")
            fout.write("itemRgb on\n")
            fout.write(f"shortLabel Gene phase\n")
            fout.write(f"longLabel Gene phase and amplitude mapped in RGB space (red: 0h, yellow: 6h, green: 12h, blue: 18h)\n")
            fout.write(f"bigDataUrl {track_folder}/phase_amp/gene_phase_amp.bb\n")
            fout.write("visibility pack\n")
            fout.write(f"priority {p}\n")
            fout.write("\n")
            
            p+=1

            # Bed tracks with kalman smoothing phase
            for strand in Strands:
                fout.write(f"track gene_kalman_smoothing_phase_{strand}\n")
                fout.write("type bigBed 9\n")
                fout.write("itemRgb on\n")
                fout.write(f"shortLabel Gene kalman smoothing phase {strand}\n")
                fout.write(f"longLabel Gene phase by kalman filter {strand} strand mapped in RGB space (red: 0h, yellow: 6h, green: 12h, blue: 18h)\n")
                fout.write(f"bigDataUrl {track_folder}/phase_amp/gene_kalman_phase_R2_{strand}_1000bp.bb\n")
                fout.write("visibility dense\n")
                fout.write(f"priority {p}\n")
                fout.write("\n")

                p+=1

            # Bed tracks with bin phase
            for bin_size in [1000]:#[100,1000,10000]:
                for strand in Strands:
                    fout.write(f"track bin_phase_{bin_size}_{strand}\n")
                    fout.write("type bigBed 9\n")
                    fout.write("itemRgb on\n")
                    fout.write(f"shortLabel Bin phase {strand} {bin_size}bp\n")
                    fout.write(f"longLabel {strand} bin {bin_size}bp phase and amplitude mapped in RGB space (red: 0h, yellow: 6h, green: 12h, blue: 18h)\n")
                    fout.write(f"bigDataUrl {track_folder}/phase_amp/bin_phase_amp_{strand}_{bin_size}bp.bb\n")
                    if bin_size == 1000:
                        fout.write("visibility dense\n")
                    else:
                        fout.write("visibility hide\n")
                    fout.write(f"priority {p}\n")
                    fout.write("\n")
                    p += 1

            # Bed tracks with expressed regions
            for bin_size in [1000]:#[100,1000,10000]:
                fout.write(f"track expressed_regions_{bin_size}\n")
                fout.write("type bigBed 3\n")
                #fout.write("itemRgb on\n")
                fout.write(f"shortLabel Expressed regions {bin_size}bp\n")
                fout.write(f"longLabel bin {bin_size}bp Expressed regions\n")
                fout.write(f"bigDataUrl {track_folder}/binned_norm_coverage/expressed_regions_bin{bin_size}bp.bb\n")
                if bin_size == 1000:
                    fout.write("visibility dense\n")
                else:
                    fout.write("visibility hide\n")
                fout.write(f"priority {p}\n")
                fout.write("\n")
                p += 1

            # bed tracks with extended kalman smoothing phase and amp on expressed regions
            for bin_size in [1000]:#[100,1000,10000]:
                for strand in Strands:
                    fout.write(f"track kalman_expressed_regions_phase_amp_{strand}_{bin_size}\n")
                    fout.write("type bigBed 9\n")
                    fout.write("itemRgb on\n")
                    fout.write(f"shortLabel Kalman expressed regions phi amp {strand} {bin_size}bp\n")
                    fout.write(f"longLabel Extended Kalman smoothing phase and amplitude on expressed regions {strand} strand {bin_size}bp\n")
                    fout.write(f"bigDataUrl {track_folder}/kalman/extended_kalman_on_expressed_regions_{strand}_bin{bin_size}bp_phi_amp.bb\n")
                    if bin_size == 1000:
                        fout.write("visibility dense\n")
                    else:
                        fout.write("visibility hide\n")
                    fout.write(f"priority {p}\n")
                    fout.write("\n")
                    p += 1

            # bed tracks with extended kalman smoothing loglikelihood (transformed in 0-1000 range) on expressed regions
            for bin_size in [1000]:#[100,1000,10000]:
                for strand in Strands:
                    fout.write(f"track kalman_expressed_regions_ll_{strand}_{bin_size}\n")
                    fout.write("type bigWig\n")
                    fout.write("itemRgb on\n")
                    fout.write(f"shortLabel Kalman expressed regions ll {strand} {bin_size}bp\n")
                    fout.write(f"longLabel Extended Kalman smoothing on expressed regions {strand} strand {bin_size}bp\n")
                    fout.write(f"bigDataUrl {track_folder}/kalman/extended_kalman_on_expressed_regions_{strand}_bin{bin_size}bp_ll.bw\n")
                    if bin_size == 1000:
                        fout.write("visibility dense\n")
                    else:
                        fout.write("visibility hide\n")
                    if strand == 'forward':
                        fout.write("\tcolor 0,0,255\n")
                    elif strand == 'reverse':
                        fout.write("\tcolor 255,0,0\n")
                    fout.write("autoScale on\n")
                    fout.write("\tminLimit 0\n")
                    fout.write("maxHeightPixels 100:30:8\n") # max:default:min
                    fout.write(f"priority {p}\n")
                    fout.write("\n")
                    p += 1

        
            # BigWig composite tracks with bin expression
            for bin_size in [1000]:#[1,100,1000,10000]:
                for strand in Strands:
                    fout.write(f"track PROseq_{strand}_{bin_size}bp\n")
                    fout.write("compositeTrack on\n")
                    fout.write("subGroup1 t Time CT00=00 CT04=04 CT08=08 CT12=12 CT16=16 CT20=20 CT28=28 CT24=24 CT32=32 CT36=36 CT40=40 CT44=44\n")
                    fout.write("dimensions dimX=t\n")
                    fout.write("sortOrder t=+\n")
                    #fout.write("subGroup2 s Strand forward=forward reverse=reverse\n")
                    #fout.write("dimensions dimX=s dimY=t\n")
                    #fout.write("sortOrder s=+ t=+\n")
                    fout.write(f"shortLabel PROseq {strand} {bin_size}bp\n")
                    fout.write(f"longLabel PRO-seq data composite track {strand} {bin_size}bp (sum normed count per bin + 1, 1bp: norm count)\n")
                    fout.write("type bigWig\n")
                    if bin_size == 1000:
                        fout.write("visibility full\n")
                    else:
                        fout.write("visibility hide\n")
                    fout.write("maxHeightPixels 100:30:8\n") # max:default:min
                    fout.write("autoScale group\n")
                    fout.write("descriptionUrl PROseq.html\n")
                    fout.write(f"priority {p}\n")
                    fout.write("\n")
                    p += 1
    
                    for sample in Samples:
                        name = sample
                        time = sample[2:]

                        fout.write(f"\ttrack {name}_{strand}_{bin_size}\n")
                        fout.write(f"\tparent PROseq_{strand}_{bin_size}bp on\n")
                        #fout.write(f"\tsubGroups t={time} s={strand}\n")
                        fout.write(f"\tsubGroups t={sample}\n")
                        fout.write(f"\tshortLabel {sample}\n")
                        fout.write(f"\tlongLabel \n")
                        if bin_size == 1:
                            fout.write(f"\tbigDataUrl {track_folder}/norm_coverage/{sample}/NormCoverage_3p_{strand}.bw\n")
                        else:
                            #fout.write(f"\tbigDataUrl {track_hub_url}/tracks_bw/{sample}/Log2NormCoverage_3p_{strand}_bin{bin_size}bp.bw\n")
                            fout.write(f"\tbigDataUrl {track_folder}/binned_norm_coverage/{sample}/NormCoverage_3p_{strand}_bin{bin_size}bp.bw\n")
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

                        # make html description file
                        outfile=f'{track_hub_folder}/{genome}/{name}_{strand}.html'
                        with open(outfile,'w', encoding="utf-8") as fout2:
                            fout2.write(f"<h2>{name}_{strand}</h2>\n")




        #link tracks to track_hub folder
        #link_to_data = f'{track_hub_folder}/tracks_bw_unbinned'
        #data_folder = f'{results_folder}/norm_coverage/'
        #if not os.path.exists(link_to_data):
        #    os.system(f'ln -s {data_folder} {link_to_data}')
        
        #link_to_data = f'{track_hub_folder}/tracks_tad'
        #data_folder = f'/bigdata/jbreda/PROseq/resources/TAD/'
        #if not os.path.exists(link_to_data):
        #    os.system(f'ln -s {data_folder} {link_to_data}')

        
