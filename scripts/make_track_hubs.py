import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Make track hub')
    parser.add_argument('--track_hub_name', help='Track hub name', type=str)
    parser.add_argument('--track_hub_folder', help='Track hub folder', type=str)
    parser.add_argument('--hub', help='Track hub file', type=str)
    parser.add_argument('--genome', help='Genome version', type=str)
    parser.add_argument('--url', help='Track hub url file', type=str)
    parser.add_argument('--trackDb', help='TrackDb file', type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    args = parse_args()

    # Parameters
    Genome = ['mm10']
    T = range(0,48,4)
    Samples = [f'CT{t:02d}' for t in T]
    Strands = ['forward','reverse']
    Bin_size = {100:'100bp', 1000:'1kb', 10000:'10kb'}
    default_bin_size = 1000

    # make track hub file
    with open(args.hub,'w', encoding="utf-8") as fout:
        fout.write(f"hub {args.track_hub_name}\n")
        fout.write(f"shortLabel {args.track_hub_name}\n")
        fout.write("longLabel PRO-seq data in mouse at times 0h-44h in steps of 4h, forward and reverse stand, with infered phase\n")
        fout.write("genomesFile genomes.txt\n")
        fout.write("email jeremie.breda@epfl.ch\n")
        description_html = f"{args.track_hub_name}.html"
        fout.write(f"descriptionUrl {description_html}\n")
        # make description html file
        with open(f'{args.track_hub_folder}/{description_html}','w', encoding="utf-8") as fout2:
            fout2.write("PRO-seq data in mouse liver at times 0h-44h in steps of 4h\n")
        fout.write("\n")

    # make genomes.txt
    with open(args.genome,'w', encoding="utf-8") as fout:
        for genome in Genome:
            fout.write(f"genome {genome}\n")
            fout.write(f"trackDb {genome}/trackDb.txt\n")
            fout.write("\n")

    # define url and save in url.txt
    track_hub_url = f"https://sv-open.epfl.ch/upnae-public/sites/{args.track_hub_name}"
    with open(args.url,'w', encoding="utf-8") as fout:
        fout.write(f"{track_hub_url}/hub.txt\n")
    
    # make trackDb.txt
    for genome in Genome:

        if genome == 'mm10':
            track_folder = f"{track_hub_url}/results/GRCm38"
        elif genome == 'mm39':
            track_folder = f"{track_hub_url}/results/GRCm39"

        outfile=f"{args.track_hub_folder}/{genome}/trackDb.txt"
        with open(outfile,'w', encoding="utf-8") as fout:

            p = 1
            # Gene phase and amplitude bigBed
            for strand in Strands:
                fout.write(f"track gene_phase_{strand[0]}\n")
                fout.write("type bigBed 9\n")
                fout.write("itemRgb on\n")
                fout.write(f"shortLabel Gene phase {strand[0]}\n")
                fout.write(f"longLabel Gene phase and amplitude {strand[0]} strand, mapped in RGB space (red: 0h, yellow: 6h, green: 12h, blue: 18h)\n")
                fout.write(f"bigDataUrl {track_folder}/phase_amp/gene_phase_amp_{strand}.bb\n")
                fout.write("visibility pack\n")
                #if strand == 'forward':
                #    fout.write("\tcolor 0,0,255\n")
                #elif strand == 'reverse':
                #    fout.write("\tcolor 255,0,0\n")
                fout.write(f"priority {p}\n")
                fout.write("\n")
            
                p+=1

            # kalman smoothing gene phase
            if False:
                for bin_size in [1000,10000]:
                    for strand in Strands:
                        fout.write(f"track gene_kalman_smoothing_phase_{Bin_size[bin_size]}_{strand[0]}\n")
                        fout.write("type bigBed 9\n")
                        fout.write("itemRgb on\n")
                        fout.write(f"shortLabel Gene kalman phase {strand[0]} {Bin_size[bin_size]}\n")
                        fout.write(f"longLabel Gene phase by kalman filter {strand} strand mapped in RGB space (red: 0h, yellow: 6h, green: 12h, blue: 18h)\n")
                        fout.write(f"bigDataUrl {track_folder}/phase_amp/gene_kalman_phase_R2_{strand}_{bin_size}bp.bb\n")
                        if bin_size == default_bin_size:
                            fout.write("visibility dense\n")
                        else:
                            fout.write("visibility hide\n")
                        if strand == 'forward':
                            fout.write("\tcolor 0,0,255\n")
                        elif strand == 'reverse':
                            fout.write("\tcolor 255,0,0\n")
                        fout.write(f"priority {p}\n")
                        fout.write("\n")

                        p+=1

            # Bin mean log2 expression bigwig
            for bin_size in [100,1000,10000]:
                for strand in Strands:
                    fout.write(f"track bin_mu_{Bin_size[bin_size]}_{strand[0]}\n")
                    fout.write("type bigWig\n")
                    fout.write("itemRgb on\n")
                    fout.write(f"shortLabel Bin mu {strand[0]} {Bin_size[bin_size]}\n")
                    fout.write(f"longLabel Bin mean log2 {strand} {Bin_size[bin_size]}\n")
                    fout.write(f"bigDataUrl {track_folder}/phase_amp/bin_mu_{strand}_{bin_size}bp.bw\n")
                    #fout.write("viewLimits 0:24\n")
                    fout.write("autoScale on\n")
                    if bin_size == default_bin_size:
                        fout.write("visibility dense\n")
                    else:
                        fout.write("visibility hide\n")
                    #if strand == 'forward':
                    #    fout.write("\tcolor 0,0,255\n")
                    #elif strand == 'reverse':
                    #    fout.write("\tcolor 255,0,0\n")
                    fout.write(f"priority {p}\n")
                    fout.write("\n")
                    p += 1

            # Bin phase and amplitude bigBed
            for bin_size in [100,1000,10000]:
                for strand in Strands:
                    fout.write(f"track bin_phase_{Bin_size[bin_size]}_{strand[0]}\n")
                    fout.write("type bigBed 9\n")
                    fout.write("itemRgb on\n")
                    fout.write(f"shortLabel Bin phi {strand[0]} {Bin_size[bin_size]}\n")
                    fout.write(f"longLabel {strand} bin {Bin_size[bin_size]} phase and amplitude mapped in RGB space (red: 0h, yellow: 6h, green: 12h, blue: 18h)\n")
                    fout.write(f"bigDataUrl {track_folder}/phase_amp/bin_phase_amp_{strand}_{bin_size}bp.bb\n")
                    if bin_size == default_bin_size:
                        fout.write("visibility dense\n")
                    else:
                        fout.write("visibility hide\n")
                    #if strand == 'forward':
                    #    fout.write("\tcolor 0,0,255\n")
                    #elif strand == 'reverse':
                    #    fout.write("\tcolor 255,0,0\n")
                    fout.write(f"priority {p}\n")
                    fout.write("\n")
                    p += 1

            # bed tracks with extended kalman smoothing phase and amp on expressed regions
            for bin_size in [100,1000,10000]:
                for strand in Strands:
                    fout.write(f"track cont_ext_kalman_phase_amp_{Bin_size[bin_size]}_{strand[0]}\n")
                    fout.write("type bigBed 9\n")
                    fout.write("itemRgb on\n")
                    fout.write(f"shortLabel Kalman phi amp {strand[0]} {Bin_size[bin_size]}\n")
                    fout.write(f"longLabel Extended Kalman smoothing phase and amplitude {strand} strand {Bin_size[bin_size]}\n")
                    fout.write(f"bigDataUrl {track_folder}/kalman/extended_kalman_on_chromosomes_{strand}_bin{bin_size}bp_phi_amp.bb\n")
                    if bin_size == default_bin_size:
                        fout.write("visibility dense\n")
                    else:
                        fout.write("visibility hide\n")
                    #if strand == 'forward':
                    #    fout.write("\tcolor 0,0,255\n")
                    #elif strand == 'reverse':
                    #    fout.write("\tcolor 255,0,0\n")
                    fout.write(f"priority {p}\n")
                    fout.write("\n")
                    p += 1

            # bed tracks with extended kalman smoothing loglikelihood (transformed in 0-1000 range) on expressed regions
            for bin_size in [100,1000,10000]:
                for strand in Strands:
                    fout.write(f"track kalman_chromosomes_ll_{strand[0]}_{Bin_size[bin_size]}\n")
                    fout.write("type bigWig\n")
                    fout.write("itemRgb on\n")
                    fout.write(f"shortLabel Kalman LL {strand[0]} {Bin_size[bin_size]}\n")
                    fout.write(f"longLabel Extended Kalman smoothing on expressed regions {strand} strand {Bin_size[bin_size]}\n")
                    fout.write(f"bigDataUrl {track_folder}/kalman/extended_kalman_on_chromosomes_{strand}_bin{bin_size}bp_ll.bw\n")
                    if bin_size == default_bin_size:
                        fout.write("visibility dense\n")
                    else:
                        fout.write("visibility hide\n")
                    #if strand == 'forward':
                    #    fout.write("\tcolor 0,0,255\n")
                    #elif strand == 'reverse':
                    #    fout.write("\tcolor 255,0,0\n")
                    fout.write("autoScale on\n")
                    fout.write("\tminLimit 0\n")
                    fout.write("maxHeightPixels 100:30:8\n") # max:default:min
                    fout.write(f"priority {p}\n")
                    fout.write("\n")
                    p += 1

        
            # BigWig composite tracks with bin expression
            for bin_size in [100,1000,10000]:
                for strand in Strands:
                    fout.write(f"track PROseq_{strand}_{Bin_size[bin_size]}\n")
                    fout.write("compositeTrack on\n")
                    fout.write("subGroup1 t Time CT00=00 CT04=04 CT08=08 CT12=12 CT16=16 CT20=20 CT28=28 CT24=24 CT32=32 CT36=36 CT40=40 CT44=44\n")
                    fout.write("dimensions dimX=t\n")
                    fout.write("sortOrder t=+\n")
                    #fout.write("subGroup2 s Strand forward=forward reverse=reverse\n")
                    #fout.write("dimensions dimX=s dimY=t\n")
                    #fout.write("sortOrder s=+ t=+\n")
                    fout.write(f"shortLabel PROseq {strand[0]} {Bin_size[bin_size]}\n")
                    fout.write(f"longLabel PRO-seq data composite track {strand[0]} {Bin_size[bin_size]} (sum normed count per bin + 1, 1bp: norm count)\n")
                    fout.write("type bigWig\n")
                    if bin_size == default_bin_size:
                        fout.write("visibility full\n")
                    else:
                        fout.write("visibility hide\n")
                    fout.write("maxHeightPixels 100:30:8\n") # max:default:min
                    fout.write("autoScale group\n")
                    fout.write(f"descriptionUrl {args.track_hub_name}.html\n")
                    fout.write(f"priority {p}\n")
                    fout.write("\n")
                    p += 1
    
                    for sample in Samples:
                        name = sample
                        time = sample[2:]

                        fout.write(f"\ttrack {name}_{strand}_{Bin_size[bin_size]}\n")
                        fout.write(f"\tparent PROseq_{strand}_{Bin_size[bin_size]} on\n")
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

                        description_html = f"{genome}/{name}_{strand}.html"
                        fout.write(f"\tdescriptionUrl {description_html}\n")
                        fout.write(f"\n")

                        # make html description file
                        with open(f'{args.track_hub_folder}/{description_html}','w', encoding="utf-8") as fout2:
                            fout2.write(f"<h2>{name}_{strand}</h2>\n")
            
            # TAD
            fout.write(f"track TAD\n")
            fout.write("type bigBed 3\n")
            #fout.write("itemRgb off\n")
            fout.write(f"shortLabel TADs\n")
            fout.write(f"longLabel TADs\n")
            fout.write(f"bigDataUrl {track_folder}/tad/TADMap_scaffold_mm10.bb\n")
            fout.write("visibility squish\n")
            description_html = f"{genome}/tad.html"
            fout.write(f"\tdescriptionUrl {description_html}\n")
            fout.write("\n")

            # make description
            with open(f'{args.track_hub_folder}/{description_html}','w', encoding="utf-8") as fout2:
                fout2.write("<h2>TAD map from https://cb.csail.mit.edu/cb/tadmap/</h2>\n")
                fout2.write("\n")
