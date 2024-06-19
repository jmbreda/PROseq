import os

if __name__ == '__main__':
    track_hub_name = 'TAD_map'

    genome = 'mm39'
    track_hub_folder = f'/data/web/sites/{track_hub_name}'
    if not os.path.exists(track_hub_folder):
        os.makedirs(track_hub_folder)
    track_hub_url = f'http://upnaesrv1.epfl.ch/{track_hub_name}'
    
    link_to_data = f'{track_hub_folder}/tracks_tad'
    data_folder = f'/bigdata/jbreda/PROseq/resources/TAD/'
    if not os.path.exists(link_to_data):
        os.system(f'ln -s {data_folder} {link_to_data}')

    # make track hub
    outfile=f'{track_hub_folder}/hub.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"hub {track_hub_name}\n")
        fout.write(f"shortLabel {track_hub_name}\n")
        fout.write("longLabel TAD from https://cb.csail.mit.edu/cb/tadmap/\n")
        fout.write("genomesFile genomes.txt\n")
        fout.write("email jeremie.breda@epfl.ch\n")
        fout.write(f"descriptionUrl {track_hub_name}.html\n")
        fout.write("\n")

    # make genomes.txt
    outfile=f'{track_hub_folder}/genomes.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"genome {genome}\n")
        fout.write(f"trackDb {genome}/trackDb.txt\n")
        fout.write("\n")

    # make {track_hub_name}.html
    outfile=f'{track_hub_folder}/{track_hub_name}.html'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("TAD map from https://cb.csail.mit.edu/cb/tadmap/\n")

    # make url.txt
    outfile=f'{track_hub_folder}/url.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"{track_hub_url}/hub.txt\n")

    # make trackDb.txt
    os.makedirs(f"{track_hub_folder}/{genome}", exist_ok=True)
    outfile=f"{track_hub_folder}/{genome}/trackDb.txt"
    with open(outfile,'w', encoding="utf-8") as fout:

        # Bed tracks with TADs
        fout.write(f"track TAD\n")
        fout.write("type bigBed 3\n")
        #fout.write("itemRgb off\n")
        fout.write(f"shortLabel TADs\n")
        fout.write(f"longLabel TADs\n")
        fout.write(f"bigDataUrl {track_hub_url}/tracks_tad/TADMap_scaffold_mm39.bb\n")
        fout.write("visibility squish\n")
        fout.write(f"\tdescriptionUrl {genome}/TAD.html\n")
        fout.write("\n")
    
    # make description
    outfile=f"{track_hub_folder}/{genome}/TAD.html"
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("TAD map from https://cb.csail.mit.edu/cb/tadmap/\n")
        fout.write("\n")