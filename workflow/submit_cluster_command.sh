#!/bin/bash
snakemake -s workflow/Snakefile \
          -j 500 \
          --use-conda \
          --cluster-config config/cluster.json \
          --cluster "sbatch --job-name {cluster.name} \
                            --qos {cluster.qos} \
                            --time {cluster.time} \
                            --mem {cluster.mem} \
                            --nodes {cluster.nodes} \
                            --ntasks {cluster.ntasks} \
                            --cpus-per-task {cluster.cpus-per-task} \
                            --output {cluster.stdout} \
                            --error {cluster.stderr}" 
                            #--rulegraph | dot -Tpdf > workflow/rulegraph.pdf
                            # --dag | dot -Tpdf > workflow/dag.pdf
                            
                            