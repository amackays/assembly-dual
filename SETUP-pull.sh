curl -O --output-dir scripts/ https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz
# not adding download instructions for fcs-gx since I'm not using it, but would be required along with its own database if desired
curl -O --output-dir scripts/ https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-adaptor.sif
singularity pull --dir scripts/ remiallio_default_mitofinder.sif library://remiallio/default/mitofinder
