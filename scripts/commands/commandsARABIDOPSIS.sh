### Transfer results

#### Save CRAM files with sample name
cut -f1,2 -d "," input/arabidopsis_test_config.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/add_rb/${ID}*.cram CRAM/${SAMPLE}.cram  ; done 
cut -f1,2 -d "," input/arabidopsis_test_config.csv  | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/add_rb/${ID}*.cram.crai CRAM/${SAMPLE}.cram.crai  ; done 

#### Save NEAT CRAM files with sample name
cut -f1,2 -d "," input/arabidopsis_test_config.csv  | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.neat.cram NEAT_CRAM/${SAMPLE}.neat.cram  ; done 
cut -f1,2 -d "," input/arabidopsis_test_config.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.neat.cram.crai NEAT_CRAM/${SAMPLE}.neat.cram.crai  ; done 
