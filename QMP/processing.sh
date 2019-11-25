
### join paired end reads

#for f in fastqs/*_1.fastq
#do
#  base=`basename $f _1.fastq`
#  join_paired_ends.py -f $PWD/fastqs/$base"_1.fastq" -r $PWD/fastqs/$base"_2.fastq" -o $PWD/join/$base \
#      -j 140
#done


### collect ###
#for f in `find join/ -name *join.fastq`
#do
#   acc=$(basename `dirname $f`)
#   cp $f joinfastqs/$acc.fastq
#done

## convert to fasta ##
#parallel bash fastq2fasta.sh -- `find $PWD/joinfastqs -type f`

## add sample ids to fasta files
#add_qiime_labels.py -m mapping_file_check/mapping_file_corrected.txt -i fasta -c joinfile

# Pick OTUs command
#pick_otus.py -i /home/ubuntu/temp/QMP/combined_seqs.fna -o /home/ubuntu/temp/QMP/working_dir/sortmerna_picked_otus \
#     -r /dev/shm/gg_13_8_otus/rep_set/97_otus.fasta -m sortmerna --sortmerna_coverage 0.97 \
#     --otu_picking_method sortmerna --sortmerna_max_pos 10000 --similarity 0.97 --suppress_new_clusters --threads 3


#pick_rep_set.py -f /home/ubuntu/temp/QMP/combined_seqs.fna \
#	-i /home/ubuntu/temp/QMP/working_dir/sortmerna_picked_otus/combined_seqs_otus.txt \
#	-m most_abundant -o /home/ubuntu/temp/QMP/working_dir/rep_set.fasta

#assign_taxonomy.py -i /home/ubuntu/temp/QMP/working_dir/rep_set.fasta -r /dev/shm/gg_13_8_otus/rep_set/97_otus.fasta \
#                   -t /dev/shm/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -m sortmerna


# make_otu_table.py -i /home/ubuntu/temp/QMP/working_dir/sortmerna_picked_otus/combined_seqs_otus.txt \
#    -t /home/ubuntu/temp/QMP/working_dir/sortmerna_assigned_taxonomy/rep_set_tax_assignments.txt \
#    -o /home/ubuntu/temp/QMP/working_dir/pick_otus/otu_table.biom \
#    -m /home/ubuntu/temp/QMP/mapping_file_check/mapping_file_corrected.txt
