
# Snake project

*Tutorial prepared by Deise Gonçalves in collaboration with Drew Larson and Peter Cerda*

Here we present all the steps to extract mitochondrion genes from transcriptome data for phylogenetic inference. Note that this tutorial can also be used to extract any other gene or sets of genes from different types of genomic datasets.

To run these analyses we will need the following software:
  
> [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) - used to trim the reads. 

> [Hybpiper](https://github.com/mossmatters/HybPiper) - used to extract sequences from the data providing a target file. 

> [MAFFT](https://mafft.cbrc.jp/alignment/software/) - used to align the sequences. This software has different versions for Mac OS X, Windows and Linux; download the appropriate version for your working machine.
  
> [Phyx](https://github.com/FePhyFoFum/phyx) - used to prepare a supermatrix and a partition file for phylogenetic inference. This software has multiple extremely helpful programs to deal with alignments and trees.  
  
> [IQ-TREE](http://www.iqtree.org) - a great tool for phylogenetic inference and conflict analyses. Software is easy to install and has great documentation.

You can add these software/scripts to the path or just give their absolute path when running the steps below.


The target file that we will use was prepared with mitogenomes of snakes that are available on GenBank (NCBI):
A fasta file with the sequences: ``Hybpiper_snake_targetfile.fa``  
A spreadsheet with information about GenBank data: ``Snake_existing_mitogenomes.xlsx``
  
> Quick summary of methods for making the target file
>> 1.	Downloaded the published mitochondria in the table from GenBank in GenBank format.
>> 2.	In Geneious Prime 2020.2.3, exported all annotations
>> 3.	Organized the sequences, reconciled name formatting differences, manually fixed a few issues, removed sequences for annotations that were not for rRNA, tRNA, or CDS.
>> 4.	Combined all sequences into a hybpiper format Target File
>> 5.	Arbitrarily renamed the two versions of Leucine and Serine tRNAs with As and Bs.

  
---

We will also use the following custom scripts to extract the sequences:
``snake_generate_trimmomatic_commands_PE.py``  -  this script will prepare the bash commands to run Trimmomatic
``combined_unpaired_dg.py`` - this script will combine unpaired reads
``generate_snake_hybpiper_commands.py``  -  this script will prepare the bash commands to run Hybpiper.   
*Note that you will have to change edit these scripts to customize the paths of your computer.*
  
We prepared two other scripts to rename trees and to extrac subsets of samples. Check them out at the end of this tutorial.


## Summary of how to run Trimmomatic and Hybpiper:
1. Put all the forward and reverse reads in a folder (they can be gzipped, whether they are or not has to be specified in the trimming script discussed below)  

2. Edit the generate read trimming script (``snake_generate_trimmomatic_commands_PE.py``) to point to the right folders and locations of Trimmomatic and the directory you want the trimmed reads in. You should also double check whether your files end with "fastq.gz" or not and change "ending" on your script. You can also adjust the trimming settings in this script. If you use nano or other text editor like:
```
nano snake_generate_trimmomatic_commands_PE_dg.py
```  

Check out what you will have to customize:
```
#This generates bash commands to run trimommatic for Paired end reads

#edit the line below with the path to the directory where you have your reads
reads_dir="/home/brlab/deisejpg/documents/snakes/0_raw/"

#edit the line below to add the path to the directory where you want your output files
output_dir="/home/brlab/deisejpg/documents/snakes/1_trimmomatic/"

#edit the line below if you have files with a different ending name
ending=".fastq.gz"
new_ending=".trim.fastq"
paired_end=True

#edit the lines belows with the path to the .jar file of Trimmomatic and the adaptor file which is within the Trimmomatic folder
trimmomatic_path="/home/brlab/deisejpg/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"
adaptor_file="/home/brlab/deisejpg/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
threads="30"
forward_pattern="_1."
reverse_pattern="_2."

import os,sys,re

#Example for single end reads
#java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Generates list of all forward reads to process in the reads directory
for fil in os.listdir(reads_dir):
        #Find only the files that are reads
        if fil[-len(ending):]==ending: #Matches if the extension is correct
                if forward_pattern in fil: #checks that it's a forward read so it doesn't get done twice
                        output_file_forward=output_dir+fil.replace(ending,new_ending)
                        output_file_reverse=output_file_forward.replace(forward_pattern,reverse_pattern)
                        forward_unpaired=output_file_forward+".unpaired"
                        reverse_unpaired=output_file_reverse+".unpaired"

                        if paired_end==True:
                                forward=reads_dir+fil
                                reverse=forward.replace(forward_pattern,reverse_pattern)
                                trim_cmd= "java -jar "+trimmomatic_path+" PE -threads "+threads+" "+forward+" "+reverse+" "+output_file_forward+" "+forward_unpaired+" "+output_file_reverse+" "+reverse_unpaired+" ILLUMINACLIP:"+adaptor_file+":2:30:10:2:TRUE SLIDINGWINDOW:5:20 LEADING:20 TRAILING:20 MINLEN:36"###Change this if you want to alter the settings

                        print(trim_cmd)
```  
  
3. Run the generate read trimming script, save the output as a bash script  
```
python snake_generate_trimmomatic_commands_PE.py &> run_trimmomatic_commands.sh

```
  
Each line of your bash script should have something similar to the following command:
```
java -jar /home/brlab/deisejpg/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 30 /home/brlab/deisejpg/documents/snakes/0_raw/Atco0910Liv_1.fastq.gz /home/brlab/deisejpg/documents/snakes/0_raw/Atco0910Liv_2.fastq.gz /home/brlab/deisejpg/documents/snakes/1_trimmomatic/Atco0910Liv_1.trim.fastq /home/brlab/deisejpg/documents/snakes/1_trimmomatic/Atco0910Liv_1.trim.fastq.unpaired /home/brlab/deisejpg/documents/snakes/1_trimmomatic/Atco0910Liv_2.trim.fastq /home/brlab/deisejpg/documents/snakes/1_trimmomatic/Atco0910Liv_2.trim.fastq.unpaired ILLUMINACLIP:/home/brlab/deisejpg/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:TRUE SLIDINGWINDOW:5:20 LEADING:20 TRAILING:20 MINLEN:36
ja
```  

4. Run the trimming bash script  
```
bash run_trimmomatic_commands.sh

#if you are running in a cluster, use the following to avoid  stopping the run
nohup bash run_trimmomatic_commands.sh &
```

5. After trimming run FASTQC on the trimmed files, may not want to do all of them since there will be at least 4x the number of samples. There will be overrepresented sequences, lots of duplicate sequences, etc. Mostly want to make sure we aren’t seeing adapters or overrepresented sequences that are Illumina indexes.  

6. Run the ``combined_unpaired_dg.py`` script in the folder with the trimmed reads (cd into the trimmed reads folder). You may need to edit this script to make the file names work properly with hybpiper, but it should be ok as is. You may or may not need to run this script, if the unpaired files are too small compared to the PE files then it wouldn't be a problem to skip this step and use only PE reads.

7. Edit the ``generate_snake_hybpiper_commands.py`` script to point to all the right paths and adjusting any settings just like you did for the Trimmomatic script.
```
import os,re

#adjust the path on the line below:
trimmed_reads="/home/brlab/deisejpg/documents/snakes/1_trimmomatic/"

#confirm the following matches the ending of your files:
ending=".fastq"
forward_pattern="_1.trim.fastq"
reverse_pattern=""

#edit the paths on the two lines below:
hybpiper_location="/home/brlab/deisejpg/bin/HybPiper/"
target_file="/home/brlab/deisejpg/documents/snakes/2_hybpiper/Hybpiper_snake_targetfile.fa"

for item in os.listdir(trimmed_reads):
    if forward_pattern == item[-len(forward_pattern):]:
            forward=trimmed_reads+item
            reverse=forward.replace("_1.","_2.")
            #name=re.findall("[A-Za-z]+_[0-9]+",item)[0]
            name=re.findall(".+(?=_1\.trim\.fastq)",item)[0]
            #unpaired=trimmed_reads+name+".combined_unpaired.fastq"
            hybpiper_cmd=hybpiper_location+"reads_first.py --cov_cutoff 8 -r "+forward+" "+reverse+" -b "+target_file+" --prefix "+name+" --timeout 1000 --ins_length 250 --bwa"
            print hybpiper_cmd
            print "mv "+name+" hybpiper_output"
#            print "rm "+forward+" "+reverse+" "+unpaired
```
  
8. Run generate_snake_hybpiper_commands.py, save the output as a bash script, and run the bash script 
```
#generating the commands:
python generate_snake_hybpiper_commands.py &> run_hybpiper_commands.sh

#running the bash script:
bash run_hybpiper_commands.sh
#or
nohup bash run_hybpiper_commands.sh &
```
    
For each sample your bash script should look like this:
```
/home/brlab/deisejpg/bin/HybPiper/reads_first.py --cov_cutoff 8 -r /home/brlab/deisejpg/documents/snakes/1_trimmomatic/Atco0910Liv_1.trim.fastq /home/brlab/deisejpg/documents/snakes/1_trimmomatic/Atco0910Liv_2.trim.fastq -b /home/brlab/deisejpg/documents/snakes/2_hybpiper/Hybpiper_snake_targetfile.fa --prefix Atco0910Liv --timeout 1000 --ins_length 250 --bwa
mv Atco0910Liv hybpiper_output
```

9.Retrieve hybpiper sequences:  
```
#change the path to the python script to reflect the path on your machine
#put the target file on the current directory or add the path to it
#the "." means your are pointing to the current directory it means that 
#you are working on the directory where you have all hybpiper output files
python ~/apps/HybPiper/retrieve_sequences.py ../Hybpiper_snake_targetfile.fa . dna

python /home/brlab/deisejpg/bin/HybPiper/retrieve_sequences.py path_to_target/Hybpiper_snake_targetfile.fa . dna
```
  

10. Generate a file with assembly stats. It is important to detect possible problems with the sequences or with the runs. To do this do something like:  

>1.	Generate a file that has the average lengths of the gene targets:  
```
python ~/apps/HybPiper/get_seq_lengths.py ../Hybpiper_snake_targetfile.fa length_of_gene_targets.txt
```  

>2. Generate a list of the samples:  
```
ls > list_of_samples.txt
```  

>3. Run  
```
python ~/apps/HybPiper/hybpiper_stats.py length_of_gene_targets.txt list_of_samples.txt &> stats.txt
```
Depending on the number of samples it may be easier to split up to the samples and then concatenate the results.


>11. Add the published sequences back into their respective fasta files to include in the downstream analyses. You can either prepare the outgroup samples from GenBank by hand or you can use the script ``extracting_subset.py``. Check out info about this script below in *Tips*. After preparing files with all outgroup samples for each gene of interest, you can use ``cat`` to concatenate ingoup and outgroup species.
```
cat sample_outgroup.fa sample_subsetIngroup.fasta >> sample_rRNA_cat.fasta

#for loop to generate bash commands:
for file in *outgroup.fa; do echo cat $file ${file//_outgroup.fa/_subsetIngroup.fasta} '>>' ${file//_outgroup.fa/_cat.fasta}; done &> run_cat.sh
```  

>12. Now we have sequences to align and infer phylogenies.

## Aligning the sequences and inferring the phylogenies using IQ-TREE
  
1. Align gene sequences using MAFFT
```
#use the following loop to prepare the code for all your gene sequences
mkdir 1_mafft_in
cd mafft_in
mkdir ../2_mafft_out
for file in *fasta; do echo mafft --thread 12 --maxiterate 1000 --ep 0.123 --genafpair $file '>>' ../2_mafft_out/${file//_subset.fasta/.aln.fa}; done &> run_mafft.sh

# run_mafft.sh should have lines looking like this:
mafft --thread 12 --maxiterate 1000 --ep 0.123 --genafpair sample_rRNA_cat.fasta >> ../2_mafft_out/sample_rRNA.aln.fa

#now run mafft:
bash mafft.sh
```

Once you have your sequences aligned, use phix to concatenate the gene sequences into a supermatrix keeping information about the gene partition.
```
#the code below considers you are in the directory where you saved your aligned sequences
#Parts.txt will have the partition in the format for IQ-TREE
pxcat -s *fa -p Parts.txt -o Supermatrix.fa
```

2. Infer phylogenies using IQ-TREE  
```
iqtree2 -s Supermatrix.fa -p Parts.txt --prefix concat_notrim -B 1000 -m TEST -mset raxml
```

  
## *Tips:*
Depending on the number of samples you have, run Trimmomatic and Hybpiper for subsets of samples and then combine the results at the end. Both software can generate really big files leading the computer to run out of storage.

If you need to rename your trees or if you want to use just a subset of the data extracted after running Hybpiper, you can easily do that using the scripts ``rename_snake_tips.py`` and ``extracting_subset.py``.


To run ``rename_snake_tips.py`` prepare a csv file with two columns, the first with 'old' names and the second with 'new' names, as follows:

AB008539	Lycodon_semicarinatus_AB008539   
DQ523162	Pantherophis_slowinskii_DQ523162  
DQ343648	Naja_naja_DQ343648               


To rename the trees you will also need the tree file! Use the code below to run python:
```
python rename_snake_tips.py sample_names.csv oldNames.tre newNames.tre
```

To extract a subset of samples, prepare a csv file with one sample name per line. Use the following code:
```
python extracting_subset.py fullFasta_file.fasta csv_file.csv subsetFasta_file.fasta
```

