
=======
# LAB 5 
We start with an already made BLAST database. The first thing to do is to download the protein sequence the protein of interest, XP_001621078.2.

    ncbi-acc-download -F fasta -m protein XP_001621078.2

Next, perform a BLAST search with the CELSR2 protein. This will align our query sequence with each sequence in the database to find and identify similar sequences (regions of local similarity).

    blastp -db ~/data/blast/allprotein.fas -query XP_001621078.2.fa -outfmt 0 -max_hsps 1 > XP_001621078.2.blastp.typical.out

This next command will specify an output format 

    blastp -db ~/data/blast/allprotein.fas -query XP_001621078.2.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out XP_001621078.2.blastp.detail.out

To filter the BLAST output for putative homologs of the protein, we use the awk command and set the e-value threshold to be 1e-46, meaning any homologs with E-values at or lower than 1e-46 will be included.

    awk '{if ($3>=260 && $6<0.0000000000000000000000000000000000000000000001)print $1 }' ~/labs/lab5/lab5-ellieteng/XP_001621078.2.blastp.detail.out > ~/labs/lab5/lab5-ellieteng/XP_001621078.2.blastp.detail.filtered.out

To see the set of putative homologs, we can use this command:

    wc -l XP_001621078.2.blastp.detail.filtered.out

This seqkit grep command will yield a file of the CELSR2 protein, its 19 putative homologs, and their amino acid sequences. The name of this file is "XP_001621078.2.blastp.detail.filtered.fas" 

    seqkit grep --pattern-file XP_001621078.2.blastp.detail.filtered.out ~/data/blast/allprotein.fas > XP_001621078.2.blastp.detail.filtered.fas

Looking into the file, we can see that the names of the proteins are long and some are missing the full name of the species. I shortened these definition lines and added the full species names to the ones lacking them by using:

    nano XP_001621078.2.blastp.detail.filtered.fas

The muscle command is used to perform a global multiple sequence alignment:

    muscle -in XP_001621078.2.blastp.detail.filtered.fas -out XP_001621078.2.blastp.detail.filtered.aligned.fas

The t_coffee command is used to view statistics about the alignment:

    t_coffee -other_pg seq_reformat -in XP_001621078.2.blastp.detail.filtered.aligned.fas -output sim

We can also view the alignment using alv. In the second alv command, the -majority tag will only color the column(s) where the most common amino acids are found in 50% of the sequences:

    alv -k XP_001621078.2.blastp.detail.filtered.aligned.fas | less -r

    alv -kli --majority XP_001621078.2.blastp.detail.filtered.aligned.fas | less -RS

The t_coffee command is used again, this time to remove highly gapped positions in the alignment. The resulting FASTA file is named "myhomologs.aligned.r50.fa":

    t_coffee -other_pg seq_reformat -in XP_001621078.2.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out myhomologs.aligned.r50.fa

We can use the alv command again to view the length of the alignment after highly gapped positions have been removed:

    alv -kli --majority myhomologs.aligned.r50.fa | less -RS

# LAB 6

We start with the Newick Utilities software already installed. To get started, copy the filtered & aligned file from lab5:

    cp ~/labs/lab5/lab5-ellieteng/myhomologs.aligned.r50.fa .
    
Use the software IQ-TREE to construct a phylogenetic tree from sequence data:

    iqtree -s myhomologs.aligned.r50.fa -nt 2

To view a summary of the substitution model of the tree, look into the iqtree file:

    less myhomologs.aligned.r50.fa.iqtree

The unrooted tree can be viewed in the iqtree file (the command above), but the nw_display command will display only the tree. This unrooted tree is shown with an arbitrary root:

    nw_display myhomologs.aligned.r50.fa.treefile

gotree is used to view the phylogeny unrooted: 

    gotree draw png -w 1000 -i myhomologs.aligned.r50.fa.treefile  -r -o  myhomologs.aligned.r50.fa.png

Two types of rooting can be performed. This first one will be a midpoint root, where the root is specified at the midpoint of the longest branch.

    gotree reroot midpoint -i myhomologs.aligned.r50.fa.treefile -o myhomologs.aligned.r50.fa.midpoint.treefile

To view the midpoint rooted tree: 

    nw_display   myhomologs.aligned.r50.fa.midpoint.treefile

    nw_display -s  myhomologs.aligned.r50.fa.midpoint.treefile -w 1000 -b 'opacity:0' >  myhomologs.aligned.r50.fa.midpoint.svg
    
The second type of rooting is a literature based outgroup root. From the literature, I have deducted that Nematostella_XP_032228413.1, an uncharacterized protein, is the oldest/least related, however, the evidence is not concrete so this root will not be used in the final analysis of the CELSR2 query protein. Still, it is interesting to see that this outgroup root is very similar to the midpoint root.

    nw_reroot myhomologs.aligned.r50.fa.treefile Nematostella_XP_032228413.1 >myhomologs.aligned.r50.fa.MendozaRoot.treefile

To view the outgroup rooted tree:

    nw_display -s  myhomologs.aligned.r50.fa.MendozaRoot.treefile -w 1000 -b 'opacity:0' >  myhomologs.aligned.r50.fa.MendozaRoot.svg

A cladogram can be created with this command. This may be useful when compared to the previous phylograms:

    nw_topology myhomologs.aligned.r50.fa.MendozaRoot.treefile | nw_display -s  -w 1000 > myhomologs.aligned.r50.fa.MendozaRoot.cladogram.svg -

# LAB 7

Copy the FASTA file and the treefile from the previous lab to this lab:

    cp ~/labs/lab6/lab6-ellieteng/myhomologs.aligned.r50.fa .
    cp ~/labs/lab6/lab6-ellieteng/myhomologs.aligned.r50.fa.treefile .

Obtain bootstrap support values by running an ultrafast bootstrap. 

    iqtree -s myhomologs.aligned.r50.fa -bb 1000 -nt 2 --prefix myhomologs.r50.ufboot

Midpoint root the IQ-TREE output:

    gotree reroot midpoint -i myhomologs.r50.ufboot.treefile -o myhomologs.aligned.r50.fa.midpoint.treefile

View bootstrap support:

    nw_display   myhomologs.aligned.r50.fa.midpoint.treefile
    nw_display -s  myhomologs.aligned.r50.fa.midpoint.treefile -w 1000 -b 'opacity:0' >  myhomologs.aligned.r50.fa.midpoint.svg

To reconcile the gene tree with the species tree, first make a batch file in nano:

    nano mygenebatch.txt
    
On the first line, put the name of species tree (species.tre) & on the second line, put the name of the gene tree (myhomologs.aligned.r50.fa.midpoint.treefile)

Use Notung to reconcile the gene and species tree:

    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b mygenebatch.txt --reconcile --speciestag prefix  --savepng --treestats --events  --phylogenomics

To view the gene within species tree, generate a RecPhyloXML. This will be pushed into the GitHub repository where its contents will be pasted into http://phylariane.univ-lyon1.fr/recphyloxml/recphylovisu which will display a tree:

    python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g myhomologs.aligned.r50.fa.midpoint.treefile.reconciled --include.species

These commands will allow us to analyze the gene copies and their changes:

    column -t -s $'\t' mygenebatch.txt.species.tre.geneCount.txt | less -S
    column -t -s $'\t' mygenebatch.txt.species.tre.geneGainLoss.txt | less -S

Use Notung again to root the tree to minimize the number of duplication and loss events:

    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b mygenebatch.txt --root --speciestag prefix  --savepng --treestats --events  --phylogenomics 

# LAB 8

First, copy the necessary files from the previous lab(s) into the current lab directory:

    cp ~/labs/lab7/lab7-ellieteng/myhomologs.aligned.r50.fa.midpoint.treefile .
    cp ~/labs/lab7/lab7-ellieteng/myhomologs.aligned.r50.fa .
    cp ~/labs/lab7/lab7-ellieteng/myhomologs.r50.ufboot.treefile .
    cp ~/labs/lab7/lab7-ellieteng/species.tre .
    cp ~/labs/lab7/lab7-ellieteng/mygenebatch.txt .

Reconcile and rearrange the tree. We will place the output in a new directory named "mygenereconcileRearrange":

    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b mygenebatch.txt --rearrange --speciestag prefix  --savepng --treestats --events  --outputdir mygenereconcileRearrange --edgeweights name --threshold 90 

Move to the "mygenereconcileRearrange" directory:

    cd mygenereconcileRearrange

This command will allow us to view the rearranged tree within the species tree:

    python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g myhomologs.aligned.r50.fa.midpoint.treefile.rearrange.0 --include.species
    
Change back to the previous directory (lab8-ellieteng):

    cd ../

or

    cd ~/labs/lab8/lab8-ellieteng

Change the format of the rearranged tree from Notung to Newick.

    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -g mygenereconcileRearrange/myhomologs.aligned.r50.fa.midpoint.treefile.rearrange.0   -s species.tre --reconcile --speciestag prefix  --treeoutput newick --nolosses

Unroot the Notung rearranged tree: 

    gotree unroot -i myhomologs.aligned.r50.fa.midpoint.treefile.rearrange.0.reconciled -o myhomologs.r50.ufboot.unrooted.treefile.rearrange

The two trees in question are the optimal tree from IQ-TREE (myhomologs.r50.ufboot.treefile) and the Notung rearranged tree found in the  "mygenereconcileRearrange" directory (myhomologs.r50.ufboot.unrooted.treefile.rearrange.0). The Notung rearranged tree must first be unrooted to be concatenated (this was done in the previous command). Concatenate these alternative trees:

    cat myhomologs.r50.ufboot.treefile myhomologs.r50.ufboot.unrooted.treefile.rearrange > myhomologs.r50.alternativetrees

Use IQ-TREE to run the topology test: 

    iqtree -s myhomologs.aligned.r50.fa -z myhomologs.r50.alternativetrees -au -zb 10000 --prefix mygene_altTrees -m LG+F+R5 -nt 2 -te myhomologs.r50.ufboot.treefile

To look at the results of the AU topology test, look in mygene_altTrees.iqtree:

    less mygene_altTrees.iqtree
    
# LAB 9

The interproscan web service is used to find pre-defined domains in our CELSR2 gene. The unaligned protein sequence is used to run iprscan5:

    iprscan5   --email ellie.teng@stonybrook.edu  --multifasta --useSeqId --sequence   XP_001621078.2.blastp.detail.filtered.fas

Concatenate the domain files (ending in .tsv.txt) for each gene. These will be compiled into a new file "my.domains.all.tsv":

    cat Nematostella_vectensis_XP_001621078_2.tsv.txt Nematostella_vectensis_XP_001622990_2.tsv.txt Nematostella_vectensis_XP_001641262_2.tsv.txt Nematostella_vectensis_XP_032221363_1.tsv.txt Nematostella_vectensis_XP_032224204_1.tsv.txt Nematostella_vectensis_XP_032227744_1.tsv.txt Nematostella_vectensis_XP_032228413_1.tsv.txt Nematostella_vectensis_XP_032231912_1.tsv.txt Nematostella_vectensis_XP_032242051_1.tsv.txt Nematostella_vectensis_XP_032242104_1.tsv.txt Homo_sapiens_AGRL2.tsv.txt Homo_sapiens_AGRL3.tsv.txt Homo_sapiens_CELSR2.tsv.txt Homo_sapiens_CELSR1.tsv.txt Homo_sapiens_CELSR3.tsv.txt Strongylocentrotus_purpuratus_SPU_003929.tsv.txt Strongylocentrotus_purpuratus_SPU_009215.tsv.txt Pocillopora_damicornis_pdam_00013144.tsv.txt Pocillopora_damicornis_pdam_00018939.tsv.txt > my.domains.all.tsv

Use the grep command to filter for domains defined by the Pfam database. The output will be in the file "my.domains.pfam.tsv": 

    grep Pfam my.domains.all.tsv >  my.domains.pfam.tsv

Rearrange the interproscan output:

    awk 'BEGIN{FS="\t"} {print $1"\t"$3"\t"$7"@"$8"@"$5}' my.domains.pfam.tsv | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g' > my.domains.pfam.evol.tsv




>>>>>>> 69e3b36deb8baa4e19da39988ae9e46967766827
# bio312-finalproject

