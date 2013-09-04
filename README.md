##### Contents
[Overview] (#overview)  
[Copy right] (#copyright)  
[How to cite BRANCH?] (#cite)  
[Short Manual] (#manual)  

<a name="overview"/>
### Overview
BRANCH is a software that extends de novo transfrags and identifies novel transfrags with DNA contigs or genes of close related species. BRANCH discovers novel exons first and then extends/joins fragmented de novo transfrags, so that the resulted transfrags are more complete.

<a name="copyright"/>
###Copy right
BRANCH is under the [Artistic License 2.0](http://opensource.org/licenses/Artistic-2.0).

<a name="cite"/>
### How to cite BRANCH?
If you use BRANCH, please cite the following paper:  
Bao E, Jiang T, Girke T (2013). BRANCH: boosting RNA-Seq assemblies with partial or related genomic sequences. Bioinformatics: [epub](http://bioinformatics.oxfordjournals.org/content/29/10/1250).

<a name="manual"/>
### Short manual
1. System requirements

   BRANCH is suitable for 32-bit or 64-bit machines with Linux operating systems. At least 4GB of system memory is recommended for assembling larger data sets. 

2. Installation

   The [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library is required to compile and run BRANCH.  
   The [BLAT](http://genome.ucsc.edu/FAQ/FAQblat.html) aligner is required to run BRANCH and the modified version (distributed with BRANCH) is highly recommended.  
   * Download the .cpp file.
   * If LEMON is already installed in your system, execute the command line: `g++ -o BRANCH BRANCH.cpp -lemon -lpthread`; otherwise, down load LEMON, compile it, and execute: `g++ -o BRANCH -I PATH2LEMON/include BRANCH.cpp -L PATH2LEMON/bin -lpthread`.
   * To use the modified BLAT, put it to your $PATH: `export PATH=PATH2BLAT:$PATH`.

3. Input
   * Single- or paired-end RNA reads in FASTA format.
   * De novo transfrags assembled by any de novo RNA assembler (Velvet/Oases, Trinity, etc.).
   * DNA contigs assembled by any de novo DNA assembler (Velvet, ABySS, etc.) or genome/gene sequences from a closely related species.

4. Using BRANCH

   ```
   BRANCH --read1 reads_1.fa --read2 reads_2.fa --transfrag transfrags.fa --contig contigs.fa --transcript transcripts.fa [--insertLow insertLow --insertHigh insertHigh --threshSize threshSize --threshCov threshCov --threshSplit threshSplit --threshConn threshConn --closeGap --noAlignment]
   ```

   Inputs:  
   --read1 is the first pair of PE RNA reads or single-end RNA reads in fasta format  
   --read2 is the second pair of PE RNA reads in fasta format  
   --transfrag is the de novo RNA transfrags to be extended  
   --contig is the reference DNA contigs  
   Output:  
   --transcript is the extended de novo transfrags  
   Options:  
   --insertLow is the lower bound of insert length (highly recommended; default: 0)  
   --insertHigh is the upper bound of insert length (highly recommended; default: 99999)  
   --threshSize is the minimum size of a genome region that could be identified as an exon (default: 2 bp)  
   --threshCov is the minimum coverage of a genome region that could be identified as an exon (default: 2)  
   --threshSplit is the minimum upstream and downstream junction coverages to split a genome region into more than one exons (default: 2)  
   --threshConn is the minimum connectivity of two exons that could be identified as a splice junction (default: 2)  
   --closeGap closes sequencing gaps using PE read information (default: none)  
   --noAlignment skips the initial time-consuming alignment step, if all the alignment files have been provided in tmp directory (default: none)

5. Output

   BRANCH outputs the transfrag file in FASTA format. It contains all the improved transfrags.

6. Important things to note
   * Single-end reads should have the same length and are not recommended, since the quality of single-end alignment is hard to be kept.
   * It is better to use related gene sequences rather than related genome sequences to greatly reduce run time and memory usage.
   * Though --insertLow and --insertHigh are options, they should always be specified to generate meaning result. Suppose the insert length is I, insertLow = I - 20 and insertHigh = I + 20 would be fine.
