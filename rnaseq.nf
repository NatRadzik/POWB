
process TRIM_GALORE {
    publishDir "output/TRIMMED", mode:'copy'
    conda 'envs/powb-next.yml'

    input:
    	tuple val(sampleid), path(reads)

    output:
        path "*"
        path "*trimmed*.fq.gz", emit:trimmed

    script:
        """
        trim_galore -q 20 --paired -q 20 --gzip --basename ${sampleid}_trimmed ${reads}
        """
}

process QC{

    publishDir "output/QC_REPORT", mode:'copy'
    conda 'envs/powb-next.yml'

    input:
    	path(reads)

    output:
    	path "*"

    script:
	    """
	    # Uruchomienie FastQC dla każdego pliku wejściowego
	    fastqc ${reads}

	    # MultiQC zbiera wyniki z FastQC i tworzy jeden zbiorczy raport
	    multiqc *fastqc*

	    # Porządki: przeniesienie pojedynczych raportów do podkatalogu
	    mkdir FASTQC
	    mv *fastqc* FASTQC
	    """
}

process STAR_INDEX {
    publishDir "output/INDEX", Mode:'copy'
    conda 'envs/powb-next.yml'
    
    input:
	    path(fasta) 
	    path(gtf)   

    output:
    	path "*", emit:index

    script:
	    """
	    # Tworzenie indeksu algorytmem STAR:
	    # --runMode genomeGenerate: Tryb generowania indeksu
	    # --genomeDir index: Folder wyjściowy
	    # --sjdbGTFfile: Plik z informacją o intronach/eksonach (GTF)
	    # --genomeSAindexNbases 12: Parametr optymalizacji (12 jest dobre dla małych genomów)
	    STAR --runThreadN 8 \\
	    --runMode genomeGenerate \\
	    --genomeDir index \\
	    --genomeFastaFiles ${fasta} \\
	    --sjdbGTFfile ${gtf} \\
	    --genomeSAindexNbases 12
	    """
}

process STAR_MAPPING {
    publishDir "output/MAPPING", mode:'copy'
    conda 'envs/powb-next.yml'

    cpus params.maxCpus

    input: 
    	tuple val(sampleid), path(read1), path(read2), path(index)

    output:
	    path "*"
	    path "*.bam", emit:bams

    script:
	    """
	    # Mapowanie odczytów:
	    # --genomeDir: Ścieżka do indeksu stworzonego w poprzednim kroku
	    # --outSAMtype BAM SortedByCoordinate: Wynik od razu w formacie BAM, posortowany
	    # --readFilesCommand zcat: Komenda do czytania skompresowanych plików wejściowych (.gz)
	    STAR --runThreadN 4 --genomeDir ${index} \\
	    --readFilesIn ${read1} ${read2} \\
	    --outSAMtype BAM SortedByCoordinate \\
	    --outFileNamePrefix ${sampleid} \\
	    --readFilesCommand zcat
	    """
}


process FEATURECOUNT {
    publishDir "output/FEATURECOUNT", mode:'copy'
    conda 'envs/powb-next.yml'

    input:
	    path(bams)   
	    path(gtf)    
	    val(strand)  

    output:
    	path "*"    

    script:
	    """
	    # Zliczanie programem featureCounts:
	    # -T 8: Użyj 8 wątków
	    # -p: Licz fragmenty (pary), a nie pojedyncze odczyty
	    # -t exon: Licz odczyty wpadające w eksony
	    # -g gene_id: Grupuj wyniki po ID genu
	    # -Q 10: Ignoruj mapowania o niskiej jakości (<10)
	    featureCounts -T 2 -s ${strand} -p --countReadPairs -t exon \\
	    -g gene_id -Q 10 -a ${gtf} -o gene_count ${bams}

	    # Generuj raport jakości zliczania
	    multiqc gene_count*
	    """
}

process VISUALIZE {
    publishDir "output/VISUALIZATION", mode: 'copy'
    conda 'envs/powb-python.yml'

    input:
    	path counts_file 

    output:
    	path "*.png"     

    script:
	    """
	    scripts/make_heatmap.py ${counts_file}
	    """
}

process DIFFERENTIAL_EXPRESSION {

    publishDir "output/DIFFERENTIAL_EXPRESSION", mode: 'copy'
    conda 'envs/deseq2.yml'

    input:
        path(counts)
        path(metadata)

    output:
        path "deseq2_results.tsv"
        path "*.png"

    script:
        """
        Rscript scripts/deseq2_analysis.R ${counts} ${metadata}
        """
}

// --- GŁÓWNY PRZEPŁYW DANYCH (WORKFLOW) ---
workflow {

    ref_fasta = Channel.fromPath(params.ref_fasta)
    ref_gtf   = Channel.fromPath(params.ref_gtf)
    
    // Wczytuje pary plików (np. sample_1.fq i sample_2.fq) i tworzy z nich krotki
    fastq_ch  = Channel.fromFilePairs(params.reads)
    
    strand    = Channel.of(params.strand)

    TRIM_GALORE(fastq_ch).set{trimmed}


    raw_fastq = fastq_ch.map{items -> items[1]}.flatten().collect()
    trimmed_fastq = trimmed.trimmed.flatten().collect()
    
    raw_fastq.mix(trimmed_fastq).collect() | QC

    STAR_INDEX(ref_fasta, ref_gtf).set{star_index}

    trimmed.trimmed.map{read1, read2 -> tuple("${read1.getFileName()}".split("_trimmed")[0], read1, read2) }
    | combine(star_index.index) 
    | STAR_MAPPING              
    | set{bams}                 

    bams.bams.collect().set{finalbams}

    fc_results = FEATURECOUNT(finalbams, ref_gtf, strand)

    fc_results
        .flatten()                          
        .filter { it.name == 'gene_count' }
        .set { count_table }                

    //VISUALIZE(count_table)
    
    metadata = Channel.fromPath(params.metadata)

    DIFFERENTIAL_EXPRESSION(
        count_table,
        metadata
    )

}
