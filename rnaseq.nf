// 1 proces - przycinanie odczytów (trimming)
// usuwa adaptery i słabej jakości końcówki sekwencji, w celu poprawienia jakości mapowania

process TRIM_GALORE {
    publishDir "output/TRIMMED", mode:'copy'
    conda 'envs/powb-next.yml'

    input:
    	tuple val(sampleid), path(reads)

    output:
        path "*" //zwracanie wszystkich plików 
        path "*trimmed*.fq.gz", emit:trimmed //emitowanie plików .fq.gz do kanału 'trimmed'

    script:
        """
        # -q 20 -> prztnij jeśli Phred < 20
        # --paired -> tryb dla odczytów sparowanych
        # --gzip -> kompresja do .gz
        # --basename -> podstawa nazwy plików wyjściowych
        trim_galore -q 20 --paired -q 20 --gzip --basename ${sampleid}_trimmed ${reads}
        """
}

// 2 proces - kontrola jakości
// sprawdza jakość danych surowych i przyciętych oraz generuje raport HTML

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

// 3 proces - indeksowanie genomu

process STAR_INDEX {
    publishDir "output/INDEX", Mode:'copy'
    conda 'envs/powb-next.yml'
    
    input:
	    path(fasta) //sekwencja genomu
	    path(gtf)   // adnotacja genów

    output:
    	path "*", emit:index // emitowanie do kanału index

    script:
	    """
	    # --runMode genomeGenerate -> tryb generowania indeksu
	    # --genomeDir index -> folder wyjściowy
	    # --sjdbGTFfile -> plik z informacją o intronach/eksonach (GTF)
	    # --genomeSAindexNbases 12 -> parametr optymalizacji (12 jest dobre dla małych genomów)
	    STAR --runThreadN 8 \\
	    --runMode genomeGenerate \\
	    --genomeDir index \\
	    --genomeFastaFiles ${fasta} \\
	    --sjdbGTFfile ${gtf} \\
	    --genomeSAindexNbases 12
	    """
}

// 4 proces - mapowanie

process STAR_MAPPING {
    publishDir "output/MAPPING", mode:'copy'
    conda 'envs/powb-next.yml'
    cpus params.maxCpus

    input: 
    	tuple val(sampleid), path(read1), path(read2), path(index)

    output:
	    path "*"
	    path "*.bam", emit:bams // emitowane plików .bam do kanału bams

    script:
	    """
	    # --genomeDir -> ścieżka do indeksu stworzonego w poprzednim kroku
	    # --outSAMtype BAM SortedByCoordinate -> wynik od razu w formacie BAM, posortowany
	    # --readFilesCommand zcat -> komenda do czytania skompresowanych plików wejściowych (.gz)
	    STAR --runThreadN 4 --genomeDir ${index} \\
	    --readFilesIn ${read1} ${read2} \\
	    --outSAMtype BAM SortedByCoordinate \\
	    --outFileNamePrefix ${sampleid} \\
	    --readFilesCommand zcat
	    """
}

// 5 proces - zliczanie odczytów

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
	    # -T -> liczba wątków
	    # -p -> licz fragmenty (pary), a nie pojedyncze odczyty
	    # -t exon -> licz odczyty wpadające w eksony
	    # -g gene_id -> grupuj wyniki po ID genu
	    # -Q 10 -> ignoruj mapowania o niskiej jakości (<10)
	    featureCounts -T 2 -s ${strand} -p --countReadPairs -t exon \\
	    -g gene_id -Q 10 -a ${gtf} -o gene_count ${bams}

	    # generowanie raportu jakości zliczania
	    multiqc gene_count*
	    """
}

// proces 6 - wizualizacja za pomocą heatmapy

process VISUALIZE {
    publishDir "output/VISUALIZATION", mode: 'copy'
    conda 'envs/powb-python.yml'

    input:
    	path counts_file 

    output:
    	path "*.png"     

    script:
	    """
	    make_heatmap.py ${counts_file}
	    """
}

// 7 proces - analiza różnicowej ekspresji danych

process DIFFERENTIAL_EXPRESSION {

    publishDir "output/DIFFERENTIAL_EXPRESSION", mode: 'copy'
    conda 'envs/deseq2.yml' //jeśli nie działa (błąd 'no such command like Rscript') użyć środowiska deseq2_vol2.yml

    input:
        path(counts)    //zliczenia
        path(metadata) //przynależność do kategorii

    output:
        path "deseq2_results.tsv", emit: deseq_results
        path "*.png"

    script: //jak ma problemy bo nie widzi pliku - sprawdzić uprawnienia lub/i dać pełną ścieżkę do pliku
        """
        deseq2_analysis.R ${counts} ${metadata}
        """
}

// porces 

process FUNCTIONAL_ENRICHMENT {

    publishDir "output/FUNCTIONAL_ENRICHMENT", mode: 'copy'
    conda 'envs/clusterprofiler.yml'

    input:
        path(deseq_results)

    output:
        path "*.tsv"
        path "*.png"

    script:
        """
        functional_enrichment.R ${deseq_results}
        """
}


workflow {

    ref_fasta = channel.fromPath(params.ref_fasta)
    ref_gtf   = channel.fromPath(params.ref_gtf)
    
    // Wczytuje pary plików (np. sample_1.fq i sample_2.fq) i tworzy z nich krotki
    fastq_ch  = channel.fromFilePairs(params.reads)
    
    strand    = channel.of(params.strand)

    TRIM_GALORE(fastq_ch).set{trimmed}

    raw_fastq = fastq_ch.map{items -> items[1]}.flatten().collect()
    trimmed_fastq = trimmed.trimmed.flatten().collect()
    
    raw_fastq.mix(trimmed_fastq).collect() | QC

    STAR_INDEX(ref_fasta, ref_gtf).set{star_index}

    // bierzemy nazwę pliku i ucinamy końcówkę _trimmed aby odzyskać oryginalne ID próbki
    trimmed.trimmed.map{read1, read2 -> tuple("${read1.getFileName()}".split("_trimmed")[0], read1, read2) }
    | combine(star_index.index) // doklejanie indeksu do każdej próbki
    | STAR_MAPPING              // uruchamianie procesu mapowania
    | set{bams}                 // zapisanie w kanale bams

    bams.bams.collect().set{finalbams}

    fc_results = FEATURECOUNT(finalbams, ref_gtf, strand)

    fc_results
        .flatten()                          //rozbicie listy plików na pojedyncze elementy
        .filter { it.name == 'gene_count' } //szukanie pliku o konkretnej nazwie
        .set { count_table }                //zapisanie do kanału

    VISUALIZE(count_table)

    metadata = channel.fromPath(params.metadata) //wczytanie informacji na temat przynależności próbki do danej kategorii

    de_results = DIFFERENTIAL_EXPRESSION(
        count_table,
        metadata
    )

    FUNCTIONAL_ENRICHMENT(
        de_results.deseq_results
    )

}
