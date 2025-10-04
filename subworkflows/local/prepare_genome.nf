//
// Uncompress and prepare reference genome files
//

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_BLACKLIST } from '../../modules/nf-core/modules/gunzip/main'

include {
    UNTAR as UNTAR_BOWTIE2_INDEX } from '../../modules/nf-core/modules/untar/main'

include { GFFREAD                } from '../../modules/nf-core/modules/gffread/main'
include { CUSTOM_GETCHROMSIZES   } from '../../modules/nf-core/modules/custom/getchromsizes/main'
include { BOWTIE2_BUILD          } from '../../modules/nf-core/modules/bowtie2/build/main'
include { GTF2BED                } from '../../modules/local/gtf2bed'
include { GENOME_BLACKLIST_REGIONS } from '../../modules/local/genome_blacklist_regions'

workflow PREPARE_GENOME {
    take:
    prepare_tool_index // string  : tool to prepare index for
    mito_name          // string  : mitochondrial chromosome name (e.g., 'MT', 'chrM', 'Mt')

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map{ it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(params.fasta))
    } 

    // // Make fasta file available if reference saved or IGV is run
    // if (params.save_reference || !params.skip_igv) {
    //     file("${params.outdir}/genome/").mkdirs()
    //     ch_fasta.copyTo("${params.outdir}/genome/")
    // }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [:], params.gff ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf      = GFFREAD ( ch_gff ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress blacklist file if required
    //
    ch_blacklist = Channel.empty()
    if (params.blacklist) {
        if (params.blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST ( [ [:], params.blacklist ] ).gunzip.map{ it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_BLACKLIST.out.versions)
        } else {
            ch_blacklist = Channel.fromPath(file(params.blacklist))
        }
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //

    // If --gtf is supplied along with --genome
    // Make gene bed from supplied --gtf instead of using iGenomes one automatically
    def make_bed = false
    if (!params.gene_bed) {
        make_bed = true
    } else if (params.genome && params.gtf) {
        if (params.genomes[ params.genome ].gtf != params.gtf) {
            make_bed = true
        }
    }

    if (make_bed) {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    } else {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], params.gene_bed ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    }

    //
    // Create chromosome sizes file
    //
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES ( ch_fasta ).sizes
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Prepare genome intervals for filtering by removing regions in blacklist file
    //
    ch_genome_filtered_bed = Channel.empty()

    GENOME_BLACKLIST_REGIONS (
        CUSTOM_GETCHROMSIZES.out.sizes,
        ch_blacklist.ifEmpty([]),
        mito_name,
        params.keep_mito,
        params.keep_blacklist
    )
    ch_genome_filtered_bed = GENOME_BLACKLIST_REGIONS.out.bed
    ch_versions = ch_versions.mix(GENOME_BLACKLIST_REGIONS.out.versions)


    //
    // Uncompress Bowtie2 index or generate from scratch if required
    //
    ch_bowtie2_index = Channel.empty()
    if (prepare_tool_index == 'bowtie2') {
        if (params.bowtie2_index) {
            if (params.bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( [ [:], params.bowtie2_index ] ).untar.map { meta, index -> index }
                ch_versions      = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
            } else {
                ch_bowtie2_index = Channel.value(file(params.bowtie2_index))
            }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta.map { fasta -> [ [:], fasta ] } ).index.map { meta, index -> index }
            ch_versions      = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }

    emit:
    fasta          = ch_fasta                             //    path: genome.fasta
    gtf            = ch_gtf                               //    path: genome.gtf
    gene_bed       = ch_gene_bed                          //    path: gene.bed
    chrom_sizes    = ch_chrom_sizes                       //    path: genome.sizes
    filtered_bed   = ch_genome_filtered_bed               //    path: *.include_regions.bed
    blacklist      = ch_blacklist                         //    path: blacklist.bed (original)
    bowtie2_index  = ch_bowtie2_index                     //    path: bowtie2/index/

    versions       = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}