//
// create indices using bwa index and samtools faidx, depending on user-inputs
//

include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX   } from '../../modules/nf-core/samtools/faidx/main'

workflow FASTA_INDICES{
    take:
        fastaF
    main:
        if ( !params.skip_picard_seqdict ){
                fasta_meta_pre = fastaF.map{ it -> tuple([id:it.baseName], it)}
                PICARD_CREATESEQUENCEDICTIONARY( 
                    fasta_meta_pre
                    )
                picard_idx = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict.map{metaP, faiP -> faiP}
            }
        else{
                picard_idx = Channel.fromPath( params.picard_seqdict )
            }
        if ( !params.skip_samtools_faidx ){
                def meta_sam = [:]
                meta_sam.id = "assembly"
                meta_sam.single_end = false
                faidx_meta_pre = fastaF.map{ file -> tuple(meta_sam, file)  }
                SAMTOOLS_FAIDX( 
                    faidx_meta_pre
                    )
                fa_idx = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> fai}
            }
        else{
                fa_idx = Channel.fromPath( params.samtools_faidx)
            }
        
    
    emit:
        picard_idx                                     // channel: [ val(meta), [ bwa_idx_dir ] ]
        fa_idx                              // channel: path(fa_idx)
        picard_seqdict_version = params.skip_picard_seqdict?' ': PICARD_CREATESEQUENCEDICTIONARY.out.versions.first()
        samtools_faidx_version = params.skip_samtools_faidx?' ': SAMTOOLS_FAIDX.out.versions.first()
}
