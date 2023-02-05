//
// create indices using bwa index and samtools faidx, depending on user-inputs
//

include { GATK4_HAPLOTYPECALLER } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_HAPLOTYPECALLER_REG } from '../../modules/local/gatk/haplotypecaller_reg/main'

workflow CREATE_GVCF{
    take:
        tuple_meta_bam
        fastaF
        fai 
        dict

    main:
    if ( !params.skip_regionwise_gvcf ){    
        Channel
            .fromPath( params.regions )
            .splitText()
            .map{ it -> it.trim() }
            .set{ reg }
        bam_region = tuple_meta_bam.combine(reg)
        bam_region.view()
        }
    else{
        regionF = Channe.fromPath( params.regions )
        bam_region = tuple_meta_bam.combine( regionF )
        }

    bam_fasta = bam_region.combine( fastaF )
    bam_fasta_fai = bam_fasta.combine( fai )
    bam_fasta_fai_seqdict = bam_fasta_fai.combine( dict )

    input_arg1 = bam_fasta_fai_seqdict.map{ meta, bam, idx, rg, fas, fai, dict -> tuple( meta, bam, idx, rg, [] ) }
    input_arg2 = bam_fasta_fai_seqdict.map{ meta, bam, idx, rg, fas, fai, dict -> fas }
    input_arg3 = bam_fasta_fai_seqdict.map{ meta, bam, idx, rg, fas, fai, dict -> fai }
    input_arg4 = bam_fasta_fai_seqdict.map{ meta, bam, idx, rg, fas, fai, dict -> dict }

    input_arg1.view()


    //
    //MODULE GATK4_HAPLOTYPECALLER

    if ( !params.skip_regionwise_gvcf ){    
        GATK4_HAPLOTYPECALLER_REG(
            input_arg1,
            input_arg2,
            input_arg3,
            input_arg4,
            [],
            []
            )
        }
    else{
        GATK4_HAPLOTYPECALLER(
            input_arg1,
            input_arg2,
            input_arg3,
            input_arg4,
            [],
            []
            )
        }
}
