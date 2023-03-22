//
// create indices using bwa index and samtools faidx, depending on user-inps
//

include { GATK4_GENOMICSDBIMPORT } from '../../modules/nf-core/gatk4/genomicsdbimport/main'

workflow RUN_GENOMICSDBIMPORT{
    take:
    mergeGvcfTbi

    main:
    
    //here it is assumed that database folder name is same as that of the chromosome name --> for our group 
    //

    Channel
        .fromPath( params.regions )
        .splitText()
        .map{ it -> it.trim() }
        .set{ chrom }

    mergeGvcfTbiReg = mergeGvcfTbi.combine( chrom )


    updateMergeGvcfTbiReg = mergeGvcfTbiReg.map{ meta, vcf, tbi, chrom -> tuple( [id:chrom], vcf, tbi )}


    if ( params.genomic_db != null ){
        inp1_gatk4_genomicsdbimprt = updateMergeGvcfTbiReg.groupTuple().map{meta, vcf, tbi -> tuple(meta, vcf, tbi, [], meta.id, file(params.genomic_db+"/"+meta.id)) }
        run_updatewspace = true
        inp1_gatk4_genomicsdbimprt.view()
    }
    
    else{
       inp1_gatk4_genomicsdbimprt = updateMergeGvcfTbiReg.groupTuple().map{meta, vcf, tbi -> tuple(meta, vcf, tbi, [], meta.id, []) }
        run_updatewspace = false
    }

    GATK4_GENOMICSDBIMPORT(
        inp1_gatk4_genomicsdbimprt,
        false,
        run_updatewspace,
        false
        )
}
