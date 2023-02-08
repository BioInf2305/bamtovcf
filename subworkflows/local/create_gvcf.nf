//
// create indices using bwa index and samtools faidx, depending on user-inps
//

include { GATK4_HAPLOTYPECALLER } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK_HAPLOTYPECALLER_REG } from '../../modules/local/gatk/haplotypecaller_reg/main'
include { PICARD_SORTVCF        } from '../../modules/local/picard/sortvcf/main'
include { GATK4_GENOMICSDBIMPORT } from '../../modules/nf-core/gatk4/genomicsdbimport/main'

workflow CREATE_GVCF{
    take:
    tup_meta_bam
    fastaF
    fai 
    dict

    main:
    
    //skip_regionwise_gvcf is the parameter for deciding whether or not to create gvcf file in parallel using file
    //provided as the argument for params.regions; this file should end with .list or .intervals. 

    if ( !params.skip_regionwise_gvcf ){    
        Channel
            .fromPath( params.regions )
            .splitText()
            .map{ it -> it.trim() }
            .set{ reg }
        bam_reg = tup_meta_bam.combine(reg)
        }
    else{
        regionF = Channel.fromPath( params.regions )
        bam_reg = tup_meta_bam.combine( regionF )
        }

    //prepare inp channels for the module GATK?_HAPLOTYPECALLER*

    bam_reg_fas = bam_reg.combine( fastaF )
    bam_reg_fas_fai = bam_reg_fas.combine( fai )
    bam_reg_fas_fai_dic = bam_reg_fas_fai.combine( dict )


    if( params.dragstr_model == null ){
            //bam_reg_fas_fai_dic.view()
            bam_reg_fas_fai_dic_drg = bam_reg_fas_fai_dic.combine(["none"])
        }
    else{
            bam_reg_fas_fai_dic_drg = bam_reg_fas_fai_dic.combine( params.dragstr_model )
        }


    if ( params.dbsnp == null ){
            bam_reg_fas_fai_dic_drg_dbs_dbi = bam_reg_fas_fai_dic_drg.combine([["none","none"]])
        }
    else{
            Channel
                .fromFilePairs( params.dbsnp)
                .map{ pfx, dbsnp -> tuple( dbsnp[0], dbsnp[1] ) }
                .set{ snp }

            bam_reg_fas_fai_dic_drg_dbs_dbi = bam_reg_fas_fai_dic_drg.combine( snp )
        }

    
    inp_arg1 = bam_reg_fas_fai_dic_drg_dbs_dbi.map{ meta, bam, idx, rg, fas, fai, dict, drg, dbs, dbi -> tuple( meta, bam, idx, rg, drg == "none" ? [] : drg )}
    inp_arg2 = bam_reg_fas_fai_dic_drg_dbs_dbi.map{ meta, bam, idx, rg, fas, fai, dict, drg, dbs, dbi -> fas }
    inp_arg3 = bam_reg_fas_fai_dic_drg_dbs_dbi.map{ meta, bam, idx, rg, fas, fai, dict, drg, dbs, dbi -> fai }
    inp_arg4 = bam_reg_fas_fai_dic_drg_dbs_dbi.map{ meta, bam, idx, rg, fas, fai, dict, drg, dbs, dbi -> dict }
    inp_arg5 = bam_reg_fas_fai_dic_drg_dbs_dbi.map{ meta, bam, idx, rg, fas, fai, dict, drg, dbs, dbi -> dbs == "none"? [] :dbs }
    inp_arg6 = bam_reg_fas_fai_dic_drg_dbs_dbi.map{ meta, bam, idx, rg, fas, fai, dict, drg, dbs, dbi -> dbi == "none"? [] : dbi}
    


    //
    //if skip_regionwise_gvcf is not false --> do not run gatk haplotypecaller parallel for each region MODULE GATK4_HAPLOTYPECALLER a separate module
    //"GATK_HAPLOTYPECALLER_REG" was created --> because nf-core module for haplotypecaller uses "file" as input channel for the region while here it could also be "val"
    //


    if ( !params.skip_regionwise_gvcf ){    
        GATK_HAPLOTYPECALLER_REG(
            inp_arg1,
            inp_arg2,
            inp_arg3,
            inp_arg4,
            inp_arg5,
            inp_arg6
            )
            
            vcfF = GATK_HAPLOTYPECALLER_REG.out.vcf.groupTuple()
            vcfF_fastaF = vcfF.combine(fastaF)
            vcf_fastaF_dict = vcfF_fastaF.combine(dict)

            inp1_picard_sortvcf = vcf_fastaF_dict.map{ meta, vcfFiles, refFile, dictFile -> tuple(meta, vcfFiles) }
            inp2_picard_sortvcf = vcf_fastaF_dict.map{ meta, vcfFiles, refFile, dictFile -> refFile }
            inp3_picard_sortvcf = vcf_fastaF_dict.map{ meta, vcfFiles, refFile, dictFile -> dictFile }

            //
            // merge and sort chromosome-wise gvcf files emited by the previous module
            //

            PICARD_SORTVCF(
                inp1_picard_sortvcf,
                inp2_picard_sortvcf,
                inp3_picard_sortvcf
            )
            
            mergeGvcf = PICARD_SORTVCF.out.vcf
            mergeTbi  = PICARD_SORTVCF.out.tbi
            mergeGvcfTbi = mergeGvcf.combine( mergeTbi, by: 0 )

        }

    else{
        GATK4_HAPLOTYPECALLER(
            inp_arg1,
            inp_arg2,
            inp_arg3,
            inp_arg4,
            inp_arg5,
            inp_arg6
            )
         mergeGvcf = GATK4_HAPLOTYPECALLER.out.vcf
         mergeTbi  = GATK4_HAPLOTYPECALLER.out.tbi
         mergeGvcfTbi = mergeGvcf.combine( mergeTbi, by: 0 )
        }

    //
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
