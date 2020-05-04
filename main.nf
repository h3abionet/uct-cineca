#!/usr/bin/env nextflow
nextflow.preview.dsl=2


/* 
 * Get sample in vcf
 */
process get_samples {

    input:
    tuple val(dataset), val(size), file(csv), chrm, file(vcf)

    output:
    tuple val(dataset), file(csv_sample), chrm, file(vcf)
    tuple val(size), file(vcf_sample)

    script:
    csv_sample = "${csv.baseName}.samples"
    vcf_sample = "${vcf.baseName}.samples"
    """
    awk -F';' 'NR>1 {print \$3}' ${csv} > ${csv_sample}
    bcftools query -l ${vcf} > ${vcf_sample}
    """
}


/* 
 * Select sample
 */
process random_select_sample {

    input:
    tuple sample_size, file(vcf_sample_file)

    output:
    file(vcf_sample_reduced)

    script:
    base = file(vcf_sample_file.baseName).baseName
    vcf_sample_reduced = "${vcf_sample_file}_reduced.samples"
    """
    sort -R ${vcf_sample_file} | tail -n ${sample_size} > ${vcf_sample_reduced}
    """
}

process extract_sample {

    input:
    tuple dataset, file(pheno_sample_file), chrm, file(vcf_file), file(vcf_sample_reduced)

    output:
    tuple dataset, file(pheno_sample_file), chrm, file(vcf_reduced)

    script:
    base = file(vcf_file.baseName).baseName
    vcf_reduced = "${base}_reduced.vcf.gz"
    """
    bcftools view \
        --samples-file ${vcf_sample_reduced} \
        --min-af 0.2 \
        -m2 -M2 -v snps \
        ${vcf_file} | \
    bcftools +fill-tags | \
         bgzip -c > ${vcf_reduced}
    """
}


process extract_variant {

    input:
    tuple dataset, file(pheno_sample_file), chrm, file(vcf_reduced), variant_size

    output:
    tuple dataset, file(pheno_sample_file), chrm, file(vcf_reduced_resized)

    script:
    base = file(vcf_reduced.baseName).baseName
    bed_vcf_reduced = "${base}_reduced_${variant_size}.bed"
    vcf_reduced_resized = "${base}_reduced_${variant_size}.vcf.gz"
    """
    tabix ${vcf_reduced}
    bcftools view \
         --min-af 0.1 \
         -m2 -M2 -v snps \
         ${vcf_reduced} |
    bcftools query \
       -f '%CHROM\\t%POS\\n' |
    sort -R | tail -n ${variant_size} | \
    awk '{print \$1"\\t"\$2"\\t"\$2+1}' > ${bed_vcf_reduced}
    bcftools view \
        --regions-file ${bed_vcf_reduced} \
        ${vcf_reduced} | \
    bcftools sort \
        -Oz -o ${vcf_reduced_resized}
    """
}


process replace_sample {

    input:
    tuple dataset, file(pheno_sample_file), chrm, file(vcf_reduced_resized)

    output:
    tuple dataset, file(vcf_reduced)

    script:
    vcf_reduced = "${dataset}_chr${chrm}.vcf.gz"
    """
    zcat ${vcf_reduced_resized} | \
    bcftools reheader \
        --samples ${pheno_sample_file} | \
    bgzip -c > ${vcf_reduced}
    """
}

/* 
 * Concatenate chromosome vcfs
 */
process concat_vcfs {

    input:
    tuple dataset, vcfs

    output:
    tuple dataset, file(vcf_out)

    script:
    vcf_out = "${dataset}.vcf.gz"
    """
    bcftools concat \
        ${vcfs.join(' ')} | \
    bcftools sort |
    bcftools +fill-tags | \
         bgzip -c > ${vcf_out}
    """
}




workflow {
    // check if ref files exist
    datas = []
    params.chrms.each { chrm ->
        vcf = sprintf(params.vcf, chrm)
        if(!file(vcf).exists()) exit 1, "File ${vcf} not found. Please check your config file."
        params.csvs.each{ dataset, csv ->
            if(!file(csv).exists()) exit 1, "File ${csv} not found. Please check your config file."
            datas << [dataset, file(csv).countLines(), file(csv), chrm, file(vcf)]
        }
    }
    
    // Step 1
    csv_data = Channel.from(datas)   
    get_samples(csv_data)
    // Step 2 Random sampling
    random_select_sample(get_samples.out[1].first())
    // Step 3
    extract_sample(get_samples.out[0].combine(random_select_sample.out))
    // Step 4
    extract_variant(extract_sample.out.combine([params.variant_size]))
    // Step 5
    replace_sample(extract_variant.out)
    // Step 6
    concat_vcfs(replace_sample.out.groupTuple(by:0)).view()
}