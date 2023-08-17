import os
import math

configfile: "config.yml"
#singularity: "topmed-methy.sif"

mem_step_size = 6400

rule invnorm_betas:
    input: 
        config["rds_file"]
    output:
        "phenotypes/betas.sorted_samples.invnorm.tsv.gz"
    threads: 1
    resources:
        mem_mb = 128000
    shell:
        """
        set -eu

        process-methy-rds.R {input} {output}
        """

rule invnorm_betas_with_positions:
    input:
        betas = rules.invnorm_betas.output,
        position_map = config["cpg_pos_map"]
    output:
        "phenotypes/betas.sorted_samples.invnorm.ids_fixed.samples_filtered.with_positions.tsv.gz"
    threads: 1
    resources:
        mem_mb = 64000
    shell:
        """
        set -euo pipefail

        merge-cpg-positions.sh {input.betas} {input.position_map} \
          | subset-methy.R \
          | gzip > {output}
        """

rule bed_file:
    input:
        rules.invnorm_betas_with_positions.output
    output:
        "bed_files/betas.sorted_samples.invnorm.ids_fixed.samples_filtered.{chrom}.bed.gz"
    threads: 1
    resources:
        mem_mb = 6400
    shell:
        """
        set -euo pipefail
        (echo '#chr' start end probe_id $(zcat {input} | head -n1 | cut -f5-) | tr " " "\t"; zcat {input} | awk -F'\t' 'BEGIN {{OFS="\t"}} {{if ($1=="{wildcards.chrom}") {{$3=$2+1; print}}}}' | sort -k2,2 --numeric-sort) | bgzip > {output}
        tabix {output}
        """


rule genotype_subset:
    input:
        vcf = config["genotype_expression"],
        id_convert_map = config["sample_id_convert_map"]
    output:
        "genotypes/{chrom}.pass_variants.bed"
    threads: 2
    resources:
        mem_mb = lambda wc, attempt:  mem_step_size + mem_step_size * attempt
    shell:
        """
        set -uo pipefail
        tmp_dir=`mktemp -d`

        plink_prefix_path=${{tmp_dir}}/$(basename {output} .bed)
        bcftools view -S <(cut -f1 {input.id_convert_map} | grep -v "^NA$" | grep -v "NWD359245\|NWD439884\|NWD550834\|NWD582434\|NWD609471\|NWD741494") {input.vcf} --include "FILTER='PASS'" -Ou \
          | bcftools view --include 'AC>0&&AC!=AN' -Ou \
          | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou \
          | bcftools reheader --samples {input.id_convert_map} /dev/stdin \
          | bcftools view -S <(cut -f2 {input.id_convert_map} | sort -k2,2) --force-samples -Ob -o ${{plink_prefix_path}}.bcf
        bcftools index ${{plink_prefix_path}}.bcf

        #bcftools view --include '(AC/AN)>0.01' ${{plink_prefix_path}}.bcf -Ob -o ${{plink_prefix_path}}.min_af_0.01.bcf
        #bcftools index ${{plink_prefix_path}}.min_af_0.01.bcf

        plink2 --make-bed --output-chr chrM --bcf ${{plink_prefix_path}}.bcf --threads 4 --out $plink_prefix_path
        mv ${{plink_prefix_path}}* genotypes/
        """


rule all_genotype_subsets:
    input:
        ["genotypes/chr" + c + ".pass_variants.bed" for c in [str(i) for i in range(1, 23)] + ["X"]]


rule null_genotype_subset:
    input:
        rules.genotype_subset.output
    output:
        "genotypes_null/{chrom}.pass_variants.null.bed"
    threads: 2
    resources:
        mem_mb = lambda wc, attempt:  mem_step_size + mem_step_size * attempt
    shell:
        """
        set -uo pipefail
        plink_prefix_path=$(dirname {output})/$(basename {output} .bed)
        input_bcf=$(dirname {input})/$(basename {input} .bed).bcf
        bcftools reheader --samples <(bcftools query -l $input_bcf | sort --random-sort) $input_bcf -o ${{plink_prefix_path}}.bcf
        bcftools index --force ${{plink_prefix_path}}.bcf
        plink2 --make-bed --output-chr chrM --bcf ${{plink_prefix_path}}.bcf --threads 4 --out $plink_prefix_path
        """

rule pcs:
    input:
        rules.invnorm_betas_with_positions.output
    output:
        "pcs/betas.sorted_samples.invnorm.ids_fixed.samples_filtered.loadings.tsv",
        "pcs/betas.sorted_samples.invnorm.ids_fixed.samples_filtered.pcs_variance_explained.txt"
    params:
        prefix = "pcs/betas.sorted_samples.invnorm.ids_fixed.samples_filtered"
    threads: 2
    resources:
        mem_mb = 64000
    shell:
        """
        set -eu
       
        generate-methy-pcs.R {input} {params.prefix}
        tmp={params.prefix}.loadings.tsv.tmp
        (printf "sample_id\t$(head -n1 $tmp | cut -f2-)\n"; tail -n+2 $tmp) > {output[0]}
        rm $tmp
        """


rule covariates:
    input:
        sample_sheet = "MISSING.csv",
        geno_pcs = "MISSING.tsv",
        methy_pcs = rules.pcs.output[0]
    output:
        "covariates/sex_gpcs_mpcs.tsv", "covariates/sex_gpcs_mpcs.transposed.tsv"
    threads: 1
    resources:
        mem_mb = mem_step_size
    shell:
        """
        set -euo pipefail

        out1=$(echo {output[0]} | rev | cut -f3- -d'.' | rev).tmp
        head -n1 {input.sample_sheet} | cut -f1,3 -d',' | tr "," "\t" | join --header -t$'\t' <(head -n1 {input.geno_pcs}) /dev/stdin > $out1
        tail -n+2 {input.sample_sheet} | cut -f1,3 -d',' | tr "," "\t" | awk -F'\t' 'BEGIN{{OFS="\t"}} {{if ($2=="F") {{ $2=1 }} else {{$2=0}}}} {{print}}' | sort -k1,1 | join -t$'\t' <(tail -n+2 {input.geno_pcs} | sort -k1,1) /dev/stdin >> $out1

        join --header -t$'\t' <(head -n1 $out1) <(head -n1 {input.methy_pcs}) > {output[0]}
        join -t$'\t' <(tail -n+2 $out1) <(tail -n+2 {input.methy_pcs}) >> {output[0]}

        cat {output[0]} | transpose_tsv.R > {output[1]}
        """


rule methy_pcs_plot:
    input:
        rules.covariates.output[0]
    output:
        "plots/methy_pcs.png"
    shell:
        """
        Rscript plot-mpcs.R {input} {output}
        """

rule methy_pcs_variance_explained_plot:
    input:
        rules.pcs.output[1]
    output:
        "plots/methy_pcs_variance_explained.png"
    shell:
        """
        set -euo pipefail
        paste <(seq 1 50) <(head -n 50 {input}) | gnuplot --persist -e 'set term png size 1080,720; set xrange [0:51]; set xlabel "PC"; set ylabel "% Variance Explained"; set style fill solid noborder; set boxwidth 0.8; plot "-" u 1:2 w boxes notitle' > {output}
        """

rule null_covariates:
    input:
        vcf = rules.genotype_subset.output[0].format(chrom="chr20"),
        vcf_null = rules.null_genotype_subset.output[0].format(chrom="chr20"),
        cov = "covariates/joint_invnorm_loadings.ids_fixed.sex_study_pcs_1_100.recode.with_geno_pcs.tsv"
    output:
        "covariates_null/joint_invnorm_loadings.ids_fixed.sex_study_pcs_1_100.recode.with_geno_pcs.null2.tsv",
        "covariates_null/joint_invnorm_loadings.ids_fixed.sex_study_pcs_1_100.recode.with_geno_pcs.null2.transposed.tsv"
    threads: 2
    resources:
        mem_mb = lambda wc, attempt:  mem_step_size + mem_step_size * attempt
    shell:
        """
        set -euo pipefail

        plink_prefix=$(dirname {input.vcf})/$(basename {input.vcf} .bed)
        diff <(cut -f2 ${{plink_prefix}}.fam) <(tail -n+2 {input.cov} | cut -f1)

        plink_prefix_null=$(dirname {input.vcf_null})/$(basename {input.vcf_null} .bed)
        tmp_dir=covariates_null
        paste <(cut -f2 ${{plink_prefix_null}}.fam) <(tail -n+2 {input.cov} | cut -f 2-11) | sort -k1,1 > $tmp_dir/null_pcs.tsv
        (head -n1 {input.cov}; paste $tmp_dir/null_pcs.tsv <(tail -n+2 {input.cov} | cut -f12-)) > {output[0]}
        cat {output[0]} | ./transpose_tsv.R > {output[1]}        

        rm $tmp_dir/null_pcs.tsv 
        
        """
