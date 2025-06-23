# notes : to generate repbase length list, use extract_length_from_repbase
import scripts.getGCForGlobalGenome

wildcard_constraints:
    status="full|normal|noTE|shuffle",
    context="nCpG|CpG",
    intra="\d+|Inf",
    inter="\d+|Inf"

def getGCgenome(wildcards):
    return scripts.getGCForGlobalGenome.getGCAndGenomeSize("data_" + wildcards.species + "/genome.fa")["gc"]

def getSpecies(wildcards):
    species = ""
    if wildcards["species"] == "hg":
        species = "Homo_sapiens"

    return species

rule all:
    input:
        "results_{species}/equilibrium_gc_full.png",
        "results_{species}/muts_nCpG_full.png",
        "results_{species}/polyA_tes_full.png",
        "results_{species}/gc_tes_full.png"

rule repeatToBed:
    input:
        "data_{species}/Repeatmasker.ucsc"
    output:
        tes = "data_{species}/tes_normal.bed",
        families = "data_{species}/families.tsv"
    shell:
        "python3 scripts/convertAndAssignToFamily.py {input} {output}"

rule repeatFull:
    input:
        "data_{species}/Repeatmasker.ucsc"
    output:
        tes = "data_{species}/tes_full.bed"
    shadow: "shallow"
    shell:
        """
        Rscript scripts/selectComplete.R {input} tmp.ucsc
        python3 scripts/convertAndAssignToFamily.py tmp.ucsc {output} whatever
        """

rule sort:
    input:
        "data_{species}/{anything}.bed"
    output:
        "data_{species}/{anything}_sorted.bed"
    shell:
        "bedtools sort -i {input} > {output}"

rule complement:
    input:
        tes = "data_{species}/tes_normal_sorted.bed",
        genome = "data_{species}/genome.genome"
    output:
        "data_{species}/tes_noTE_fam.bed"
    shell:
        """
        bedtools complement -g {input.genome} -i {input.tes} > {output}
        sed -i "s/$/\tnoTE\t.\t+/g" {output}
        """

rule teToFam:
    input:
        "data_{species}/tes_{status}.bed",
        "data_{species}/families.tsv",
        "data_{species}/Repeatmasker.ucsc"
    output:
        temp("data_{species}/tes_{status}_fam.bed")
    shadow: "shallow"
    wildcard_constraints:
        status="noTE|normal|shuffle"
    shell:
        """
        python3 scripts/tetypeToFam.py {input} tmp.bed
        grep -Pv "\t\t" tmp.bed > {output}
        """

# rule teToFamFull:
#     input:
#         "data_{species}/tes_normal.bed",
#         "data_{species}/families.tsv",
#         "data_{species}/repbase_length_ids.tsv"
#     output:
#         temp("data_{species}/tes_full_fam.bed")
#     shadow: "shallow"
#     params:
#         te_len = "full",
#         species = getSpecies
#     shell:
#         """
#         python3 scripts/tetypeToFam.py {input} tmp.bed {params}
#         grep -Pv "\t\t" tmp.bed > {output}
#         """

rule filterFam:
    input:
        "data_{species}/tes_{status_fam}.bed",
        "data_{species}/families.tsv"
    output:
        temp("data_{species}/tes_{status_fam}_filtered.bed")
    params: "{status_fam}"
    shell:
        "Rscript scripts/filter_fams.R {input} {output} {params}"

rule createGenomeFile:
    input:
        "data_{species}/genome.fa"
    output:
        "data_{species}/genome.genome"
    shell:
        "python3 scripts/createGenomeFile.py {input} {output}"

rule filter_intern_nieb:
    input:
        niebs = "data_{species}/niebs.bed",
        tes = "data_{species}/tes_normal_sorted.bed",
        aoe = "data_{species}/niebs.aoe"
    output:
        "data_{species}/no_intern_niebs.aoe"
    shadow: "shallow"
    shell:
        """
        bedtools intersect -v -a {input.niebs} -b {input.tes} -f 1 > tmp.niebs
        bedtools intersect -u -a {input.aoe} -b tmp.niebs > {output}
        """

rule select_niebs:
    input:
        "data_{species}/no_intern_niebs.aoe"
    output:
        "data_{species}/{intra}_{inter}_niebs.aoe"
    params:
        intra = "{intra}",
        inter = "{inter}"
    shell:
        "Rscript scripts/selectNIEBs.R {input} {params} {output}"

rule getSequences:
    input:
        bed = "data_{species}/tes_full_fam.bed",
        genome = "data_{species}/genome.fa"
    output:
        "data_{species}/seq.fa"
    shell:
        "bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output} -s"

rule getPolyA:
    input:
        "data_{species}/seq.fa"
    output:
        "data_{species}/polyA.bed"
    threads: workflow.cores
    shell:
        "python3 scripts/annotate2.py {input} {output} {threads}"

rule renamePolByTE:
    input:
        "data_{species}/polyA.bed",
        "data_{species}/tes_full_fam.bed"
    output:
        temp("data_{species}/polyA.renamed.bed")
    shell:
        "Rscript scripts/renamePolyA.R {input} {output}"

rule getPosPolyA:
    input:
        aoe = "data_{species}/niebs.aoe",
        polyA = "data_{species}/polyA.renamed.bed"
    output:
        "data_{species}/polyA.tsv"
    shell:
        """
        bin/countFeatures.bin +a {input.aoe} +b {input.polyA} +o {output} +m +i 1000000 +p start +d
        """

rule getPosTESpecific:
    input:
        aoe = "data_{species}/{intra}_{inter}_niebs.aoe",
        tes = "data_{species}/tes_{status_fam}.bed"
    output:
        "data_{species}/count_{status_fam}_{intra}_{inter}.tsv"
    shell:
        """
        bin/countFeatures.bin +a {input.aoe} +b {input.tes} +o {output} +p stop +s +d +k hit
        """ # final output is pos val type

rule getNCpG:
    input:
        "data_{species}/anc_genome.fa"
    output:
        CpG = "data_{species}/CpG.bed",
        nCpG = "data_{species}/nCpG.bed"
    shadow: "shallow"
    shell:
        "bin/getPattern.bin +f {input} +1 {output.CpG} +2 {output.nCpG} +p CG"

rule getRefCounts:
    input:
        "data_{species}/{intra}_{inter}_niebs.aoe"
    output:
        "data_{species}/ref_counts_{intra}_{inter}.tsv"
    shell:
        "bin/countFeatures.bin +a {input} +b {input} +o {output}"

rule preIntersect:
    input:
        tes = "data_{species}/{category}_{status}_fam_filtered.bed",
        mask = "data_{species}/{context}.bed"
    output:
       "data_{species}/preinter_{context}_{category}_{status}_fam_filtered.bed"
    wildcard_constraints:
        category = "first_pos|tes"
    shell:
        "bin/intersectKeepingNames.bin +a {input.tes} +b {input.mask} +o {output}"

rule count_bases:
    input:
        aoe = "data_{species}/{intra}_{inter}_niebs.aoe",
        genome = "data_{species}/anc_genome.fa",
        mask = "data_{species}/preinter_{context}_{category}_{status}_fam_filtered.bed"
    output:
        "data_{species}/bases_{context}_{category}_{status}_{intra}_{inter}.tsv"
    shell:
        "bin/countBases.bin +f {input.genome} +a {input.aoe} +b {input.mask} +o {output} +i +m hit"

rule count_bases_full:
    input:
        aoe = "data_{species}/{intra}_{inter}_niebs.aoe",
        genome = "data_{species}/genome.fa",
        mask = "data_{species}/{category}_{status}_fam_filtered.bed"
    output:
        "data_{species}/bases_{category}_{status}_{intra}_{inter}.tsv"
    wildcard_constraints:
         category="first_pos|tes"
    shell:
        "bin/countBases.bin +f {input.genome} +a {input.aoe} +b {input.mask} +o {output} +i +m hit"

rule make_AOE_from_te:
    input:
        "data_{species}/tes_full.bed"
    output:
        "data_{species}/tes_full.aoe"
    shell:
        "Rscript scripts/extract_all_starts.R {input} {output}"

rule count_niebs_on_te:
    input:
        tes = "data_{species}/tes_full.aoe",
        niebs = "data_{species}/niebs.bed"
    output:
        "data_{species}/niebs_on_tes.tsv"
    shell:
        "./bin/countFeatures.bin +a {input.tes} +b {input.niebs} +o {output} +d +k source +s +p stop +f"

rule count_ref_te:
    input:
        "data_{species}/tes_full.aoe"
    output:
        "data_{species}/ref_count_tes.tsv"
    shell:
        "./bin/countFeatures.bin +a {input} +b {input} +o {output} +d +k source +s +p whole +f"

rule plot_nieb_on_te:
    

rule getProfileGC:
    input:
        "data_{species}/bases_{category}_{status}_{intra}_{inter}.tsv"
    output:
        "data_{species}/gc_{category}_{status}_{intra}_{inter}.tsv"
    shell:
        "Rscript scripts/getGCFromCountPerID.R {input} {output}"

rule count_mutations:
    input:
        aoe = "data_{species}/{intra}_{inter}_niebs.aoe",
        muts = "data_{species}/muts.vcf",
        mask = "data_{species}/preinter_{context}_{category}_{status}_fam_filtered.bed"
    output:
        "data_{species}/muts_{context}_{category}_{status}_{intra}_{inter}.tsv"
    shell:
        "bin/countMuts.bin +v {input.muts} +a {input.aoe} +b {input.mask} +o {output} +j +m hit"

rule getRate:
    input:
        "data_{species}/bases_{context}_{category}_{status}_{intra}_{inter}.tsv",
        "data_{species}/muts_{context}_{category}_{status}_{intra}_{inter}.tsv"
    output:
        "data_{species}/rates_{context}_{category}_{status}_{intra}_{inter}.tsv"
    shell:
        """
        Rscript scripts/rateSW.R {input} {output}
        """

rule gillespiage:
    input:
        CpG_rates = "data_{species}/rates_CpG_{category}_{status}_{intra}_{inter}.tsv",
        nCpG_rates = "data_{species}/rates_nCpG_{category}_{status}_{intra}_{inter}.tsv"
    output:
        "data_{species}/equilibrium_{category}_{status}_{intra}_{inter}.tsv"
    shadow: "shallow"
    threads: workflow.cores
    shell:
        """
        Rscript scripts/conformMutsV2.R {input.CpG_rates} muts_CpG.tsv
        Rscript scripts/conformMutsV2.R {input.nCpG_rates} muts_nCpG.tsv
        Rscript scripts/simplGillespie.R muts_nCpG.tsv muts_CpG.tsv {output} {threads}
        """

rule plotCov:
    input:
        "data_{species}/count_{status_fam}_{intra}_{inter}.tsv",
        "data_{species}/ref_counts_{intra}_{inter}.tsv",
        "data_{species}/families.tsv"
    output:
        directory("results_{species}/cov_{status_fam}_{intra}_{inter}"),
        "results_{species}/cov_{status_fam}_savestate_{intra}_{inter}.csv"
    threads: workflow.cores
    shell:
        "Rscript scripts/plotPos.R {input} {output} {threads}"

rule plotPolyA:
    input:
        "data_{species}/polyA.tsv",
        "data_{species}/ref_counts.tsv"
    output:
        "results_{species}/polyA.png",
        "results_{species}/polyA_savestate.csv"
    shell:
        "Rscript scripts/plotPolyA.R {input} {output}"

rule plotGC:
    input:
        "data_{species}/gc_{category}_{status}_{intra}_{inter}.tsv"
    output:
        "results_{species}/gc_{category}_{status}_{intra}_{inter}.png",
        "results_{species}/gc_{category}_{status}_savestate_{intra}_{inter}.csv"
    params:
        getGCgenome
    shell:
        "Rscript scripts/plotGC.R {input} {params} {output}"
        
rule plotMuts:
    input:
        "data_{species}/rates_{context}_{status}_{intra}_{inter}.tsv"
    output:
        "results_{species}/muts_{context}_{status}_{intra}_{inter}.png",
        "results_{species}/muts_{context}_{status}_{intra}_{inter}_savestate.csv"
    shell:
        "Rscript scripts/plotMuts.R {input} {output}"

rule plotEquilibrium:
    input:
        "data_{species}/equilibrium_{category}_{status}_{intra}_{inter}.tsv"
    output:
        "results_{species}/equilibrium_gc_{category}_{status}_{intra}_{inter}.png",
        "results_{species}/equilibrium_{category}_{status}_{intra}_{inter}_savestate.csv"
    shell:
        "Rscript scripts/plotEquilibrium.R {input} {output}"
