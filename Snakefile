# notes : to generate repbase length list, use extract_length_from_repbase
import math
import scripts.getGCForGlobalGenome

wildcard_constraints:
    aggstatus="full|full_fam",
    context="nCpG|CpG",
    intra="\d+|Inf",
    inter="\d+|Inf",
    types = "tes|niebs"

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

rule repeatFamComplete:
    input:
        "data_{species}/Repeatmasker.ucsc"
    output:
        tes = "data_{species}/tes_full_fam.bed",
        sizes = "data_{species}/selected_sizes.csv"
    shadow: "shallow"
    shell:
        """
        Rscript scripts/agg_by_fam.R {input} tmp.ucsc {output.sizes}
        python3 scripts/convertAndAssignToFamily.py tmp.ucsc {output.tes} whatever
        """

rule summarizeTEs:
    input:
        "data_{species}/tes_{aggstatus}.bed"
    output:
        'data_{species}/summary_{aggstatus}.csv'
    shell:
        "Rscript scripts/tes_per_fam.R {input} {output}"

rule sort:
    input:
        "data_{species}/{anything}.bed"
    output:
        "data_{species}/{anything}_sorted.bed"
    shell:
        "bedtools sort -i {input} > {output}"

# rule complement:
#     input:
#         tes = "data_{species}/tes_normal_sorted.bed",
#         genome = "data_{species}/genome.genome"
#     output:
#         "data_{species}/tes_noTE_fam.bed"
#     shell:
#         """
#         bedtools complement -g {input.genome} -i {input.tes} > {output}
#         sed -i "s/$/\tnoTE\t.\t+/g" {output}
#         """

# rule filterFam:
#     input:
#         "data_{species}/tes_{status_fam}.bed",
#         "data_{species}/families.tsv"
#     output:
#         temp("data_{species}/tes_{status_fam}_filtered.bed")
#     params: "{status_fam}"
#     shell:
#         "Rscript scripts/filter_fams.R {input} {output} {params}"

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
        "data_{species}/niebs.aoe"
        # "data_{species}/no_intern_niebs.aoe"
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

# rule getPolyA:
#     input:
#         "data_{species}/seq.fa"
#     output:
#         "data_{species}/polyA.bed"
#     threads: workflow.cores / 3
#     shell:
#         "python3 scripts/annotate2.py {input} {output} {threads}"

# rule renamePolByTE:
#     input:
#         "data_{species}/polyA.bed",
#         "data_{species}/tes_full_fam.bed"
#     output:
#         temp("data_{species}/polyA.renamed.bed")
#     shell:
#         "Rscript scripts/renamePolyA.R {input} {output}"

# rule getPosPolyA:
#     input:
#         aoe = "data_{species}/niebs.aoe",
#         polyA = "data_{species}/polyA.renamed.bed"
#     output:
#         "data_{species}/polyA.tsv"
#     shell:
#         """
#         bin/countFeatures.bin +a {input.aoe} +b {input.polyA} +o {output} +m +i 1000000 +p start +d
#         """

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

rule getPosCopySpecific:
    input:
        aoe = "data_{species}/{intra}_{inter}_niebs.aoe",
        tes = "data_{species}/tes_full_fam.bed"
    output:
        temp("data_{species}/pos_copy_specific_{intra}_{inter}.tsv")
    shell:
        """
        Rscript scripts/idToChrStartEnd.R {input.tes} tmp.tes.bed
        bin/countFeatures.bin +a {input.aoe} +b tmp.tes.bed +o {output} +d +k hit +s +p stop
        """

rule selectNIEBs:
    input:
        tes_copy = "data_{species}/pos_copy_specific_{intra}_{inter}.tsv",
        niebs = "data_{species}/niebs.bed",
        aoes = "data_{species}/{intra}_{inter}_niebs.aoe"
    output:
        "data_{species}/first_pos_{intra}_{inter}.aoe"
    shadow: "shallow"
    shell:
        """
        echo start
        grep -P '\\-[1-50]\\t1' {input.tes_copy} | grep -P "\\-\\t\\+|\\+\\t\\-" | cut -f1 | sed -E 's/_/\\t/g' > tmp.copy.bed
        echo intersect
        bedtools intersect -b tmp.copy.bed -a {input.niebs} -u > tmp.niebs.bed
        echo last intersect
        bedtools intersect -b tmp.niebs.bed -a {input.aoes} -wa > {output}
        """

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

rule plotPos:
    input:
        "data_{species}/count_{status_fam}_{intra}_{inter}.tsv",
        "data_{species}/ref_counts_{intra}_{inter}.tsv",
    output:
        directory("results_{species}/cov_{status_fam}_{intra}_{inter}"),
        "results_{species}/cov_{status_fam}_savestate_{intra}_{inter}.csv",
        "results_{species}/cov_{status_fam}_stats_{intra}_{inter}.csv"
    threads: math.floor(workflow.cores / 3)
    shell:
        "Rscript scripts/plotPos.R {input} {output} {threads}"

rule preIntersect:
    input:
        tes = "data_{species}/tes_{aggstatus}.bed", # correct this later to use aggregated values at fam level?
        mask = "data_{species}/{context}.bed"
    output:
       "data_{species}/preinter_{context}_{aggstatus}.bed"
    shell:
        "bin/intersectKeepingNames.bin +a {input.tes} +b {input.mask} +o {output}"

rule count_bases:
    input:
        aoe = "data_{species}/first_pos_{intra}_{inter}.aoe",
        genome = "data_{species}/anc_genome.fa",
        mask = "data_{species}/preinter_{context}_{aggstatus}.bed"
    output:
        "data_{species}/bases_{context}_{aggstatus}_{intra}_{inter}.tsv"
    threads: min(workflow.cores, 5)
    shell:
        "bin/countBases.bin +f {input.genome} +a {input.aoe} +b {input.mask} +o {output} +i +m hit +t {threads}"

rule count_bases_tes:
    input:
        aoe = "data_{species}/tes_{aggstatus}.aoe",
        genome = "data_{species}/anc_genome.fa",
        mask = "data_{species}/{context}.bed"
    output:
        "data_{species}/bases_tes_{context}_{aggstatus}_{intra}_{inter}.tsv"
    threads: min(workflow.cores, 5)
    shell:
        "bin/countBases.bin +f {input.genome} +a {input.aoe} +b {input.mask} +o {output} +i +m source +t {threads}"

rule count_gc_niebs_full:
    input:
        aoe = "data_{species}/first_pos_{intra}_{inter}.aoe",
        genome = "data_{species}/genome.fa",
        mask = "data_{species}/tes_{aggstatus}.bed"
    output:
        "data_{species}/bases_niebs_{aggstatus}_{intra}_{inter}.tsv"
    threads: min(workflow.cores, 5)
    shell:
        "bin/countBases.bin +f {input.genome} +a {input.aoe} +b {input.mask} +o {output} +i +m hit +g +t {threads}"

rule count_gc_tes_full:
    input:
        # mask = "data_{species}/niebs.bed",
        genome = "data_{species}/genome.fa",
        aoe = "data_{species}/tes_{aggstatus}.aoe"
    output:
        "data_{species}/bases_tes_{aggstatus}_{intra}_{inter}.tsv"
    threads: min(workflow.cores, 5)
    shell:
        "bin/countBases.bin +f {input.genome} +a {input.aoe} +o {output} +i +m source +g +t {threads}"

rule select_spaced_tes:
    input:
        "data_{species}/tes_{aggstatus}.bed"
    output:
        # temp("data_{species}/tmp.ids.{aggstatus}.txt")
        temp("data_{species}/tes_{aggstatus}_nooverlap.bed")
    shadow: "shallow"
    shell:
        """ # python3 scripts/extend_tes.py {input} tmp.bed
        bedtools intersect -a {input} -b {input} -c > tmp.bed
        grep 1$ tmp.bed | cut -f1-6 | bedtools sort -i - > {output}
        """

rule make_AOE_from_te:
    input:
        "data_{species}/tes_{aggstatus}_nooverlap_sorted.bed"
    output:
        "data_{species}/tes_{aggstatus}.aoe"
    shell:
        "Rscript scripts/extract_all_starts.R {input} {output}"

rule count_niebs_on_te:
    input:
        tes = "data_{species}/tes_{aggstatus}.aoe",
        niebs = "data_{species}/niebs.bed"
    output:
        "data_{species}/niebs_on_tes_{aggstatus}.tsv"
    shell:
        "./bin/countFeatures.bin +a {input.tes} +b {input.niebs} +o {output} +d +k source +p whole"

rule count_nucs_on_te:
    input:
        read = "data_{species}/reads_nucs.bed",
        aoe = "data_{species}/tes_{aggstatus}.aoe"
    output:
        "data_{species}/reads_counts_{aggstatus}.tsv"
    shell:
        "bin/countFeatures.bin +a {input.aoe} +b {input.read} +o {output} +d +k both +p mid"

rule count_ref_te:
    input:
        "data_{species}/tes_{aggstatus}.aoe"
    output:
        "data_{species}/ref_count_tes_{aggstatus}.tsv"
    shell:
        "./bin/countFeatures.bin +a {input} +b {input} +o {output} +d +k source +p whole"

rule plot_nieb_on_te:
    input:
        "data_{species}/niebs_on_tes_{aggstatus}.tsv",
        "data_{species}/ref_count_tes_{aggstatus}.tsv"
    output:
        directory("results_{species}/plot_niebs_on_tes_{aggstatus}"),
        "results_{species}/savestate_nieb_on_te_{aggstatus}.csv"
    shell:
        "Rscript scripts/plot_te_counts_on_niebs.R {input} {output} {threads}"

rule plot_nuc_on_te:
    input:
        "data_{species}/reads_counts_{aggstatus}.tsv",
        "data_{species}/ref_count_tes_{aggstatus}.tsv"
    output:
        directory("results_{species}/plot_nucs_on_tes_{aggstatus}"),
        "results_{species}/savestate_nuc_on_te_{aggstatus}.csv"
    shell:
        "Rscript scripts/plot_nucs_counts_on_niebs.R {input} {output}"

rule getProfileGC:
    input:
        "data_{species}/bases_{types}_{aggstatus}_{intra}_{inter}.tsv"
    output:
        "data_{species}/gc_{types}_{aggstatus}_{intra}_{inter}.tsv"
    shell:
        "Rscript scripts/getGCFromCountPerID.R {input} {output}"

rule count_mutations:
    input:
        aoe = "data_{species}/first_pos_{intra}_{inter}.aoe",
        muts = "data_{species}/muts.vcf",
        mask = "data_{species}/preinter_{context}_{aggstatus}.bed"
    output:
        "data_{species}/muts_{context}_{aggstatus}_{intra}_{inter}.tsv"
    shell:
        "bin/countMuts.bin +v {input.muts} +a {input.aoe} +b {input.mask} +o {output} +j +m hit"

rule count_mutations_tes:
    input:
        aoe = "data_{species}/tes_{aggstatus}.aoe",
        muts = "data_{species}/muts.vcf",
        mask = "data_{species}/{context}.bed"
    output:
        "data_{species}/muts_tes_{context}_{aggstatus}_{intra}_{inter}.tsv"
    shell:
        "bin/countMuts.bin +v {input.muts} +a {input.aoe} +b {input.mask} +o {output} +j +m source"

rule getRate:
    input:
        "data_{species}/bases_{spe_context}_{aggstatus}_{intra}_{inter}.tsv",
        "data_{species}/muts_{spe_context}_{aggstatus}_{intra}_{inter}.tsv"
    output:
        "data_{species}/rates_{spe_context}_{aggstatus}_{intra}_{inter}.tsv"
    shell:
        """
        Rscript scripts/rateSW.R {input} {output}
        """

rule gillespiage:
    input:
        CpG_rates = "data_{species}/rates_CpG_{aggstatus}_{intra}_{inter}.tsv",
        nCpG_rates = "data_{species}/rates_nCpG_{aggstatus}_{intra}_{inter}.tsv"
    output:
        "data_{species}/equilibrium_{aggstatus}_{intra}_{inter}.tsv"
    shadow: "shallow"
    threads: workflow.cores
    shell:
        """
        Rscript scripts/conformMutsV2.R {input.CpG_rates} muts_CpG.tsv
        Rscript scripts/conformMutsV2.R {input.nCpG_rates} muts_nCpG.tsv
        Rscript scripts/simplGillespie.R muts_nCpG.tsv muts_CpG.tsv {output} {threads}
        """

rule simplConform:
    input:
        CpG_rates = "data_{species}/rates_tes_CpG_{aggstatus}_{intra}_{inter}.tsv",
        nCpG_rates = "data_{species}/rates_tes_nCpG_{aggstatus}_{intra}_{inter}.tsv"
    output:
        CpG = "data_{species}/conform_tes_CpG_{aggstatus}_{intra}_{inter}.tsv",
        nCpG = "data_{species}/conform_tes_nCpG_{aggstatus}_{intra}_{inter}.tsv"
    shadow: "shallow"
    shell:
        """
        Rscript scripts/conformMutsV2.R {input.CpG_rates} {output.CpG}
        Rscript scripts/conformMutsV2.R {input.nCpG_rates} {output.nCpG}
        """

# rule plotPolyA:
#     input:
#         "data_{species}/polyA.tsv",
#         "data_{species}/ref_counts.tsv"
#     output:
#         "results_{species}/polyA.png",
#         "results_{species}/polyA_savestate.csv"
#     shell:
#         "Rscript scripts/plotPolyA.R {input} {output}"

rule plotGC:
    input:
        "data_{species}/gc_{types}_{aggstatus}_{intra}_{inter}.tsv"
    output:
        directory("results_{species}/gc_{types}_{aggstatus}_{intra}_{inter}"),
        "results_{species}/gc_{types}_{aggstatus}_savestate_{intra}_{inter}.csv"
    params:
        getGCgenome
    shell:
        "Rscript scripts/plotGC.R {input} {params} {output}"
        
rule plotMuts:
    input:
        "data_{species}/rates_{context}_{aggstatus}_{intra}_{inter}.tsv"
    output:
        directory("results_{species}/muts_{context}_{aggstatus}_{intra}_{inter}"),
        "results_{species}/muts_{context}_{aggstatus}_{intra}_{inter}_savestate.csv"
    shell:
        "Rscript scripts/plotMuts.R {input} {output}"

rule plotEquilibrium:
    input:
        "data_{species}/equilibrium_{aggstatus}_{intra}_{inter}.tsv"
    output:
        directory("results_{species}/equilibrium_gc_{aggstatus}_{intra}_{inter}"),
        "results_{species}/equilibrium_{aggstatus}_{intra}_{inter}_savestate.csv"
    shell:
        "Rscript scripts/plotEquilibrium.R {input} {output}"
