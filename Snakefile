configfile: "config.yaml"
wildcard_constraints:
    te_type="[\w-]+",
    status="[\w-]+"

rule repeatToBed:
    input:
        config["DFAM"]
    output:
        tes = "data/tes_normal.bed"
    shell:
        "python3 scripts/convertDfamAnnotToBed.py {input} {output}"

rule getFamilies:
    input:
        config["DFAM_CLASS"],
        config["DFAM_LENGTH"]
    output:
        "data/dfam_families.tsv"
    shell:
        "python3 scripts/createFamilyFile.py {input} {output}"

rule teToFam:
    input:
        "data/tes_normal.bed",
        config["DFAM"],
        "data/dfam_families.tsv"
    output:
        temp("data/tes_{status}_fam.bed")
    params:
        "{status}" # either full or partial
    shell:
        "python3 scripts/tetypeToFam.py {input} {output} {params}"

rule createGenomeFile:
    input:
        config["GENOME"]
    output:
        "data/genome.genome"
    shell:
        "python3 scripts/createGenomeFile.py {input} {output}"

# rule select_tes:
#     input:
#         "data/tes_{status}_fam.bed"
#     output:
#         "data/tes_{status}_fam_{te_type}.bed"
#     params:
#         "{te_type}"
#     shell:
#         "Rscript scripts/selectTEsByName.R {input} <(echo {params}) {output}"

rule getSequences:
    input:
        bed = "data/tes_{status}_fam.bed",
        genome = config["GENOME"]
    output:
        "data/seq_{status}.fa"
    shell:
        "bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output} -s"

rule getPolyA:
    input:
        "data/seq_{status}.fa"
    output:
        "data/polyA_{status}.bed"
    threads: max(workflow.cores / 2, 1)
    shell:
        "python3 scripts/annotate2.py {input} {output} {threads}"

rule renamePolByTE:
    input:
        "data/polyA_{status}.bed",
        "data/tes_{status}_fam.bed"
    output:
        temp("data/polyA_{status}.renamed.bed")
    shell:
        "Rscript scripts/renamePolyA.R {input} {output}"

rule getPosPolyA:
    input:
        aoe = config["AOE"],
        polyA = "data/polyA_{status}.renamed.bed"
    output:
        "data/polyA_{status}.tsv"
    shell:
        """
        bin/countFeatures.bin +a {input.aoe} +b {input.polyA} +o {output} +m +i 1000000 +p start +d
        """

# rule collect_polyA:
#     input:
#         expand("data/polyA_{{status}}.{te_type}.tsv", te_type = config["FAMS"])
#     output:
#         "data/polyA_{status}.tsv"
#     shell:
#         "cat {input} > {output}"
rule changeID:
    input:
        "data/tes_{status}_fam.bed"
    output:
        temp("data/tes_{status}_fam_uniqID.bed")
    shell:
        "Rscript scripts/createBedtoolsID.R {input} {output}"

rule getPosCopySpecific:
    input:
        aoe = config["AOE"],
        tes = "data/tes_{status}_fam_uniqID.bed"
    output:
        "data/cpy_pos_{status}.tsv"
    shell:
        "bin/countFeatures.bin +a {input.aoe} +b {input.tes} +o {output} +m +i 1000000 +p stop +s +d"

rule getPosTESpecific:
    input:
        aoe = config["AOE"],
        tes = "data/tes_{status}_fam.bed"
    output:
        "data/pos_{status}.tsv"
    shell:
        """
        bin/countFeatures.bin +a {input.aoe} +b {input.tes} +o {output} +m +i 1000000 +p stop +s +d
        """ # final output is pos val type

rule selectTEbyPos:
    input:
        polyA = "data/cpy_pos_{status}.tsv",
        tes = "data/tes_{status}_fam.bed"
    output:
        "data/first_pos_{status}.bed"
    shell:
        "Rscript scripts/selectAlusByPos.R {input.polyA} {input.tes} {output} -50 50"

rule getNCpG:
    input:
        config["ANC_GEN"]
    output:
        CpG = "data/CpG.bed",
        nCpG = "data/nCpG.bed"
    shadow: "shallow"
    shell:
        "bin/getPattern.bin +f {input} +1 {output.CpG} +2 {output.nCpG} +p CG"

rule preIntersect:
    input:
        tes = "data/first_pos_{status}.bed",
        mask = "data/nCpG.bed"
    output:
        "data/nCG_{status}.bed"
    shell:
        "bin/intersectKeepingNames.bin +a {input.tes} +b {input.mask} +o {output}"

rule count_bases:
    input:
        aoe = config["AOE"],
        genome = config["ANC_GEN"],
        mask = "data/nCG_{status}.bed"
    output:
        "data/bases_{status}_first_pos.tsv"
    shell:
        "bin/countBases.bin +f {input.genome} +a {input.aoe} +b {input.mask} +o {output} +i +m hit"

rule getProfileGC:
    input:
        "data/bases_nCpG.{status}.{te_type}_first_pos.tsv"
    output:
        "data/gc_{status}.{te_type}.tsv"
    params:
        "{te_type}"
    shell:
        """
        Rscript scripts/getRealGCFromCounts.R {input} {output}
        sed -i "s/$/\t{params}/g" {output}
        """

rule collect_gc:
    input:
        expand("data/gc_{{status}}.{te_type}.tsv", te_type = config["FAMS"])
    output:
        "data/gc_{status}.tsv"
    shell:
        "cat {input} > {output}"

rule count_mutations:
    input:
        aoe = config["AOE"],
        muts = config["MUTATIONS"],
        mask = "data/{mask}.{status}.{te_type}.bed"
    output:
        "data/muts_{mask}.{status}.{te_type}_first_pos.tsv"
    shell:
        "bin/countMuts.bin +v {input.muts} +a {input.aoe} +b {input.mask} +o {output}"

rule getRate:
    input:
        "data/bases_{mask}.{status}.{te_type}_first_pos.tsv",
        "data/muts_{mask}.{status}.{te_type}_first_pos.tsv"
    output:
        "data/rates_{mask}.{status}.{te_type}_first_pos.tsv"
    params:
        "{te_type}"
    shell:
        """
        Rscript scripts/rateSW.R {input} {output}
        sed -i "s/^/{params}\t/g" {output}
        """

rule gillespiage:
    input:
        CpG_bases = "data/bases_CpG.{status}.{te_type}_first_pos.tsv",
        CpG_muts = "data/muts_CpG.{status}.{te_type}_first_pos.tsv",
        nCpG_bases = "data/bases_nCpG.{status}.{te_type}_first_pos.tsv",
        nCpG_muts = "data/muts_nCpG.{status}.{te_type}_first_pos.tsv"
    output:
        "data/equilibrium_{status}.{te_type}.tsv"
    shadow: "shallow"
    threads: workflow.cores  / 2
    params:
        "{te_type}"
    shell:
        """
        Rscript scripts/conformMutsV2.R {input.CpG_bases} {input.CpG_muts} muts_CpG.tsv
        Rscript scripts/conformMutsV2.R {input.nCpG_bases} {input.nCpG_muts} muts_nCpG.tsv
        Rscript scripts/simplGillespie.R muts_nCpG.tsv muts_CpG.tsv {output} {threads}
        sed -i "s/^/{params}\t/g" {output}
        """

rule collect_gillespie:
    input:
        expand("data/equilibrium_{{status}}.{te_type}.tsv", te_type = config["FAMS"])
    output:
        "data/equilibrium_{status}.tsv"
    shell:
        "cat {input} > {output}"

rule collect_muts:
    input:
        expand("data/rates_nCpG.{{status}}.{te_type}_first_pos.tsv", te_type = config["FAMS"])
    output:
        "data/muts_{status}.tsv"
    shell:
        "cat {input} > {output}"

rule plot:
    input:
        "data/gc_{status}.tsv",
        "data/polyA_{status}.tsv",
        "data/muts_{status}.tsv",
        "data/equilibrium_{status}.tsv"
    output:
        "results/gc_tes_{status}.svg",
        "results/polyA_tes_{status}.svg",
        "results/muts_tes_{status}.svg",
        "results/equilibrium_gc_{status}.svg"
    shell:
        "Rscript scripts/plot_all.R {input} {output}"

rule plot_all:
    input:
        "results/equilibrium_gc_full.svg",
        "results/equilibrium_gc_partial.svg"
    output:
        touch("done")
