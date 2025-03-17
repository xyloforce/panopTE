wildcard_constraints:
    status="normal|noTE|shuffle",
    context="nCpG|CpG"

rule all:
    input:
        "results_{species}/equilibrium_gc_full.png",
        "results_{species}/muts_nCpG_full.png",
        "results_{species}/polyA_tes_full.png",
        "results_{species}/gc_tes_full.png"

rule repeatToBed:
    input:
        config["DFAM"]
    output:
        tes = "data_{species}/tes_normal.bed"
    shell:
        "python3 scripts/convertDfamAnnotToBed.py {input} {output}"

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

rule getFamilies:
    input:
        config["DFAM_CLASS"],
        config["DFAM_LENGTH"]
    output:
        "data_{species}/dfam_families.tsv"
    shell:
        "python3 scripts/createFamilyFile.py {input} {output}"

rule teToFam:
    input:
        "data_{species}/tes_{status}.bed",
        config["DFAM"],
        "data_{species}/dfam_families.tsv"
    output:
        temp("data_{species}/tes_{status}_fam.bed")
    shadow: "shallow"
    params: "{status}" # either full or partial
    # wildcard_constraints:
    #     status="full|partial",
    shell:
        """
        python3 scripts/tetypeToFam.py {input} tmp.bed {params}
        grep -Pv "\t\t" tmp.bed > {output}
        """

rule createGenomeFile:
    input:
        config["GENOME"]
    output:
        "data_{species}/genome.genome"
    shell:
        "python3 scripts/createGenomeFile.py {input} {output}"

rule shuffleRepeat:
    input:
        tes = "data_{species}/tes_normal.bed",
        genome = "data_{species}/genome.genome"
    output:
        "data_{species}/tes_shuffle.bed"
    shell:
        "bedtools shuffle -i {input.tes} -g {input.genome} > {output}"

rule NIEBxTE:
    input:
        tes = "data_{species}/tes_{status}_fam.bed",
        nieb = config["NIEB"]
    output:
        temp("data_{species}/niebs_x_tes_{status}.tsv")
    shell:
        "bedtools intersect -wa	-a {input.tes} -b {input.nieb} > {output}"

rule aggNIEBsPTE:
    input:
        "data_{species}/niebs_x_tes_normal.tsv",
        "data_{species}/niebs_x_tes_shuffle.tsv"
    output:
        "data_{species}/merge_tes_to_niebs.csv"
    shell:
        "Rscript scripts/mergeIntersects.R {input} {output}"

rule getSequences:
    input:
        bed = "data_{species}/tes_{status}_fam.bed",
        genome = config["GENOME"]
    output:
        "data_{species}/seq_{status}.fa"
    shell:
        "bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output} -s"

rule getPolyA:
    input:
        "data_{species}/seq_{status}.fa"
    output:
        "data_{species}/polyA_{status}.bed"
    threads: workflow.cores
    shell:
        "python3 scripts/annotate2.py {input} {output} {threads}"

rule renamePolByTE:
    input:
        "data_{species}/polyA_{status}.bed",
        "data_{species}/tes_{status}_fam.bed"
    output:
        temp("data_{species}/polyA_{status}.renamed.bed")
    shell:
        "Rscript scripts/renamePolyA.R {input} {output}"

rule getPosPolyA:
    input:
        aoe = config["AOE"],
        polyA = "data_{species}/polyA_{status}.renamed.bed"
    output:
        "data_{species}/polyA_{status}.tsv"
    shell:
        """
        bin/countFeatures.bin +a {input.aoe} +b {input.polyA} +o {output} +m +i 1000000 +p start +d
        """

rule changeID:
    input:
        "data_{species}/tes_{status}_fam.bed"
    output:
        temp("data_{species}/tes_{status}_fam_uniqID.bed")
    shell:
        "Rscript scripts/createBedtoolsID.R {input} {output}"

rule getPosCopySpecific:
    input:
        aoe = config["AOE"],
        tes = "data_{species}/tes_{status}_fam_uniqID.bed"
    output:
        "data_{species}/cpy_pos_{status}.tsv"
    shell:
        "bin/countFeatures.bin +a {input.aoe} +b {input.tes} +o {output} +m +i 1000000 +p stop +s +d"

rule getPosTESpecific:
    input:
        aoe = config["AOE"],
        tes = "data_{species}/tes_{status}_fam.bed"
    output:
        "data_{species}/pos_{status}.tsv"
    shell:
        """
        bin/countFeatures.bin +a {input.aoe} +b {input.tes} +o {output} +m +i 1000000 +p stop +s +d
        """ # final output is pos val type

rule selectTEbyPos:
    input:
        polyA = "data_{species}/cpy_pos_{status}.tsv",
        tes = "data_{species}/tes_{status}_fam.bed"
    output:
        "data_{species}/first_pos_{status}.bed"
    shell:
        "Rscript scripts/selectAlusByPos.R {input.polyA} {input.tes} {output} -50 50"

rule getNCpG:
    input:
        config["ANC_GEN"]
    output:
        CpG = "data_{species}/CpG.bed",
        nCpG = "data_{species}/nCpG.bed"
    shadow: "shallow"
    shell:
        "bin/getPattern.bin +f {input} +1 {output.CpG} +2 {output.nCpG} +p CG"

rule getRefCounts:
    input:
        config["AOE"]
    output:
        "data_{species}/ref_counts.tsv"
    shell:
        "bin/countFeatures.bin +a {input} +b {input} +o {output}"

rule preIntersect:
    input:
        tes = "data_{species}/first_pos_{status}.bed",
        mask = "data_{species}/{context}.bed"
    output:
        "data_{species}/{context}_{status}.bed"
    shell:
        "bin/intersectKeepingNames.bin +a {input.tes} +b {input.mask} +o {output}"

rule count_bases:
    input:
        aoe = config["AOE"],
        genome = config["ANC_GEN"],
        mask = "data_{species}/{context}_{status}.bed"
    output:
        "data_{species}/bases_{context}_{status}_first_pos.tsv"
    shell:
        "bin/countBases.bin +f {input.genome} +a {input.aoe} +b {input.mask} +o {output} +i +m hit"

rule count_bases_full:
    input:
        aoe = config["AOE"],
        genome = config["ANC_GEN"],
        mask = "data_{species}/first_pos_{status}.bed"
    output:
        "data_{species}/bases_{status}_first_pos.tsv"
    shell:
        "bin/countBases.bin +f {input.genome} +a {input.aoe} +b {input.mask} +o {output} +i +m hit"

rule getProfileGC:
    input:
        "data_{species}/bases_{status}_first_pos.tsv"
    output:
        "data_{species}/gc_{status}.tsv"
    shell:
        "Rscript scripts/getGCFromCountPerID.R {input} {output}"

rule count_mutations:
    input:
        aoe = config["AOE"],
        muts = config["MUTATIONS"],
        mask = "data_{species}/{context}_{status}.bed"
    output:
        "data_{species}/muts_{context}_{status}_first_pos.tsv"
    shell:
        "bin/countMuts.bin +v {input.muts} +a {input.aoe} +b {input.mask} +o {output} +j +m hit"

rule getRate:
    input:
        "data_{species}/bases_{context}_{status}_first_pos.tsv",
        "data_{species}/muts_{context}_{status}_first_pos.tsv"
    output:
        "data_{species}/rates_{context}_{status}_first_pos.tsv"
    shell:
        """
        Rscript scripts/rateSW.R {input} {output}
        """

rule gillespiage:
    input:
        CpG_rates = "data_{species}/rates_CpG_{status}_first_pos.tsv",
        nCpG_rates = "data_{species}/rates_nCpG_{status}_first_pos.tsv"
    output:
        "data_{species}/equilibrium_{status}.tsv"
    shadow: "shallow"
    threads: workflow.cores
    shell:
        """
        Rscript scripts/conformMutsV2.R {input.CpG_rates} muts_CpG.tsv
        Rscript scripts/conformMutsV2.R {input.nCpG_rates} muts_nCpG.tsv
        Rscript scripts/simplGillespie.R muts_nCpG.tsv muts_CpG.tsv {output} {threads}
        """

rule plotPos:
    input:
        "data_{species}/pos_{status}.tsv",
        "data_{species}/ref_counts.tsv"
    output:
        "results_{species}/pos_{status}_by_fam.png",
        "results_{species}/pos_{status}_savestate.csv"
    shell:
        "Rscript scripts/plotPos.R {input} {output}"

rule plotPolyA:
    input:
        "data_{species}/polyA_{status}.tsv",
        "data_{species}/ref_counts.tsv"
    output:
        "results_{species}/polyA_tes_{status}.png",
        "results_{species}/polyA_{status}_savestate.csv"
    shell:
        "Rscript scripts/plotPolyA.R {input} {output}"

rule plotGC:
    input:
        "data_{species}/gc_{status}.tsv"
    output:
        "results_{species}/gc_tes_{status}.png",
        "results_{species}/gc_{status}_savestate.csv"
    shell:
        "Rscript scripts/plotGC.R {input} {output}"
        
rule plotMuts:
    input:
        "data_{species}/rates_{context}_{status}_first_pos.tsv"
    output:
        "results_{species}/muts_{context}_{status}.png",
        "results_{species}/muts_{context}_{status}_savestate.csv"
    shell:
        "Rscript scripts/plotMuts.R {input} {output}"

rule plotEquilibrium:
    input:
        "data_{species}/equilibrium_{status}.tsv"
    output:
        "results_{species}/equilibrium_gc_{status}.png",
        "results_{species}/equilibrium_{status}_savestate.csv"
    shell:
        "Rscript scripts/plotEquilibrium.R {input} {output}"
