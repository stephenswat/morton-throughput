CLANGPP = "clang++-15"
GCC = "g++-13"
UARCHS = ["haswell", "znver3"]
COMPILERS = ["gcc", "clang"]
DIMS = range(1, 11)
LAYOUTS = ["Morton", "Canonical"]
FILES = ["main"]

OSACA_UARCH_NAME = {
  "haswell": "HSW",
  "znver3": "ZEN3"
}

rule all:
  input:
    "out/plots/main.llvm-mca.pdf",
    "out/plots/main.osaca.pdf",

rule plot:
  input:
    table = "out/tables/{file_name}.csv"
  output:
    "out/plots/{file_name}.llvm-mca.pdf",
    "out/plots/{file_name}.osaca.pdf",
  script:
    "scripts/make_plot.py"

rule table:
  input:
    llvm_mca = expand("out/analysis/{file_name}.llvm-mca-{uarch}-{compiler}.json", file_name=FILES, compiler=COMPILERS, uarch=UARCHS),
    osaca = expand("out/analysis/{file_name}.osaca-{uarch}-{compiler}-{layout}-{dim}.yaml", file_name=FILES, compiler=COMPILERS, dim=DIMS, layout=LAYOUTS, uarch=UARCHS)
  output:
    table = "out/tables/{file_name}.csv"
  script:
    "scripts/make_table.py"

rule asm_clang:
  input:
    cxx = "{file_name}.cpp"
  output:
    asm = "out/asm/{file_name}.clang-{uarch}.s"
  shell:
    "{CLANGPP} -mtune={wildcards.uarch} -march={wildcards.uarch} -std=c++17 -O2 -mbmi2 -c -S {input.cxx} -o {output.asm}"

rule gcc_clang:
  input:
    cxx = "{file_name}.cpp"
  output:
    asm = "out/asm/{file_name}.gcc-{uarch}.s"
  shell:
    "{GCC} -mcpu={wildcards.uarch} -mtune={wildcards.uarch} -march={wildcards.uarch} -std=c++17 -O2 -mbmi2 -c -S {input.cxx} -o {output.asm}"

rule llvm_mca_result:
  input:
    asm = "out/asm/{file_name}.{compiler}-{uarch}.s"
  output:
    analysis = "out/analysis/{file_name}.llvm-mca-{uarch}-{compiler}.json"
  shell:
    "llvm-mca-15 --march=x86-64 --mcpu={wildcards.uarch} -json {input.asm} > {output.analysis}"

rule osaca_result:
  input:
    asm = "out/asm/{file_name}.{compiler}-{uarch}.s"
  output:
    analysis = "out/analysis/{file_name}.osaca-{uarch}-{compiler}-{layout}-{dim}.yaml"
  params:
    uarch_name = lambda wildcards, output: OSACA_UARCH_NAME[wildcards.uarch]
  shell:
    "sed -En \"/LLVM-MCA-BEGIN {wildcards.layout}\\\\\\${wildcards.dim}$/,/LLVM-MCA-END {wildcards.layout}\\\\\\${wildcards.dim}$/p\" {input.asm} | osaca --arch {params.uarch_name} --yaml-out {output.analysis} --out {output.analysis}.txt -"
