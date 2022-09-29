#%%
import subprocess
from pathlib import Path


def bash(cmd: str, cwd: Path) -> list[str]:
    """Run a bash command and return the output"""
    out = []
    with subprocess.Popen(
        ["/bin/bash", "-c", cmd],
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ) as process:
        for line in process.stdout:  # type: ignore
            s = line.decode("utf-8")
            print(s)
            out.append(s)
    return out


# %%
geo = "GSM1939675"
currDir = Path.cwd()
dataDir = currDir / "data" / geo

raw = dataDir / "raw"
raw.mkdir(parents=True, exist_ok=True)

#%%
bash(f"ffq -o {geo}.json {geo}", cwd=raw)
print("Gotten metadata")

# %%
bash(
    f"""jq -r ' .[].samples | .[].experiments | .[].runs | .[].files.ftp | .[].url ' {geo}.json |
    head -2 |
    xargs -n 1 -P 2 wget -c""",
    cwd=raw,
)
# %%
print("Getting sample metadata...")
bash(
    f"""jq -r ' .[].samples | .[] | [.accession, (.experiments | .[].runs | .[] | .files.ftp | .[] | .[]), (.attributes | .[])] | @tsv' {geo}.json |
    head -2 >{geo}.tsv""",
    cwd=raw,
)
print("Gotten sample metadata")

#%%
bash("echo $(ls *.fastq.gz) | xargs -n 1 -P 4 fastqc", cwd=raw)

#%%
kal = dataDir / "metadata" / "kallisto_index"
kal.mkdir(parents=True, exist_ok=True)

bash(
    f"""IDX_URL="https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/mus_musculus.tar.gz"
wget -c $IDX_URL
tar -xzvf mus_musculus.tar.gz
""",
    cwd=kal,
)

#%%
results = dataDir / geo / "results"
results.mkdir(parents=True, exist_ok=True)

for fastq in raw.iterdir():
    if fastq.name.endswith(".fastq.gz"):
        sample_name = fastq.stem  # filename without extension
        print(f"Performing kallisto pseudoalignment for {sample_name}.")
        bash(
            f"kallisto quant \
                -i {kal}/mus_musculus/transcriptome.idx \
                -o {sample_name} \
                --single -l 200 -s 20 \
                --threads=8 \
                {fastq} \
                2> >(tee {sample_name}.log)",
            cwd=results,
        )
#%%
bash("multiqc .", cwd=dataDir)

# %%
