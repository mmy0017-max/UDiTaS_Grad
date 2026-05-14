import os
import gzip
import subprocess
import regex
from collections import Counter
from Bio import SeqIO
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

### Configuration ###
import argparse

parser = argparse.ArgumentParser(description="UDiTaS Analysis Pipeline")
parser.add_argument(
    "-u",
    "--umi",
    default="MN124/PSP_Mix2_S6_L001_R2_001.fastq",
    help="Path to UMI FASTQ file (e.g., R2)",
)
parser.add_argument(
    "-t",
    "--target",
    default="MN124/PSP_Mix2_S6_L001_R3_001.fastq",
    help="Path to Target FASTQ file (e.g., R3)",
)
parser.add_argument(
    "-x",
    "--genome",
    default="GRCh38_clean/hg38_clean",
    help="Path to Bowtie2 genome index prefix",
)
parser.add_argument(
    "-o", "--outdir", default="uditas_output", help="Output directory path"
)
parser.add_argument(
    "-p", "--threads", type=int, default=8, help="Number of threads for Bowtie2"
)
parser.add_argument(
    "--primer", default="GATCTTGTGGGACCCCTCTGTCCAG", help="Primer sequence"
)
parser.add_argument("--origin", default="CCCAGCCTGGGTGTGCATCT", help="Origin sequence")
parser.add_argument(
    "--chr", default="Auto", help="Target chromosome (e.g., chr1) or Auto"
)

# Parse arguments (sys.argv)
args, unknown = parser.parse_known_args()

FASTQ_UMI_FILE = args.umi
FASTQ_TARGET_FILE = args.target
GENOME_INDEX = args.genome

# Settings
PRIMER_SEQ = args.primer
ORIGIN_SEQ = args.origin  # Set to "" if there is no origin sequence

TARGET_CHR = args.chr  # Set to Chr# for target chromosome, Auto to detect from SAM file
MAX_MISMATCH_PRIMER = 1
LARGE_DELETION_THRESHOLD = 45
SOFT_CLIP_BUFFER = 5  # Added to above value to calculate soft-clip ignore value
SOFT_CLIP_IGNORE_THRESHOLD = LARGE_DELETION_THRESHOLD + SOFT_CLIP_BUFFER
FASTQ_WRITE_BUFFER = 10000  # Process I/O in chunks

# Bowtie2 Settings
BOWTIE2_CMD = "bowtie2"
BOWTIE2_THREADS = args.threads  # Match the number of cores available

# Output paths
OUTPUT_DIR = args.outdir
OUTPUT_FASTQ = os.path.join(OUTPUT_DIR, "trimmed.fastq")
OUTPUT_SAM = os.path.join(OUTPUT_DIR, "trimmed.sam")
OUTPUT_TSV = os.path.join(OUTPUT_DIR, "align_summary.tsv")
OUTPUT_PNG = os.path.join(OUTPUT_DIR, "mutation_visualization.png")


def ensure_dirs():
    os.makedirs(OUTPUT_DIR, exist_ok=True)


### Functions ###


def fastq_read_generator(umi_read_path, target_read_path):
    if not os.path.exists(umi_read_path):
        print(f"error: UMI read file not found -> {umi_read_path}")
        return None
    if not os.path.exists(target_read_path):
        print(f"error: target read file not found -> {target_read_path}")
        return None

    open_func_umi = gzip.open if str(umi_read_path).endswith(".gz") else open
    open_func_target = gzip.open if str(target_read_path).endswith(".gz") else open

    try:
        with open_func_umi(umi_read_path, "rt") as h_umi, open_func_target(
            target_read_path, "rt"
        ) as h_target:
            umi_iter = SeqIO.parse(h_umi, "fastq")
            target_iter = SeqIO.parse(h_target, "fastq")
            for rec_umi, rec_target in zip(umi_iter, target_iter):
                yield rec_umi, rec_target
    except Exception as e:
        print(f"error: input files parsing failed. {e}")
        return None


def process_umi(umi_seq_str):
    mid = len(umi_seq_str) // 2
    return umi_seq_str[:mid], umi_seq_str[mid:]


def trim_origin(target_record, umi_record, pattern_primer, origin_str):
    seq = str(target_record.seq).upper()
    umi_seq = str(umi_record.seq)

    match = pattern_primer.search(seq)
    if match:
        primer_end = match.end()
        # If origin_str is "", len(origin_str) is 0, so it trims safely at the primer.
        origin_end = primer_end + len(origin_str)

        if len(seq) >= origin_end:
            origin_seq = seq[primer_end:origin_end]
            mismatches = sum(1 for a, b in zip(origin_seq, origin_str) if a != b)
            mut_status = "WT" if mismatches == 0 else f"Mut_{origin_seq}"

            # Ensure we only try to parse UMIs if the sequence is long enough, otherwise generic UMI
            if len(umi_seq) > 4:
                inline_idx, umi = process_umi(umi_seq)
            else:
                inline_idx, umi = "sim", umi_record.id.replace(":", "").replace("_", "")

            new_id = f"{target_record.id}_{mut_status}_{inline_idx}_{umi}"

            trimmed_rec = target_record[origin_end:]
            trimmed_rec.id = new_id
            trimmed_rec.name = ""
            trimmed_rec.description = ""

            return trimmed_rec, mut_status
    return None, None


def run_bowtie2(trimmed_fastq, output_sam, genome_index, threads):
    print(f"aligning bait sequences with Bowtie2. (reference: {genome_index})")
    cmd = f'{BOWTIE2_CMD} -x "{genome_index}" -U "{trimmed_fastq}" -S "{output_sam}" -p {threads} --local'
    try:
        subprocess.run(cmd, check=True, shell=True)
        print("alignment complete.")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error running Bowtie2 alignment: {e}")
    except FileNotFoundError:
        raise RuntimeError("Error: Bowtie2 executable not found. Check installation")


# Analyze bait align res. & mut. classification (+ UMI & SV detection)
def analyze_bait_align_res(
    sam_file,
    summary_file,
    large_del_threshold=LARGE_DELETION_THRESHOLD,
    soft_clip_ignore_threshold=SOFT_CLIP_IGNORE_THRESHOLD,
):
    print(
        f"analyzing alignment results... (large deletion threshold: >{large_del_threshold}bp, soft-clip ignore: <={soft_clip_ignore_threshold}bp)"
    )

    stats = {
        "Substitutions": 0,
        "Small_Deletions": 0,
        "Large_Deletions": 0,
        "Insertions": 0,
        "Translocations": 0,
        "Inversions": 0,
        "WT": 0,
    }

    seen_molecules = set()
    pcr_duplicates = 0
    target_chr = TARGET_CHR
    chr_counts = Counter()

    with open(sam_file, "r") as f_sam, open(summary_file, "w") as f_out:
        f_out.write(
            "read_header\torigin_mut\tbait_align_info\tbait_mut\tmut_summary\tmut_category\n"
        )

        lines_buffer = []

        for line in f_sam:
            if line.startswith("@"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue

            flag = int(fields[1])
            if flag & 4:
                continue  # Unmapped

            # Single-pass auto-detect chromosome logic
            if target_chr == "Auto":
                chr_counts[fields[2]] += 1
                lines_buffer.append(fields)
                if len(lines_buffer) > 5000:
                    target_chr = chr_counts.most_common(1)[0][0]
                    print(f"auto-detected target chromosome: {target_chr}")
                    for buf_fields in lines_buffer:
                        _process_sam_line(
                            buf_fields,
                            target_chr,
                            stats,
                            seen_molecules,
                            f_out,
                            large_del_threshold,
                            soft_clip_ignore_threshold,
                        )
                    lines_buffer = []
                continue
            else:
                if lines_buffer:
                    for buf_fields in lines_buffer:
                        _process_sam_line(
                            buf_fields,
                            target_chr,
                            stats,
                            seen_molecules,
                            f_out,
                            large_del_threshold,
                            soft_clip_ignore_threshold,
                        )
                    lines_buffer = []
                _process_sam_line(
                    fields,
                    target_chr,
                    stats,
                    seen_molecules,
                    f_out,
                    large_del_threshold,
                    soft_clip_ignore_threshold,
                )

        # Flush remaining buffer if there were < 5000 mapped reads
        if target_chr == "Auto" and lines_buffer:
            if chr_counts:
                target_chr = chr_counts.most_common(1)[0][0]
                print(f"auto-detected target chromosome: {target_chr}")
            for buf_fields in lines_buffer:
                _process_sam_line(
                    buf_fields,
                    target_chr,
                    stats,
                    seen_molecules,
                    f_out,
                    large_del_threshold,
                    soft_clip_ignore_threshold,
                )

    print(f"tsv saved. deduplication removed {pcr_duplicates} PCR duplicates.")
    print(f"final molecular counts: {stats}")
    return stats


def _process_sam_line(
    fields,
    target_chr,
    stats,
    seen_molecules,
    f_out,
    large_del_threshold,
    soft_clip_ignore_threshold,
):
    header = fields[0]
    chrom = fields[2]
    pos = fields[3]
    cigar = fields[5]

    umi = header.split("_")[-1]
    molecule_signature = (chrom, pos, umi)

    if molecule_signature in seen_molecules:
        return  # PCR Duplicate
    seen_molecules.add(molecule_signature)

    origin_mut = "Mut" if "_Mut_" in header else "WT"
    mut_category = "WT"
    bait_mut = "None"

    if chrom != target_chr:
        mut_category = "Translocations"
        bait_mut = f"Mapped_to_{chrom}"
    else:
        deletions = list(map(int, regex.findall(r"(\d+)D", cigar)))
        max_del = max(deletions) if deletions else 0

        soft_clips = list(map(int, regex.findall(r"(\d+)S", cigar)))
        meaningful_clips = [s for s in soft_clips if s > soft_clip_ignore_threshold]
        max_clip = max(meaningful_clips) if meaningful_clips else 0

        if max_del > large_del_threshold or max_clip > large_del_threshold:
            mut_category = "Large_Deletions"
            bait_mut = f"{max(max_del, max_clip)}bp_Del_or_SV"
        elif "I" in cigar:
            mut_category = "Insertions"
            bait_mut = "Insertion_Detected"
        elif "D" in cigar:
            mut_category = "Small_Deletions"
            bait_mut = f"{max_del}bp_Deletion"
        elif origin_mut == "Mut":
            mut_category = "Substitutions"
            bait_mut = "Mismatch_Only"

    stats[mut_category] += 1
    mut_summary = f"{origin_mut}+{mut_category}"
    f_out.write(
        f"{header}\t{origin_mut}\t{cigar}\t{bait_mut}\t{mut_summary}\t{mut_category}\n"
    )


# Visualization (Pie Chart + Schematic Image)
def visualize_results(stats, output_plot, origin_str="TARGET"):
    print("creating pie chart and sequence alignment schematic.")
    if not stats:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 8), gridspec_kw={'width_ratios': [1.2, 1.5]})

    # Pie Chart
    labels = [k for k, v in stats.items() if v > 0]
    sizes = [v for v in stats.values() if v > 0]
    colors_map = {
        "WT": "#A5D6A7",
        "Small_Deletions": "#90CAF9",
        "Large_Deletions": "#64B5F6",
        "Insertions": "#FFCC80",
        "Substitutions": "#EF9A9A",
        "Translocations": "#CE93D8",
        "Inversions": "#FFAB91",
    }
    colors = [colors_map.get(l, "#CCCCCC") for l in labels]

    if sizes:
        ax1.pie(
            sizes,
            labels=[l.replace("_", " ") for l in labels],
            autopct="%1.1f%%",
            startangle=140,
            colors=colors,
            textprops={"fontsize": 12},
        )
    else:
        ax1.text(0.5, 0.5, "No Data Available", ha="center", va="center", fontsize=14)
    ax1.set_title(
        "Mutation Categories (Unique Molecules)", fontsize=16, fontweight="bold", pad=30
    )

    # Sequence Alignment Schematic
    ax2.axis("off")
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_title(
        "Alignment Patterns & Cleavage Position", fontsize=16, fontweight="bold", pad=30
    )

    safe_origin = origin_str if len(origin_str) >= 5 else "N" * 20
    ref_seq = f"{safe_origin}-TGCTTGGTCGGCACTGATAG"[:40]
    wt_seq = ref_seq
    sm_del_seq = ref_seq[:10] + "-" * 10 + ref_seq[20:]
    lg_del_seq = ref_seq[:8] + "-" * 18 + ref_seq[26:]
    ins_seq = ref_seq[:20] + "A" + ref_seq[20:]
    sub_seq = ref_seq[:20] + "-Ta" + ref_seq[23:]

    x_start = 0.05
    x_step = 0.015
    y_start = 0.75
    y_step = 0.12

    cut_idx = len(safe_origin)
    cut_x = x_start + cut_idx * x_step

    ax2.axvline(
        x=cut_x, ymin=0.05, ymax=0.88, color="gray", linestyle="--", lw=2, zorder=0
    )
    ax2.text(
        cut_x,
        0.90,
        "Predicted cleavage position",
        ha="center",
        va="bottom",
        fontsize=11,
        color="#333333",
    )

    def draw_sequence(ax, seq, y, is_ref=False, highlight_idx=[]):
        for i, char in enumerate(seq):
            if i * x_step > 0.9:
                break  # Prevent running off the image
            x = x_start + i * x_step
            c = "black"
            fw = "normal"

            if is_ref:
                c = "#333333"
            else:
                if char == "-" and (i < len(ref_seq) and ref_seq[i] != "-"):
                    c = "red"
                elif i in highlight_idx:
                    c = "red"
                    fw = "bold"
                    ax.add_patch(
                        patches.Rectangle(
                            (x - x_step * 0.4, y - 0.04),
                            x_step * 0.8,
                            0.08,
                            color="#FFE4E1",
                            zorder=1,
                        )
                    )

            ax.text(
                x,
                y,
                char,
                family="monospace",
                fontsize=13,
                color=c,
                ha="center",
                va="center",
                fontweight=fw,
                zorder=2,
            )

    y_pos = y_start
    draw_sequence(ax2, ref_seq, y_pos, is_ref=True)
    ax2.text(
        x_start + len(ref_seq) * x_step + 0.02,
        y_pos,
        "- Reference",
        va="center",
        fontsize=12,
        fontweight="bold",
        color="#333333",
    )

    sg_start_idx = max(0, len(safe_origin) - 17)
    sg_end_idx = len(safe_origin)
    ax2.hlines(
        y_pos - 0.05,
        x_start + sg_start_idx * x_step,
        x_start + sg_end_idx * x_step,
        color="red",
        lw=2,
    )
    ax2.text(
        x_start + (sg_start_idx + sg_end_idx) / 2 * x_step,
        y_pos - 0.08,
        "sgRNA",
        color="red",
        va="top",
        ha="center",
        fontsize=11,
        fontweight="bold",
    )

    categories = [
        ("WT", wt_seq, []),
        ("Small_Deletions", sm_del_seq, []),
        ("Large_Deletions", lg_del_seq, []),
        ("Insertions", ins_seq, [cut_idx]),
        ("Substitutions", sub_seq, [cut_idx + 1]),
    ]

    y_pos -= y_step + 0.02

    for cat_name, seq, highlights in categories:
        count = stats.get(cat_name, 0)
        draw_sequence(ax2, seq, y_pos, is_ref=False, highlight_idx=highlights)

        display_name = cat_name.replace("_", " ")
        ax2.text(
            x_start + min(len(seq), 50) * x_step + 0.02,
            y_pos,
            f"{display_name} ({count} reads)",
            va="center",
            fontsize=12,
            color="#333333",
        )

        y_pos -= y_step

    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    print(f"visualization result saved: {output_plot}")


### Main Execution ###
if __name__ == "__main__":
    ensure_dirs()

    print(
        f"analyzing target sequences from '{FASTQ_TARGET_FILE}' and UMI from '{FASTQ_UMI_FILE}'"
    )

    read_gen = fastq_read_generator(FASTQ_UMI_FILE, FASTQ_TARGET_FILE)
    re_primer = regex.compile(f"({PRIMER_SEQ}){{e<={MAX_MISMATCH_PRIMER}}}")

    total_reads = 0
    trimmed_reads = 0
    records_buffer = []

    if read_gen:
        with open(OUTPUT_FASTQ, "w") as out_handle:
            for umi_rec, target_rec in read_gen:
                total_reads += 1
                valid_rec, status = trim_origin(
                    target_rec, umi_rec, re_primer, ORIGIN_SEQ
                )

                if valid_rec:
                    records_buffer.append(valid_rec)
                    trimmed_reads += 1

                    if len(records_buffer) >= FASTQ_WRITE_BUFFER:
                        SeqIO.write(records_buffer, out_handle, "fastq")
                        records_buffer = []

                if total_reads % 100000 == 0:
                    print(f"processed {total_reads} reads.")

            if records_buffer:
                SeqIO.write(records_buffer, out_handle, "fastq")

    if total_reads == 0:
        print("no valid reads processed. check input file.")
        exit()

    print(f"total reads: {total_reads}")
    print(f"trimmed reads: {trimmed_reads}")

    run_bowtie2(OUTPUT_FASTQ, OUTPUT_SAM, GENOME_INDEX, BOWTIE2_THREADS)
    stats_result = analyze_bait_align_res(
        OUTPUT_SAM, OUTPUT_TSV, LARGE_DELETION_THRESHOLD, SOFT_CLIP_IGNORE_THRESHOLD
    )

    if stats_result:
        visualize_results(stats_result, OUTPUT_PNG, origin_str=ORIGIN_SEQ)

    print("pipeline completed.")
