import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import re

def parse_fasta_header(header_line):
    parts = header_line.split('#')
    gene_id_full = parts[0].strip('>')
    start = int(parts[1].strip())
    end = int(parts[2].strip())
    direction = int(parts[3].strip())
    match_id = re.search(r'ID=(\d+_\d+)', parts[4])
    gene_label = match_id.group(1) if match_id else gene_id_full.split(' ')[0]
    return {
        'id': gene_id_full,
        'label': gene_label,
        'start': start,
        'end': end,
        'direction': direction
    }

def plot_genes_with_importance(fasta_file_path, importance_scores, output_image_path="gene_map.png"):
    genes_data = []
    current_header = None

    with open(fasta_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    genes_data.append(parse_fasta_header(current_header))
                current_header = line
        if current_header:
            genes_data.append(parse_fasta_header(current_header))

    if len(genes_data) != len(importance_scores):
        raise ValueError("Gene and importance score counts do not match.")

    importance_scores_np = np.array(importance_scores)
    norm_scores = (
        (importance_scores_np - importance_scores_np.min()) /
        (importance_scores_np.max() - importance_scores_np.min())
        if importance_scores_np.max() != importance_scores_np.min()
        else np.full_like(importance_scores_np, 0.5)
    )

    min_coord = min(g['start'] for g in genes_data)
    max_coord = max(g['end'] for g in genes_data)

    fig, ax = plt.subplots(figsize=(15, 4))
    y_center = 0.5
    arrow_height = 0.2
    text_offset = 0.05

    for i, gene in enumerate(genes_data):
        start, end, direction = gene['start'], gene['end'], gene['direction']
        gene_label, importance = gene['label'], norm_scores[i]
#        color = plt.cm.viridis(importance)
        color = plt.cm.jet(importance)


        length = end - start
        head_length = min(max(length * 0.15, 100), length * 0.5)
        head_width = arrow_height * 1.2

        dx = length if direction == 1 else -length
        arrow_start = start if direction == 1 else end

        arrow = patches.FancyArrowPatch(
            (arrow_start, y_center),
            (arrow_start + dx, y_center),
            arrowstyle=f"simple,head_length={head_length / length:.2f},head_width={head_width / arrow_height:.2f},tail_width=0.6",
            mutation_scale=10,
            fc=color, ec='black', linewidth=0.4, zorder=2
        )
        ax.add_patch(arrow)

        label_x = start + length / 2
        ax.text(label_x, y_center + arrow_height / 2 + text_offset, gene_label,
                ha='center', va='bottom', fontsize=9, color='black', rotation=45)

    ax.set_xlim(min_coord - 0.03 * (max_coord - min_coord), max_coord + 0.03 * (max_coord - min_coord))
    ax.set_ylim(y_center - 0.3, y_center + 0.5)
    ax.set_yticks([])
    ax.set_xlabel("Genome Position (bp)", fontsize=12)
    ax.set_title("Gene Map with Importance Scores", fontsize=14, fontweight='bold')
    ax.spines[['top', 'right', 'left', 'bottom']].set_visible(False)

    sm = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=importance_scores_np.min(), vmax=importance_scores_np.max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.02, pad=0.01)
    cbar.set_label("Importance Score", fontsize=11)

    plt.tight_layout()
    plt.savefig(output_image_path, dpi=300)
    plt.close()
