import pandas as pd

# Load SQANTI3 classification file
df = pd.read_csv("sqanti_classification.txt", sep="\t")


df['polyA_motif'] = df['polyA_motif'].astype(str).str.upper().isin(['TRUE', 'YES', '1', 'T'])
df['within_polya_site'] = df['within_polya_site'].astype(str).str.upper().isin(['TRUE', 'YES', '1', 'T'])


# Thresholds
fsm_tss_threshold = 50
ism_tss_threshold = 20

# Filter FSM and ISM for base full-length
fsm = df[df['structural_category'] == 'full-splice_match']
ism = df[df['structural_category'] == 'incomplete-splice_match']

fsm_full_length = fsm[(fsm['diff_to_TSS'] <= fsm_tss_threshold) & 
                      ((fsm['polyA_motif']) | (fsm['within_polya_site']))]

ism_full_length = ism[(ism['diff_to_TSS'] <= ism_tss_threshold) & 
                      ((ism['polyA_motif']) | (ism['within_polya_site']))]

full_length_transcripts = pd.concat([fsm_full_length, ism_full_length])

# --- Full-length counts per gene ---
full_length_transcripts['full_length_count'] = full_length_transcripts['associated_gene'].map(
    full_length_transcripts.groupby('associated_gene').size()
)

# --- Confidence score ---
def compute_confidence(row):
    score = 2 if row['structural_category'] == 'full-splice_match' else 1
    # +1 TSS proximity
    threshold = fsm_tss_threshold if row['structural_category']=='full-splice_match' else ism_tss_threshold
    if row['diff_to_TSS'] <= threshold:
        score += 1
    # +1 polyA
    if row['polyA_motif'] or row['within_polya_site']:
        score += 1
    return score

full_length_transcripts['full_length_confidence'] = full_length_transcripts.apply(compute_confidence, axis=1)

# --- CAGE confidence ---
if df['within_cage_peak'].dtype == 'object':
    full_length_transcripts['within_cage_peak'] = full_length_transcripts['within_cage_peak'].astype(str).str.upper().isin(['TRUE', 'YES', '1', 'T'])


full_length_transcripts['cage_confidence'] = full_length_transcripts['within_cage_peak'].astype(int)

# --- Total confidence ---
full_length_transcripts['total_confidence'] = full_length_transcripts['full_length_confidence'] + full_length_transcripts['cage_confidence']

# --- Ranking per gene ---
full_length_transcripts['confidence_rank'] = full_length_transcripts.groupby('associated_gene')['total_confidence'] \
   .rank(method='dense', ascending=True).astype(int)

print(f"High-confidence full-length transcripts: {len(full_length_transcripts)}")
full_length_transcripts.to_csv("full_length_FSM_ISM.csv", index=False)