'''
Looking to make a version of aegis that works on single sequences
'''
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Flatten, Reshape, Dropout
import tensorflow as tf
import matplotlib.pyplot as plt
import os

sequence = "MTSEKGPSTGDPTLRRRIEPWEFDVFYDPRELRKEACLLYEIKWGMSRKIWRSSGKNTTNHVEVNFIKKFTSERDFHPSMSCSITWFLSWSPCWECSQAIREFLSRHPGVTLVIYVARLFWHMDQQNRQGLRDLVNSGVTIQIMRASEYYHCWRNFVNYPPGDEAHWPQYPPLWMMLYALELHCIILSLPPCLKISRRWQNHLTFFRLHLQNCHYQTIPPHILLATGLIHPSVAWR"

# Find unique symbols in amino acids
# Extract unique amino acids
amino_acids = sorted(list(set(list(sequence))))
V = len(amino_acids)
aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
idx_to_aa = {i: aa for aa, i in aa_to_idx.items()}
# Sliding window parameters
WINDOW_SIZE = 20
WINDOW_STEP = 1

# Extract window sub-sequences 
sub_sequences = []
for i in range(0, len(sequence) - WINDOW_SIZE + 1, WINDOW_STEP):
    sub_sequences.append(list(sequence[i:i + WINDOW_SIZE]))

sequences = np.array(sub_sequences)
num_windows, L = sequences.shape  

# Create 3D One-Hot Encoded Dataset:
encoded_sequences = np.zeros((num_windows, L, V), dtype=np.float32)
for i, seq in enumerate(sequences):
    for j, aa in enumerate(seq):
        encoded_sequences[i, j, aa_to_idx[aa]] = 1.0

input_dim = L
num_amino_acids = V

# Split the data into train and test sets
X_train, X_test = train_test_split(encoded_sequences, test_size=0.25, random_state=42)

latent_dim = sequences.shape[1]  # same as input_dim

input_layer = Input(shape=(input_dim, num_amino_acids))
x = Flatten()(input_layer)
x = Dense(latent_dim, activation='relu')(x)
x = Dropout(0.30, name="mc_dropout")(x)   # <-- Stochastic latent for MC

output_layer = Dense(input_dim * num_amino_acids, activation='softmax')(x)
output_layer = Reshape((input_dim, num_amino_acids))(output_layer)

autoencoder = Model(inputs=input_layer, outputs=output_layer)
autoencoder.compile(optimizer='adam', loss='mse')
history = autoencoder.fit(X_train, X_train, epochs=100, batch_size=16, shuffle=True, verbose=1)

# Obtain the latent representations
encoders = Model(inputs=input_layer, outputs=x)
latent_vectors = encoders.predict(encoded_sequences, verbose=0)

T = 200                  # Number of MC passes
eps = 1e-12

mc_pass_means = []       # each item: (L, V) mean over test samples for that pass
for t in range(T):
    preds = autoencoder(X_test, training=True).numpy()     # keep dropout ON
    preds = preds / np.clip(preds.sum(axis=2, keepdims=True), 1e-12, None)
    pass_mean = preds.mean(axis=0)                         # (L, V)
    mc_pass_means.append(pass_mean)

mc_pass_means = np.stack(mc_pass_means, axis=0)            # (T, L, V)

# Aggregate across MC passes
mean_prob_mc = mc_pass_means.mean(axis=0)                  # (L, V)
var_prob_mc  = mc_pass_means.var(axis=0)                   # (L, V)
epi_var_pos  = var_prob_mc.mean(axis=1)                    # (L,) epistemic variance per position
entropy_pos  = -(mean_prob_mc * np.log(np.clip(mean_prob_mc, eps, None))).sum(axis=1)  # (L,)

ALPHA = 0.5     # Laplace/Dirichlet smoothing for transitions
LAMBDA = 0.7    # blend weight: fused ∝ p_MC^λ * p_Markov^(1-λ)

# Per-position transition matrices Tmat[j] shape (V,V): P(next=b | current=a) at position j
Tmat = np.zeros((L, V, V), dtype=float)
for j in range(L):
    prev = sequences[:-1, j]
    nxt  = sequences[1:,  j]
    for a, b in zip(prev, nxt):
        ia, ib = aa_to_idx.get(a), aa_to_idx.get(b)
        if ia is not None and ib is not None:
            Tmat[j, ia, ib] += 1.0
    # smooth and row-normalize
    Tmat[j] = Tmat[j] + ALPHA
    Tmat[j] = Tmat[j] / (Tmat[j].sum(axis=1, keepdims=True) + eps)

def stationary(P, iters=2000, tol=1e-12):
    # left stationary distribution via power iteration
    x = np.ones(P.shape[0], dtype=float) / P.shape[0]
    for _ in range(iters):
        x_new = x @ P
        if np.linalg.norm(x_new - x, 1) < tol:
            break
        x = x_new
    return x / (x.sum() + eps)

# Stationary distribution per position (L,V)
pi_markov = np.stack([stationary(Tmat[j]) for j in range(L)], axis=0)

# Fuse with MC mean distribution computed above
p_mc = mean_prob_mc / (mean_prob_mc.sum(axis=1, keepdims=True) + eps)   # (L,V)
p_mk = pi_markov  / (pi_markov.sum(axis=1, keepdims=True) + eps)        # (L,V)
p_fused = np.power(p_mc + eps, LAMBDA) * np.power(p_mk + eps, 1.0 - LAMBDA)
p_fused = p_fused / (p_fused.sum(axis=1, keepdims=True) + eps)

# Rank fused distribution by small margin + high entropy
order = np.argsort(-p_fused, axis=1)
top1_f = order[:, 0]; top2_f = order[:, 1]
p1f = p_fused[np.arange(L), top1_f]; p2f = p_fused[np.arange(L), top2_f]
margin_f = p1f - p2f
entropy_f = -(p_fused * np.log(np.clip(p_fused, eps, None))).sum(axis=1)

# Normalize & score
def minmax(x):
    lo, hi = float(np.min(x)), float(np.max(x))
    return np.zeros_like(x) if hi <= lo + 1e-12 else (x - lo) / (hi - lo)

score_fused = (1.0 - minmax(margin_f)) * (entropy_f / np.log(V))

TOPK = min(20, L)
rank_idx = np.argsort(-score_fused)[:TOPK]

# Make results directory
os.makedirs("results", exist_ok=True)

fused_rows = []
for j in rank_idx:
    fused_rows.append({
        "Position ": int(j),
        "Top1 (fused)": idx_to_aa[int(top1_f[j])],
        "Top2 (fused)": idx_to_aa[int(top2_f[j])],
        "Top1 Prob (fused)": float(p1f[j]),
        "Top2 Prob (fused)": float(p2f[j]),
        "Margin (fused)": float(margin_f[j]),
        "Entropy (fused)": float(entropy_f[j]),
        "Score (fused)": float(score_fused[j]),
    })

aegis_fused_df = pd.DataFrame(fused_rows)
aegis_fused_df.to_csv("results/aegis_mc_markov_fused_hotspots.csv", index=False)
L, V = mean_prob_mc.shape
pos_axis = np.arange(L)

#oredering amino acide by groups for better understanding
aa_order = [
    # Hydrophobic
    'A','V','I','L','M','F','W','Y',
    # uncharged
    'S','T','N','Q',
    # Charged
    'D','E','K','R','H',
    #ambigous
    'B', 'Z',
    # other
    'G','P','C'
    ]


df = pd.read_csv("results/aegis_mc_markov_fused_hotspots.csv")
top_10_df = df.head(10)
top_10_sorted = top_10_df.sort_values(by="Position ", ascending=True)
subsetA = top_10_sorted["Position "].values

#grabs probabilities
P = p_fused[subsetA, :]  
# creates an index map of amino acids
aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
#double checking only amino acids in data are used (prevents crashing)
valid_aa = [aa for aa in aa_order if aa in aa_to_idx]
#adds any amino acids not in key
extra_aa = [aa for aa in amino_acids if aa not in aa_order]
final_aa = valid_aa + extra_aa
#reordering them to group by type
reorder_idx = [aa_to_idx[aa] for aa in final_aa]
amino_acids_reordered = final_aa
#reorders probabilities
P = P[:, reorder_idx]
#normalizes inter column for better internal change understanding but losing between column comparison
P_norm = P / np.clip(P.max(axis=1, keepdims=True), 1e-12, None)

plt.figure(figsize=(10, 3.5))
plt.imshow(P.T, aspect="auto", origin="lower", interpolation="nearest")
plt.yticks(np.arange(len(amino_acids_reordered)), amino_acids_reordered)
plt.xticks(range(len(subsetA)), [str(pos_axis[i]) for i in subsetA], rotation=0)
#print("mapped positions:", [pos_axis[i] for i in subsetA])
plt.colorbar(label="Probability")
plt.xlabel("Protein residue"); plt.ylabel("Amino acid")
plt.title("Per-AA probabilities at top hotspots (MC Mean)")
plt.tight_layout(); plt.savefig("results/heatmap_probs_top_hotspots.png")

#looking within row variation
plt.figure(figsize=(10, 3.5))
plt.imshow(P_norm.T, aspect="auto", origin="lower", interpolation="nearest")
plt.yticks(np.arange(len(amino_acids_reordered)), amino_acids_reordered)
plt.xticks(range(len(subsetA)), [str(pos_axis[i]) for i in subsetA], rotation=0)
plt.colorbar(label="Relative probability")
plt.xlabel("Protein residue"); plt.ylabel("Amino acid")
plt.title("Per-AA probabilities at top hotspots (MC Mean Normalized)")

plt.tight_layout(); plt.savefig("results/heatmap_probs_normalized.png")
print("\nSaved: aegis_mc_markov_fused_hotspots.csv")