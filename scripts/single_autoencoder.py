'''
Looking to make a version of aegis that works on single sequences
'''
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import lightning as L
from torch.utils.data import DataLoader, TensorDataset

sequence = "MTSEKGPSTGDPTLRRRIEPWEFDVFYDPRELRKEACLLYEIKWGMSRKIWRSSGKNTTNHVEVNFIKKFTSERDFHPSMSCSITWFLSWSPCWECSQAIREFLSRHPGVTLVIYVARLFWHMDQQNRQGLRDLVNSGVTIQIMRASEYYHCWRNFVNYPPGDEAHWPQYPPLWMMLYALELHCIILSLPPCLKISRRWQNHLTFFRLHLQNCHYQTIPPHILLATGLIHPSVAWR"

# Find unique symbols in amino acids
# Should be 22 symbols with X included
amino_acids = np.unique(list(sequence))

# Map symbols to integers
amino_acid_to_int = {aa: i for i, aa in enumerate(amino_acids)}

# Find length of sequences (Assumes all to be the same due to MSA) and number of symbols 
num_amino_acids = len(amino_acids)

WINDOW_SIZE = 10
sub_sequences = []
for i in range(len(sequence) - WINDOW_SIZE + 1):
    sub_sequences.append(sequence[i:i+WINDOW_SIZE])

num_sequences = len(sub_sequences)


# This will hold the one-hot encoded sequences
encoded_sequence = np.zeros((num_sequences, WINDOW_SIZE, num_amino_acids), dtype=np.float32)

for i, seq in enumerate(sub_sequences):
    for j, aa in enumerate(seq):
        encoded_sequence[i, j, amino_acid_to_int[aa]] = 1.0

# Reshape into flat vectors for your dense linear layers
flattened_dataset = encoded_sequence.reshape(num_sequences, -1)
input_dimension = flattened_dataset.shape[1]

# Convert to PyTorch DataLoader
dataset = TensorDataset(torch.tensor(flattened_dataset))
dataloader = DataLoader(dataset, batch_size=16, shuffle=True)


class BioAutoEncoder(L.LightningModule):
    def __init__(self, input_dim):
        super().__init__()
        
        # Encoder: 5 layers (Input -> 32 -> 32 -> 32 -> 32 -> 32)
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 32), nn.Sigmoid(), nn.Dropout(0.2),
            nn.Linear(32, 32), nn.Sigmoid(), nn.Dropout(0.2),
            nn.Linear(32, 32), nn.Sigmoid(), nn.Dropout(0.2),
            nn.Linear(32, 32), nn.Sigmoid(), nn.Dropout(0.2),
            nn.Linear(32, 32), nn.Sigmoid(), nn.Dropout(0.2)
        )
        
        # Decoder: 5 layers (32 -> 32 -> 32 -> 32 -> 32 -> Output)
        # Note: The final output shape matches the input dimension
        self.decoder = nn.Sequential(
            nn.Linear(32, 32), nn.Sigmoid(), nn.Dropout(0.2),
            nn.Linear(32, 32), nn.Sigmoid(), nn.Dropout(0.2),
            nn.Linear(32, 32), nn.Sigmoid(), nn.Dropout(0.2),
            nn.Linear(32, 32), nn.Sigmoid(), nn.Dropout(0.2),
            nn.Linear(32, input_dim)
        )

    def forward(self, x):
        z = self.encoder(x)
        return self.decoder(z)

    def training_step(self, batch, batch_idx):
        # Expecting a batch of flattened one-hot vectors
        x, = batch 
        x_hat = self(x)
        
        # MSE loss works well since inputs are one-hot encoded (0s and 1s)
        # and the final layer uses Sigmoid (outputs bounded between 0 and 1)
        loss = F.mse_loss(x_hat, x)
        
        self.log('train_loss', loss, prog_bar=True)
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=1e-3)
    

model = BioAutoEncoder(input_dimension)

# Adding this to help Sigmoid layers pass gradients effectively
def init_weights(m):
    if isinstance(m, nn.Linear):
        torch.nn.init.xavier_uniform_(m.weight)
        if m.bias is not None:
            nn.init.constant_(m.bias, 0.0)

model.apply(init_weights)

trainer = L.Trainer(max_epochs=1000, accelerator="auto")
trainer.fit(model, dataloader)

model.eval()
mutation_predictions = []
with torch.no_grad():
    # Evaluate every sliding window position across the sequence
    for window_idx in range(len(flattened_dataset)):
        sample_input = torch.tensor(flattened_dataset[window_idx:window_idx+1]) 
        reconstruction = model(sample_input)
        
        # Reshape the flat prediction tensor back to: (WINDOW_SIZE, num_amino_acids)
        reconstruction_grid = reconstruction.view(WINDOW_SIZE, num_amino_acids).numpy()
        
        # Pull the one-hot representation for the window
        actual_window_matrix = flattened_dataset[window_idx].reshape(WINDOW_SIZE, num_amino_acids)
        
        # Focus explicitly on the first amino acid position of this specific window context
        actual_aa_idx = actual_window_matrix[0].argmax()
        actual_aa_char = [k for k, v in amino_acid_to_int.items() if v == actual_aa_idx][0]
        
        # Grab the model's 1D probability distribution array for this exact position
        position_scores = reconstruction_grid[0]  # This is a 1D array of length num_amino_acids
        
        # Sort predictions to evaluate alternate choices
        top_indices = np.argsort(position_scores)[::-1]
        predicted_aa_idx = top_indices[0]
        
        # Calculate Reconstruction Error: High error = Unstable / High Mutation Liability
        # position_scores[actual_aa_idx] is now guaranteed to index into a 1D array
        position_error = np.abs(1.0 - position_scores[actual_aa_idx])
        
        # Find the top alternative option predicted by the model's environment context
        alt_aa_idx = top_indices[1] if predicted_aa_idx == actual_aa_idx else top_indices[0]
        alt_aa_char = [k for k, v in amino_acid_to_int.items() if v == alt_aa_idx][0]
        
        mutation_predictions.append({
            "Sequence_Position": window_idx,
            "Wild_Type": actual_aa_char,
            "Reconstruction_Error": float(position_error),
            "Top_Predicted_Mutation": alt_aa_char,
            "Mutation_Score": float(position_scores[alt_aa_idx])
        })

# Parse the final output into a Pandas DataFrame
df_results = pd.DataFrame(mutation_predictions)

# Sort the results so the amino acids MOST likely to mutate are at the top
df_results = df_results.sort_values(by="Reconstruction_Error", ascending=False)

df_results.to_csv("Single_Sequences_Results.csv", index=False)