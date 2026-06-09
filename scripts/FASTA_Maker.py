import pandas as pd


Stats = pd.read_csv("results/sequence_check.csv")

MostLikely = Stats["Most Likely AA"].to_list()

sequence = "".join(MostLikely)
#print(str(sequence))

with open("sequence.txt", "w") as f:
    f.write(sequence)
    f.close