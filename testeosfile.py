import numpy as np
from pathlib import Path

# # Define the relative path to the file
# relative_path =  "eos" / "eosA"

# # Convert the relative path to an absolute path (optional, but good practice)
# absolute_path = relative_path.resolve()

# Load the file using np.loadtxt
data = np.loadtxt('eos/eosA',skiprows=1)

# Print the data to verify
print(data)
from eos import EOS
