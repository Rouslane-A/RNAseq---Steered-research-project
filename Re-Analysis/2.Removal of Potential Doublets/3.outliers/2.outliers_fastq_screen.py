import pandas as pd

# Step 1: Read names from the text file
with open('output.txt', 'r') as file:
    names = [line.strip() for line in file]

# Step 2: Read data from the input file
input_file = 'raw_counts.csv'  # Change this to your input file name
output_file = 'output_data.csv'  # Change this to your output file name
data = pd.read_csv(input_file)

# Step 3-4: Keep only the columns that match the names from the text file
columns_to_keep = [data.columns[0]]  # Start with the first column
for col in data.columns[1:]:  # Start checking from the second column
    for name in names:
        if col.startswith(name):
            columns_to_keep.append(col)
            break  # Once a match is found, no need to check further

filtered_data = data[columns_to_keep]

# Step 5: Write the resulting DataFrame to a different file
filtered_data.to_csv(output_file, index=False)

