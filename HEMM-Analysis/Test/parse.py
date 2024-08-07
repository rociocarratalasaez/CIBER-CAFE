import os
import pandas as pd

# Folder where the results are stored
results_folder = '/app/Test/Results'

# Initialize an empty list to store the data
data = []

# Iterate over all files in the results folder
for file_name in os.listdir(results_folder):
    if file_name.startswith('res_') and file_name.endswith('.csv'):
        # Extract parameters from the file name
        parts = file_name.split('_')
        enc = parts[1]
        pre = parts[2]
        lib = parts[3]
        dim = parts[4]
        rep = parts[5].split('.')[0]

        # Read the content of the file
        file_path = os.path.join(results_folder, file_name)
        df = pd.read_csv(file_path, delimiter=';')

        # Assume 'ini' and 'fin' are columns in the CSV file
        if 'ini' in df.columns and 'fin' in df.columns:
            # Calculate the difference between 'fin' and 'ini'
            difference = df['fin'] - df['ini']

            # Create a new DataFrame with the necessary columns
            df_selected = pd.DataFrame({
                'Precision': pre,
                'Encryption': enc,
                'Dimension': dim,
                'Library': lib,
                'Repetition': rep,
                'Difference': difference
            })

            # Append the data to the list
            data.append(df_selected)
        else:
            print(f"File {file_name} does not contain 'ini' and 'fin' columns. Skipping.")

# Concatenate all the dataframes
combined_df = pd.concat(data, ignore_index=True)

# Save the combined dataframe to a new CSV file
output_file = os.path.join(results_folder, 'combined_results.csv')
combined_df.to_csv(output_file, index=False, sep=';')

print(f'All results have been combined into {output_file}')

