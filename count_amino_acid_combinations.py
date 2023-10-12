@author: ivay
"""

import pandas as pd
import csv


class FastaFile:
    def __init__(self, file_path):
        self.file_path = file_path
        self.sequences = {}
        self.read_file()

    def read_file(self):
        sequence = ""
        name = ""
        with open(self.file_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # Store the previous sequence and start a new one
                    if sequence and name:
                        self.sequences[name] = sequence
                    name = line.strip()[1:]
                    sequence = ""
                else:
                    sequence += line.strip()

            # Store the last sequence
            if sequence and name:
                self.sequences[name] = sequence

    def write_csv(self, csv_file):
        with open(csv_file, "w", newline='') as f:
            writer = csv.writer(f)
            for name, sequence in self.sequences.items():
                if not sequence:  # Skip sequences with no data
                    continue
                writer.writerow([name, sequence])

class AASequenceAlignment:

    def __init__(self, file_path, positions_dict):
        self.df_alignment = pd.read_csv(file_path, names=['AccNum', 'Seq'])
        self.positions = positions_dict

    def split_seq_to_chars(self):
        self.df_alignment_split = pd.concat([self.df_alignment['AccNum'],
                                             self.df_alignment['Seq'].astype(str).apply(lambda x: pd.Series(list(x))).astype(str)],
                                            axis=1)

    def rename_columns(self):
        self.df = self.df_alignment_split.rename(columns=self.positions)[['AccNum', *self.positions.values()]]

    def find_unique_sequences(self):
        df_uniques = self.df.groupby(list(self.positions.values()))['AccNum'].apply(list).reset_index()
        df_uniques['solutions'] = df_uniques['AccNum'].apply(len)
        self.df_final = df_uniques[['AccNum', *self.positions.values(), 'solutions']]

    def get_instance_names(self):
        instance_names = {}
        for _, row in self.df_final.iterrows():
            solutions = row['AccNum']
            if solutions:
                acc_num = solutions[0]  # Only consider the first instance in each solution
                if acc_num not in instance_names:
                    instance_names[acc_num] = row['AccNum']
                    instance_names[acc_num + '_criteria'] = ''.join(row[1:-1].astype(str))  # Extract criteria

        return instance_names

    def save_instance_names(self, output_txt):
        instance_names = self.get_instance_names()

        with open(output_txt, 'w') as f:
            for acc_num, solutions in instance_names.items():
                if '_criteria' in acc_num:
                    f.write(f'Criteria: {solutions}\n\n')  # Write criteria section
                else:
                    f.write(f'AccNum: {acc_num}\n')
                    f.write('Solutions:\n')  # New line for Solutions header
                    for solution in solutions:
                        f.write(f'{solution}\n')  # Write each solution on a new line
                    f.write('\n')  # Add an extra new line for separation

    def save_output_to_excel(self, output_path):
        self.df_final.to_excel(output_path, index=False)

if __name__ == "__main__":
    # Take user input for positions dictionary
    positions_dict = {}
    print("Enter amino acid names and corresponding positions in the alignment:")
    while True:
        pos_name = input("Amino acid name (i.e. W178; press ENTER to quit): ")
        if not pos_name:
            break
        pos_num = int(input("Position in the alignment: ")) - 1
        positions_dict[pos_num] = pos_name

    fasta_file = FastaFile("input.fasta")
    fasta_file.write_csv("alignment.csv")
    dna_alignment = AASequenceAlignment('alignment.csv', positions_dict)
    dna_alignment.split_seq_to_chars()
    dna_alignment.rename_columns()
    dna_alignment.find_unique_sequences()
    dna_alignment.save_output_to_excel('solutions.xlsx')
    dna_alignment.save_instance_names('instance_names.txt')
