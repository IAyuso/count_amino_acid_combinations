# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 16:29:46 2023

@author: ivay
"""

import pandas as pd
import csv
import matplotlib.pyplot as plt


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

    def __init__(self, file_path):
        self.df_alignment = pd.read_csv(file_path, names=['AccNum', 'Seq'])
        self.positions = {252: 'R178', 257: 'D180', 261: 'E183', 263: 'F185'}

    def split_seq_to_chars(self):
        self.df_alignment_split = pd.concat([self.df_alignment['AccNum'],
                                             self.df_alignment['Seq'].astype(str).apply(lambda x: pd.Series(list(x))).astype(str)],
                                            axis=1)

    def rename_columns(self):
        self.df = self.df_alignment_split.rename(columns=self.positions)[['AccNum', *self.positions.values()]]

    def find_unique_sequences(self):
        df_uniques = self.df.groupby(list(self.positions.values())).first().reset_index()
        df_uniques2 = self.df.groupby(list(self.positions.values())).size().reset_index(name='solutions')
        self.df_final = pd.merge(df_uniques, df_uniques2, on=list(self.positions.values()))
        self.df_final = self.df_final[['AccNum', *self.positions.values(), 'solutions']]

    def plot_solutions(self):
        # Create a pie chart using the 'solutions' column of df_final
        fig, ax = plt.subplots()
        ax.pie(self.df_final['solutions'], labels=self.df_final['AccNum'], autopct='%1.1f%%')
        ax.set_title('Solutions by accession number')

        # Show the plot
        plt.show()

    def save_output_to_excel(self, output_path):
        self.df_final.to_excel(output_path, index=False)


if __name__ == "__main__":
    fasta_file = FastaFile("input.fasta")
    fasta_file.write_csv("alignment.csv")
    dna_alignment = AASequenceAlignment('alignment.csv')
    dna_alignment.split_seq_to_chars()
    dna_alignment.rename_columns()
    dna_alignment.find_unique_sequences()
    dna_alignment.plot_solutions()
    dna_alignment.save_output_to_excel('solutions.xlsx')