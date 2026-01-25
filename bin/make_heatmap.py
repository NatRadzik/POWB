#!/usr/bin/env python3
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def main():
    if len(sys.argv) < 2:
        print("Uzycie: make_heatmap.py <plik_gene_count>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = "heatmap_genes.png"

    # 1. Wczytanie (FeatureCounts ma nagłówek z #, pomijamy go, ale zachowujemy nazwy kolumn)
    # Zazwyczaj drugi wiersz to nagłówki
    try:
        df = pd.read_csv(input_file, sep="\t", comment="#", index_col=0)
    except Exception as e:
        print(f"Blad wczytywania: {e}")
        sys.exit(1)

    # 2. Usuwamy kolumny techniczne (FeatureCounts zawsze je dodaje)
    cols_to_drop = ['Chr', 'Start', 'End', 'Strand', 'Length']
    # Usuwamy tylko te, które faktycznie istnieją w pliku
    df_counts = df.drop(columns=[c for c in cols_to_drop if c in df.columns], errors='ignore')

    # 3. Wybieramy TOP 30 genów o największej ekspresji
    # Żeby wykres był czytelny
    if df_counts.empty:
        print("Pusta tabela zliczen!")
        sys.exit(1)

    # Sumujemy wiersze i bierzemy 30 największych
    top_genes = df_counts.sum(axis=1).nlargest(30).index
    df_plot = df_counts.loc[top_genes]

    # 4. Logarytmizacja (log2(x+1)) dla lepszej skali kolorów
    df_log = np.log1p(df_plot)

    # 5. Rysowanie
    plt.figure(figsize=(10, 12))
    
    # Jeśli mamy więcej niż 1 próbkę -> robimy klastrowanie (drzewko)
    if df_log.shape[1] > 1:
        g = sns.clustermap(df_log, cmap="viridis", annot=True, fmt=".1f", figsize=(12, 12))
        plt.title("Top 30 Genes (Log Scale)")
        g.savefig(output_file)
    else:
        # Jak jest 1 próbka, clustermap nie zadziała, robimy zwykłą heatmapę
        sns.heatmap(df_log, cmap="viridis", annot=True, fmt=".1f")
        plt.title("Top 30 Genes (Log Scale)")
        plt.savefig(output_file)

    print(f"Wygenerowano: {output_file}")

if __name__ == "__main__":
    main()