#!/usr/bin/env python3
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def main():
    #Sprawdzenie czy podano plik wejściowy w argumentach komendy
    if len(sys.argv) < 2:
        print("Uzycie: make_heatmap.py <plik_gene_count>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = "heatmap_genes.png"

    #Wczytanie danych z tabeli featureCounts
    #Parametr comment wyklucza linie zaczynające się od płatka
    #index_col ustawia pierwszą kolumnę z nazwami genów jako identyfikatory wierszy
    try:
        df = pd.read_csv(input_file, sep="\t", comment="#", index_col=0)
    except Exception as e:
        print(f"Blad wczytywania: {e}")
        sys.exit(1)

    #Czyszczenie tabeli z informacji dodatkowych o lokalizacji genów
    #Zostawia tylko kolumny z liczbą odczytów dla poszczególnych próbek
    cols_to_drop = ['Chr', 'Start', 'End', 'Strand', 'Length']
    df_counts = df.drop(columns=[c for c in cols_to_drop if c in df.columns], errors='ignore')

    # Sprawdzenie czy po filtrowaniu w tabeli zostały jakiekolwiek dane
    if df_counts.empty:
        print("Pusta tabela zliczen!")
        sys.exit(1)

    #Wybór 30 genów o najwyższej sumarycznej ekspresji we wszystkich próbkach
    top_genes = df_counts.sum(axis=1).nlargest(30).index
    df_plot = df_counts.loc[top_genes]

    #Transformacja logarytmiczna danych wzorem log x plus 1
    #Zapobiega to dominacji wykresu przez pojedyncze geny o bardzo wysokich wartościach
    df_log = np.log1p(df_plot)

    #Ustawienie rozmiaru pola wykresu
    plt.figure(figsize=(10, 12))
    
    #Wybór typu wykresu zależnie od liczby analizowanych próbek
    if df_log.shape[1] > 1:
        #Jeśli próbek jest więcej niż jedna stosujemy klastrowanie hierarchiczne
        #clustermap rysuje heatmapę oraz drzewka podobieństwa genów i próbek
        g = sns.clustermap(df_log, cmap="viridis", annot=True, fmt=".1f", figsize=(12, 12))
        plt.title("Top 30 Genes (Log Scale)")
        g.savefig(output_file)
    else:
        #Dla pojedynczej próbki rysuje zwykłą mapę ciepła bez drzewek podobieństwa
        sns.heatmap(df_log, cmap="viridis", annot=True, fmt=".1f")
        plt.title("Top 30 Genes (Log Scale)")
        plt.savefig(output_file)

    #Informacja zwrotna o pomyślnym zakończeniu generowania obrazu
    print(f"Wygenerowano: {output_file}")

if __name__ == "__main__":
    main()
