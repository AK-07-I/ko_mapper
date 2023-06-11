import subprocess
import argparse
import sys
import pandas as pd

parser = argparse.ArgumentParser(description='Script Description')
parser.add_argument('-i', '--input', type=str, help='Input file')
parser.add_argument('-o', '--output', type=str, help='Output file')
args = parser.parse_args()

if not args.input or not args.output:
    parser.print_help()
    sys.exit(1)

input_file = args.input
output_file = args.output

#Retrieve Kegg orthology
curl_command = ['curl', 'https://rest.kegg.jp/get/br:ko00001']
result = subprocess.run(curl_command, capture_output=True, text=True)

if result.returncode == 0:
    ko00001 = result.stdout
    with open('KEGG_ORTHOLOGY.txt', 'w') as file:
        file.write(ko00001)
    print("Kegg orthology data downloaded successfully")
else:
    print("Error downloading Kegg orthology data")
    print("Error message:")
    print(result.stderr)

cola = []
colb = []
colc = []
cold = []

#Convert Hierarchial kegg orthology to tab delimited format
with open('KEGG_ORTHOLOGY.txt', 'rt') as file:
    for lines in file:
        if lines[0] == 'D':
            cold.append(lines)

with open('KEGG_ORTHOLOGY.txt', 'rt') as file:    
    for lines in file:
        if lines[0] == 'C':
            c = lines
            for linec in file:
                if linec[0] == 'D':
                    colc.append(c)
                elif linec[0] == 'C':
                    c = linec

with open('KEGG_ORTHOLOGY.txt', 'rt') as file:
    for lines in file:
        if lines[0] == 'B':
            if lines[1] == ' ':
                b = lines
                for lineb in file:
                    if lineb[0] == 'D':
                        colb.append(b)
                    elif lineb[0] == 'B':
                        b = lineb

with open('KEGG_ORTHOLOGY.txt', 'rt') as file:
    for lines in file:
        if lines[0] == 'A':
            a = lines
            for linea in file:
                if linea[0] == 'D':
                    cola.append(a)
                elif linea[0] == 'A':
                    a = linea

dfa = pd.DataFrame(cola, columns=['A'])
dfa = dfa.A.str.split(expand=True)
dfa.drop(dfa.iloc[:, 0:1],inplace=True, axis=1)
dfa['Category'] = dfa[dfa.columns[0:]].apply(lambda x: ' '.join(x.dropna().astype(str)),axis=1)

dfb = pd.DataFrame(colb, columns=['B'])
dfb = dfb.B.str.split(expand=True)
dfb.drop(dfb.iloc[:, 0:2],inplace=True, axis=1)
dfb['Subcategory'] = dfb[dfb.columns[0:]].apply(lambda x: ' '.join(x.dropna().astype(str)),axis=1)

dfc = pd.DataFrame(colc, columns=['C'])
dfc = dfc.C.str.split(expand=True)
dfc.drop(dfc.iloc[:, 0:2],inplace=True, axis=1)
dfc['Pathway'] = dfc[dfc.columns[0:]].apply(lambda x: ' '.join(x.dropna().astype(str)),axis=1)

dfd = pd.DataFrame(cold, columns=['D'])
dfd = dfd.D.str.split(expand=True)
dfd.drop(dfd.iloc[:, 0:1],inplace=True, axis=1)
dfd['K0'] = dfd[dfd.columns[0:]].apply(lambda x: ' '.join(x.dropna().astype(str)),axis=1)

df = pd.concat([dfd['K0'],dfc['Pathway'],dfb['Subcategory'],dfa['Category']], axis=1)
df.to_csv('KEGG_Orthology.tsv', sep='\t', index=False)

#Map K0 ids to Kegg
with open(input_file, 'r') as infile:
    ids = [line.strip() for line in infile]

with open('KEGG_Orthology.tsv', 'r') as KEGG, open(output_file, 'w') as outfile:
    for line in KEGG:
        for id in ids:
            if line.startswith(id):
                outfile.write(line)
                break