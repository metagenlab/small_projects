#!/usr/bin/env python

from docx.api import Document
import pandas as pd

def extract_tables(filename):
    # Load the first table from your document. In your example file,
    # there is only one table, so I just grab the first one.
    document = Document(filename)
    df_list = []
    for n, table in enumerate(document.tables):
        print(f"Parsing table {n}")
        data = []
        keys = None
        for i, row in enumerate(table.rows):
            if i == 0:
                header = [cell.text for cell in row.cells]
            else:
                row_content = [cell.text for cell in row.cells]
                data.append(row_content)

            
        df = pd.DataFrame(data)
        df.columns = header
        df_list.append(df)

    return df_list

def join_df(df_list, join_column):

    for other in df_list[1:len(df_list)]:
        df_list[0] = df_list[0].join(other, lsuffix=join_column, rsuffix=join_column)
    return df_list[0]



if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--input_docx',type=str,help="input docx file")
    parser.add_argument("-j",'--join',type=str,help="join column name")

    args = parser.parse_args()
    tables = extract_tables(args.input_docx)
    if args.join:
        tables = [join_df(tables, args.join)]
    for n, table in enumerate(tables):
        filename = os.path.basename(args.input_docx).split(".")[0]
        table.to_csv(f"{filename}_table_{n}.tsv", sep="\t", index=False)
