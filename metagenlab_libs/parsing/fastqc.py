#!/usr/bin/env python


class FastQCParser():
    
    def __init__(self, fastqc_file):
        import io
        
        if isinstance(fastqc_file, io.IOBase):
            self.modules = fastqc_file.read().split(">>END_MODULE")[0:-1]
        else:
            self.modules = open(fastqc_file, "r").read().split(">>END_MODULE")[0:-1]

        parsed_modules = [self._parse_module(n) for n in range(0, len(self.modules))]
        
        self.fastqc_modules = {i[0]:i[1] for i in parsed_modules if i is not None}
        
    def _parse_module(self, index):
        import pandas 
        module_data = self.modules[index].split("\n")
        module_name = module_data[1].replace('>>', '').split("\t")[0]
        header = False
        values = []
        # last row is always empty
        # last comment line is a table header
        for n, line in enumerate(module_data):
            if line.startswith("#") or line.startswith(">>") or len(line) == 0:
                pass
            else:
                if not header:
                    header = module_data[n-1].replace('#', '').split("\t")
                values.append(line.split("\t"))
        
        if len(values) > 0:
            df = pandas.DataFrame(values, columns=header)
            # set first column as index
            df = df.set_index(df.columns[0])
            return module_name, df
        else:
            return None