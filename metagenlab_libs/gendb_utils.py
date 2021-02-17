import pandas
from django.conf import settings
from django.db import IntegrityError
import os 

# setup django do be able to access django db models 
import GEN_database.settings as GEN_settings

try:
    settings.configure(INSTALLED_APPS=GEN_settings.INSTALLED_APPS,
                       DATABASES=GEN_settings.DATABASES)

    import django
    django.setup()
except:
    print("django setup failed-- already done?")
    pass

class DB:
    def __init__(self,):
        
        
        db_type = GEN_settings.DB_DRIVER

        if db_type != "sqlite":
            import MySQLdb
            sqlpsw = os.environ['SQLPSW']
            self.conn = MySQLdb.connect(passwd=sqlpsw,
                                        user="root",
                                        host="127.0.0.1",
                                        db="GEN_LIMS")
            self.cursor = self.conn.cursor()
            # placeholder for sql querries (differ between sqlite and mysql)
            self.spl = '%s'
        else:
            import sqlite3
            self.db_path = GEN_settings.SQLITE_DB
            self.conn = sqlite3.connect(self.db_path)
            self.cursor = self.conn.cursor()
            # placeholder for sql querries (differ between sqlite and mysql)
            self.spl = '?'


    def get_fastq_metadata(self, metric_name, index_str=True, analysis_id=False):
        
        if analysis_id:
            analysis_filter = f'and analysis_id={analysis_id}'
        else:
            analysis_filter = ''
        sql = f'''select fastq_id,value from  GEN_fastqfilesmetadata t1
                  inner join GEN_term t2 on t1.term_id=t2.id
                  where t2.name="{metric_name}" {analysis_filter}            
                  union 
                  select t3.fastq_id,t1.value from  GEN_samplemetadata t1
                  inner join GEN_term t2 on t1.term_id=t2.id
                  inner join GEN_fastqtosample t3 on t1.sample_id=t3.sample_id
                  where t2.name="{metric_name}"
                  '''

        #print(sql)
        if index_str:
            return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}
        else:
            return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}

    def get_sample_metadata(self, metric_name, index_str=True):
        sql = f'''select t1.sample_id,t1.value from GEN_samplemetadata t1
                  inner join GEN_term t2 on t1.term_id=t2.id
                  where t2.name="{metric_name}";'''
        print(sql)      
        if index_str:
            return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}
        else:
            return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}

    def get_fastq_id2run_name(self,):

        sql = '''select t1.id,run_name from GEN_fastqfiles t1
                 inner join GEN_runs t2 on t1.run_id=t2.id
        '''

        return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}

    def get_metadata_labels(self,):

        sql = '''select distinct t2.id,t2.name from GEN_fastqfilesmetadata t1
                 inner join GEN_term t2 on t1.term_id=t2.id
                 union 
                 select distinct t2.id,t2.name from GEN_samplemetadata t1
                 inner join GEN_term t2 on t1.term_id=t2.id
        '''
        
        return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}


    def get_molis_id2fastq_id(self,):

        sql = ''' 
        select fastq_id, molis_id from GEN_fastqtosample t1
        inner join GEN_sample t2 on t1.sample_id=t2.id
        '''

        return {i[0]:i[1] for i in self.cursor.execute(sql,).fetchall()}

    def get_molis_id2sample_id_list(self,):

        sql = ''' 
        select molis_id, id from GEN_sample
        '''

        molis_id2sample_list = {}
        for row in self.cursor.execute(sql,).fetchall():
            if row[0] not in molis_id2sample_list:
                molis_id2sample_list[row[0]] = [row[1]]
            else:
                molis_id2sample_list[row[0]].append(row[1])

        return molis_id2sample_list


    def get_sample_df(self,):

        sql = ''' 
        select * from GEN_sample
        '''

        return pandas.read_sql(sql, self.conn)


    def get_fastq_metadata_list(self, 
                                term_list=False, 
                                fastq_filter=None, 
                                run_name_list=False,
                                metadata_value_list=False,
                                add_molis=False):
        '''
        retrieve metadata from both sample and fastq metadata table
        '''

        res_filter = ''
        if term_list:
            print("term_list", term_list)
            term_filter = '","'.join(term_list)
            res_filter += f'and t2.name in ("{term_filter}")\n' 
        if fastq_filter:
            fastq_filter_str = ','.join([str(i) for i in fastq_filter])
            res_filter += f'and fastq_id in ({fastq_filter_str})\n' 
        if run_name_list:
            run_filter = '","'.join(run_name_list)
            res_filter += f'and run_name in ("{run_filter}")'
        if metadata_value_list:
            metadata_filter = '","'.join(metadata_value_list)
            res_filter += f'and t1.value in ("{metadata_filter}")'
        
        
        sql = f'''
            select distinct fastq_id, t2.name,t1.value, run_name from GEN_fastqfilesmetadata t1 
            inner join GEN_term t2 on t1.term_id=t2.id 
            inner join GEN_fastqfiles t3 on t1.fastq_id=t3.id 
            inner join GEN_runs t4 on t3.run_id=t4.id 
            where t3.fastq_prefix not like "Undetermined%"
            {res_filter}
            group by fastq_id, t2.name,t3.fastq_prefix,t4.run_date,t4.run_name,t4.read_length
            union 
            select distinct t3.fastq_id, t2.name,t1.value,run_name from GEN_samplemetadata t1 
            inner join GEN_term t2 on t1.term_id=t2.id 
            inner join GEN_fastqtosample t3 on t1.sample_id=t3.sample_id
            inner join GEN_fastqfiles t4 on t3.fastq_id=t4.id 
            inner join GEN_runs t5 on t4.run_id=t5.id 
            where t4.fastq_prefix not like "Undetermined%"
            {res_filter}
            group by t3.fastq_id, t2.name,t4.fastq_prefix,t5.run_date,t5.run_name,t5.read_length
            '''
        #print(sql)
        # AND t4.run_name like "%_CleanPlex"
        # 
        df = pandas.read_sql(sql, self.conn)

        if add_molis:
            print("adding molis")
            df_molis = self.get_fastq_and_sample_data(df["fastq_id"].to_list()).set_index("fastq_id")
            print(df_molis.head())
            df = df.set_index("fastq_id").join(df_molis, on="fastq_id", rsuffix='_other')
            print("------------------")
            print(df.head())
            print("----------------------")
            df = df[["molis_id", "name", "value", "run_name"]]
            print(df.head())
        return df

    def get_xslx_id2fastq_id(self,):
        
        sql = '''select t2.xlsx_sample_ID,t1.fastq_id from GEN_fastqtosample t1 
                inner join GEN_sample t2 on t1.sample_id=t2.id '''

        return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}


    def get_xslx_id2sample_id(self,):
        
        sql = '''select xlsx_sample_ID,id from GEN_sample'''

        return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}


    def count_qc_warning(self, key_str=False):
        
        df = self.get_fastq_metadata_list()[["fastq_id","value"]]

        term2n_warnings = df.query('value == "WARN"').groupby(["fastq_id"])['value'].count().to_dict()

        if key_str:
            return {str(i):term2n_warnings[i] for i in term2n_warnings}
        else:
            return term2n_warnings




    def get_fastq_metadata_stats(self, term_list):

        df = self.get_fastq_metadata_list(term_list)[["fastq_id","name","value"]]
        
        df["value"] = pandas.to_numeric(df["value"])
        term2median = df.groupby(["name"])['value'].median().round(3).to_dict()
        term2mean = df.groupby(["name"])['value'].mean().round(3).to_dict()
        term2sd = df.groupby(["name"])['value'].std().round(3).to_dict()        

        return term2median, term2mean, term2sd

    def get_fastq_and_sample_data(self, fastq_id_list):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_id_list])
        
        # left join because some fastq won't have match in the sample table
        # possible problem: fastq prefix match with multiple samples from different species
        # in that case: remove species name
        sql = f'''select distinct t1.id as fastq_id,fastq_prefix,R1,R2,species_name,molis_id from GEN_fastqfiles t1 
                left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
                left join GEN_sample t3 on t2.sample_id=t3.id 
                where t1.id in ("{fastq_list_filter}");
            '''
        print(sql,)
        return pandas.read_sql(sql, self.conn)


    def get_analysis_metadata(self, term_name):
        
        sql = f"""select t1.analysis_id,value from GEN_analysismetadata t1 
                  inner join GEN_term t2 on t1.term_id =t2.id 
                  where t2.name like '{term_name}';"""

        return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}

    def parse_CT_scoring_table(self, table_path):
        
        df = pandas.read_csv(table_path, 
                               index_col=None,
                               header=0,
                               sep="\t")
        df = df.fillna(value=False)
        term_name2data = {}
        for n, row in df.iterrows():
            term_name2data[row["name"]] = {"LCL_failed": row["LCL_failed"], "UCL_failed":row["UCL_failed"], "LCL_warning":row["LCL"], "UCL_warning":row["UCL"]}

        return term_name2data

    def qc_score(self, fastq_id_list, config):
        
        metric2scoring = self.parse_CT_scoring_table(config["SCORING"])

        df = self.get_fastq_metadata_list(term_list=list(metric2scoring.keys()), fastq_filter=fastq_id_list)
        # fastq_id, t2.name,t1.value, run_name
        # default score to 0
        fastq_id2n_fail = {i:0 for i in df["fastq_id"].to_list()}
        fastq_id2n_warn = {i:0 for i in df["fastq_id"].to_list()}
        fastq_id2metric2score = {i:{} for i in df["fastq_id"].to_list()}
                
        for n, row in df.iterrows(): 
            fastq_id = row["fastq_id"]

            if row["name"] in metric2scoring:
                # default green
                LCL_failed = metric2scoring[row["name"]]["LCL_failed"]
                UCL_failed = metric2scoring[row["name"]]["UCL_failed"]
                LCL_warning = metric2scoring[row["name"]]["LCL_warning"]
                UCL_warning = metric2scoring[row["name"]]["UCL_warning"]
                if LCL_failed:
                    if float(float(row["value"])) < float(LCL_failed):
                        fastq_id2n_fail[fastq_id] += 1
                        fastq_id2metric2score[fastq_id][row["name"]] = 'FAIL'
                if UCL_failed:
                    if float(float(row["value"])) > float(UCL_failed):
                        fastq_id2n_fail[fastq_id] += 1
                        fastq_id2metric2score[fastq_id][row["name"]] = 'FAIL'
                if row["name"] not in fastq_id2metric2score[fastq_id]:
                    if LCL_warning:
                        if float(float(row["value"])) < float(LCL_warning):
                            fastq_id2n_warn[fastq_id] += 1
                            fastq_id2metric2score[fastq_id][row["name"]] = 'WARN'
                    if UCL_warning:
                        if float(float(row["value"])) > float(UCL_warning):
                            fastq_id2n_warn[fastq_id] += 1
                            fastq_id2metric2score[fastq_id][row["name"]] = 'WARN' 

        return fastq_id2n_fail, fastq_id2n_warn, fastq_id2metric2score

    def get_run_name2run_id(self,):
        sql = 'select run_name,id from GEN_runs'
        return {i[0]:i[1] for i in self.cursor.execute(sql,).fetchall()}

    def match_fastq_to_sample(self, fastq_prefix):

        # mapping based on xlsx sample id
        # conflict possible with several projects with sample with numeric labels (1,2,3, 232,...)
        # use date filter to exclude mapping older than X on the various columns
        sql = f'select id from GEN_sample where xlsx_sample_ID="{fastq_prefix}"' 
        try:
            sample_id = self.cursor.execute(sql,).fetchall()[0][0]
        except:
            sql = f'select id from GEN_sample where sample_name="{fastq_prefix}"' 
            try:
                sample_id = self.cursor.execute(sql,).fetchall()[0][0]
            except:
                sql = f'select id from GEN_sample where alias="{fastq_prefix}"' 
                try:
                    sample_id = self.cursor.execute(sql,).fetchall()[0][0]
                except:
                    sample_id = None 
        return sample_id

    def get_sample_id(self, sample_xls_id):
        sql = f'select id from GEN_sample where xlsx_sample_ID={self.spl}'

        return self.cursor.execute(sql,(sample_xls_id,)).fetchall()[0][0]

    def add_sample_to_fastq_relation(self, fastq_id, sample_id):
        
            # ignore if relation already known
            sql2 = f'insert or ignore into GEN_fastqtosample(fastq_id, sample_id) values({self.spl},{self.spl})'
            
            self.cursor.execute(sql2, 
                               (fastq_id,sample_id))
            
            self.conn.commit()


    def insert_run(self,
                   run_name,
                   run_date,
                   assay,
                   read_length,
                   paired,
                   filearc_folder):
        
        sql = f'''INSERT into GEN_runs (run_name, run_date, assay, read_length, paired, filearc_folder) values({self.spl},{self.spl},{self.spl},{self.spl},{self.spl},{self.spl}) 
            ON CONFLICT(GEN_runs.run_name) DO UPDATE SET run_date={self.spl}, assay={self.spl}, read_length={self.spl}, paired={self.spl}, filearc_folder={self.spl};
        '''
        self.cursor.execute(sql, [run_name,
                                  run_date,
                                  assay,
                                  read_length,
                                  paired, 
                                  filearc_folder,
                                  run_date,
                                  assay,
                                  read_length,
                                  paired, 
                                  filearc_folder])
        self.conn.commit()

    def compare_samples(self, fastq_id_1, fastq_id2,print_data=False):
        
        sql = f'select * from GEN_snps where fastq_id={fastq_id_1}'
        df1 = pandas.read_sql(sql, self.conn)
        sql2 = f'select * from GEN_snps where fastq_id={fastq_id2}'
        df2 = pandas.read_sql(sql2, self.conn)
        df1["VAR"] = df1["ref"] + df1["position"].astype(str) + df1["alt"]
        df2["VAR"] = df2["ref"] + df2["position"].astype(str) + df2["alt"]

        df1_vars = set(df1["VAR"].to_list())
        df2_vars = set(df2["VAR"].to_list())
        
        shared = df1_vars.intersection(df2_vars)
        unique_df1 = df1_vars.difference(df2_vars)
        unique_df2= df2_vars.difference(df1_vars)
        union = df1_vars.union(df2_vars)
        if print_data:
            print("unique_1", df1_vars)
            print("unique_2", df2_vars)

        n_shared = len(shared)
        n_unique_df1 = len(unique_df1)
        n_unique_df2 = len(unique_df2)
        n_union = len(union)
        n_diffs = n_unique_df1 +  n_unique_df2
        return [fastq_id_1, fastq_id2, n_union, n_shared, n_diffs, n_unique_df1, n_unique_df2]

    def parwise_snps_comp(self, fastq_list):
        import itertools
        comb = itertools.combinations(fastq_list, 2)
        res = []
        complete_mat = []
        for one_pair in comb:
            vals = self.compare_samples(one_pair[0], one_pair[1],print_data=False)
            res.append(vals)
            complete_mat.append(vals)
            # reverse comp
            vals = [vals[1], vals[0]] + vals[2:]
            complete_mat.append(vals)
            
        df = pandas.DataFrame(res)
        df.columns = ["fastq1", "fastq2", "n_union", "n_shared", "n_diffs", "n_unique_1", "n_unique_2"]
        df2 = pandas.DataFrame(complete_mat)
        df2.columns = ["fastq1", "fastq2", "n_union", "n_shared", "n_diffs", "n_unique_1", "n_unique_2"]
        m = pandas.pivot(df2, index="fastq1", columns="fastq2", values="n_diffs").fillna(0)

        return df, m


    def insert_sample(self,
                      col_names,
                      values_list,
                      sample_xls_id):
        
        update_str = f'%s=%{self.spl}'

        update_str_comb = ','.join([update_str % colname for colname in col_names])
        # INSERT into GEN_sample(xlsx_sample_ID,species_name,date_received,sample_name,sample_type,analysis_type,description,molis_id,myseq_passage,run_date,date_registered,date_sample_modification,user_creation_id,user_modification_id) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(GEN_sample.xlsx_sample_ID) DO UPDATE SET xlsx_sample_ID=?,species_name=?,date_received=?,sample_name=?,sample_type=?,analysis_type=?,description=?,molis_id=?,myseq_passage=?,run_date=?,date_registered=?,date_sample_modification=?,user_creation_id=?,user_modification_id=?;
        # [29, 'Staphylococcus aureus', False, nan, 'strain', 'research', nan, 1306182530, nan, False, '2020-09-02', '2020-09-02', 2, 2, 29, 'Staphylococcus aureus', False, nan, 'strain', 'research', nan, 1306182530, nan, False, '2020-09-02', '2020-09-02', 2, 2]            
        # update all columns in case of conflict with xlsx_sample_id
        
        # NOTE: xlsx_sample_ID used as reference: if a row is updated in the xlsx table, the corresponding row is updated in the sql table
        sql_template = f'INSERT into GEN_sample(%s) values(%s)' \
                       f' ON CONFLICT(GEN_sample.xlsx_sample_ID) DO UPDATE SET %s;' % (','.join(col_names),
                                                                                      ','.join([f'{self.spl}']*len(col_names)),
                                                                                      update_str_comb)
        
        print(sql_template)
        #print(values_list)
        self.cursor.execute(sql_template, values_list + values_list)
        self.conn.commit()

        return self.get_sample_id(sample_xls_id)
  

    def get_mapped_fastq_list(self,):

        sql = 'select distinct fastq_id from GEN_fastqtosample'
        
        df = pandas.read_sql(sql, self.conn)
        if df.empty:
            return []
        else:
            fastq_id_list = df["fastq_id"].to_list()
            return fastq_id_list

    def match_sample_to_fastq(self, sample_prefix, filter_already_mapped=False):
        sql = 'select id from GEN_fastqfiles where fastq_prefix=?'
        
        try:
            fastq_id_list = [i[0] for i in self.cursor.execute(sql,(sample_prefix,)).fetchall()]
        except:
            return [] 
        if filter_already_mapped:
            return [str(i) for i in fastq_id_list if i not in filter_already_mapped]
        else:
            return fastq_id_list


    def get_fastq(self, run_name=False):
        
        sql = 'select t1.id,run_name,date_run,qc,fastq_prefix,xlsx_sample_ID,species_name,date_received,read_length,t1.id from GEN_fastqfiles t1 ' \
            ' inner join GEN_runs t2 on t1.run_id=t2.id ' \
            ' left join GEN_fastqtosample t3 on t1.id=t3.fastq_id' \
            ' left join GEN_sample t4 on t3.sample_id=t4.id'

        sql = '''
        select t1.id as fastq_id, run_name,t5.status,fastq_prefix,read_length,t3.id as sample_id,taxonomy,date_received,t4.run_date,t1.id,t4.qc_id from GEN_fastqfiles t1 
        left join GEN_fastqtosample t2 on t1.id=t2.fastq_id 
        left join GEN_sample t3 on t2.sample_id=t3.id 
        left join GEN_runs t4 on t1.run_id=t4.id
        left join GEN_analysis t5 on t4.qc_id=t5.id
        '''

        if run_name:
            sql += f'\nwhere run_name="{run_name}"'

        #print(sql)
    
        return pandas.read_sql(sql, self.conn)
      
    
    def get_run_table(self,):
        
        from GEN.models import Analysis
        
        sql = '''select run_date,run_name,read_length,filearc_folder,qc_id,qc_path,count(*) as n_fastq from GEN_runs t1
        inner join GEN_fastqfiles t2 on t1.id=t2.run_id group by run_date,run_name,read_length,filearc_folder,qc_id,qc_path
        ''' 
        data = [list(i) for i in self.cursor.execute(sql,).fetchall()]
        for n, row in enumerate(data):
            if row[4]:
                data[n][4] = Analysis.objects.filter(id=row[4])[0]
        return data

    def get_run_samples(self, run_name):
        
        sql = f'''select fastq_prefix,R1,R2 from GEN_fastqfiles t1 
                  inner join GEN_runs t2 on t1.run_id=t2.id where run_name="{run_name}"
               '''
        return self.cursor.execute(sql,).fetchall()

    def get_run_table_df(self,):
        
        sql = f'''select * from GEN_runs 
               '''
        return pandas.read_sql(sql, self.conn)
    
    
    def get_run_sample2species(self, run_name):
        # left join because some fastq won't have match in the sample table
        sql = f'''select fastq_prefix,species_name from GEN_fastqfiles t1 
                inner join GEN_runs t2 on t1.run_id=t2.id
                left join GEN_fastqtosample t3 on t1.id=t3.fastq_id
                left join GEN_sample t4 on t3.sample_id=t4.id where run_name="{run_name}";
            '''
        return {i[0]: i[1] for i in self.cursor.execute(sql,).fetchall()}
    
        
    def fastq_prefix2fastq_id(self,fastq_prefix_list):
        
        fasta_filter = '","'.join(fastq_prefix_list)
        sql = f'select fastq_prefix,id from GEN_fastqfiles where fastq_prefix in ("{fasta_filter}")'
        
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}
    
    def add_fastq_metadata(self, 
                           fastq_id,
                           term_name,
                           value,
                           analysis_id=False):
        '''
        [{"fastq_id": fastq_id,
        "term_name": <name>,
        "value": <value>,
        "analysis_id" analysis_id}]
        '''
        from GEN.models import Term
        from GEN.models import FastqFilesMetadata
        term = Term.objects.get_or_create(name=term_name)[0]
        if not analysis_id:
            m = FastqFilesMetadata(term=term, fastq_id=fastq_id, value=value)
            m.save()
        else:
            m = FastqFilesMetadata(term=term, fastq_id=fastq_id, value=value, analysis_id=analysis_id)
            m.save()

    def add_sample_metadata(self, 
                           sample_id,
                           term_name,
                           value):
        '''
        TODO: add option to update on conflict
        [{"sample_id": sample_id,
        "term_name": <name>,
        "value": <value>,
        "analysis_id" analysis_id}]
        '''
        from GEN.models import Term
        from GEN.models import SampleMetadata
        term = Term.objects.get_or_create(name=term_name)[0]

        m = SampleMetadata(term=term, sample_id=sample_id, value=value)
        m.save()


    def get_term2term_id(self, 
                         term_list):
        from GEN.models import Term
        term2term_id = {}
        for term_name in term_list:
            term = Term.objects.get_or_create(name=term_name)[0]
            term2term_id[term_name] = term.id
        return term2term_id

    def get_term_id2term_name(self, 
                              term_id_list):
        from GEN.models import Term
        term_id2term_name = {}
        for term_id in term_id_list:
            term = Term.objects.get_or_create(id=term_id)[0]
            term_id2term_name[term_id] = term.name
        return term_id2term_name

    def fastq_mutation(self, aa_change_list):
        
        change_filter = '","'.join(aa_change_list)
        sql = f'select fastq_id,aa_change from GEN_snps where aa_change in ("{change_filter}");'
        df = pandas.read_sql(sql, self.conn)
        fastq2changes = {i:[] for i in df["fastq_id"].to_list()}
        for n, row in df.iterrows():
            fastq2changes[row["fastq_id"]].append(row["aa_change"])

        return fastq2changes


    def add_QC_report(self, run_name, run_path):
        
        sql = 'update GEN_runs set qc=1 where run_name=?'
        self.cursor.execute(sql, [run_name]) 
        
        sql = 'update GEN_runs set qc_path=? where run_name=?'
        self.cursor.execute(sql, (run_path, run_name))
        self.conn.commit()
    
    
    def fastq_qc_filter(self, analysis_id, value):
        
        sql = f'''select fastq_id from GEN_fastqfilesmetadata t1 
                 inner join GEN_term t2 on t1.term_id=t2.id 
                 where t1.analysis_id={analysis_id} and t2.name="qc_status" and t1.value="{value}";
        '''
        print(sql)
        
        return [i[0] for i in self.cursor.execute(sql,)]
        
    def format_snps(self, fastq_id_list):
        
        fastq_filter = ','.join([str(i) for i in fastq_id_list])
        sql = f'select fastq_id,nucl_change, aa_change, gene from GEN_snps where fastq_id in ({fastq_filter})'

        df = pandas.read_sql(sql, self.conn).set_index("fastq_id")
        
        print("head", df.head())
        
        fastq2data = {fastq_id: {} for fastq_id in df.index.unique()}
        
        for fastq_id in list(fastq2data.keys()):
            target_fastq = df.loc[fastq_id]
            aa_changes = [f'{pos["aa_change"]} ({pos["gene"]})' for n, pos in target_fastq.iterrows() if not pandas.isna(pos["aa_change"])]
            nucl_changes = [pos["nucl_change"] for n, pos in target_fastq.iterrows()]
            fastq2data[fastq_id]["aa_changes"] = ';'.join(aa_changes)
            fastq2data[fastq_id]["nucl_changes"] = ';'.join(nucl_changes)
            
        return fastq2data   
        
        
        
    
    def insert_or_get_fastq(self, 
                            fastq_prefix, 
                            run_id, 
                            R1, 
                            R2):
        from GEN.models import FastqFiles

        if "control" in fastq_prefix:
            control = 1
        else:
            control = 0
        # insert 
        fastq = FastqFiles.objects.get_or_create(control_sample=control,
                                                 fastq_prefix=fastq_prefix,
                                                 run_id=run_id,
                                                 R1=R1,
                                                 R2=R2)[0]
        
        return fastq.id

    def get_sample2species(self, sample_list):
        
        sample_list_filter = '","'.join(sample_list)
        sql = f'''select distinct fastq_prefix,species_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where fastq_prefix in ("{sample_list_filter}");
           '''
        
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}
    
    def get_fastq_id2species(self, fastq_list):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_list])
        sql = f'''select distinct fastq_prefix,species_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where t1.id in ("{fastq_list_filter}");
           '''
        print(sql)
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}
    
    def get_fastq_id2sample_name(self, fastq_list, key_str=True):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_list])
        sql = f'''select distinct fastq_id,sample_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where t1.id in ("{fastq_list_filter}");
           '''
        print(sql)
        if key_str:
            return {str(i[0]):i[1] for i in self.cursor.execute(sql,)}
        else:
            return {int(i[0]):i[1] for i in self.cursor.execute(sql,)}

    def get_sample_table(self,):
        
        sql = '''
        select distinct A.id,A.sample_type,A.sample_name,A.taxonomy,A.date_registered, A.date_received,A.n_fastq, count(t3.subproject_id ) as n_projects from (
        select distinct t1.id,t1.sample_type,t1.sample_name,t1.taxonomy,t1.date_registered, t1.date_received, count(t2.fastq_id) as n_fastq from GEN_sample t1 
        left join GEN_fastqtosample t2 on t1.id=t2.sample_id 
        group by t1.id) A  
        left join GEN_subprojectsample t3 on A.id=t3.sample_id
        group by A.id 
        '''

        return self.cursor.execute(sql,).fetchall()
    
    def get_analysis_fastq_list(self,analisis_id):
        
        sql = f'''
        select fastq_id from GEN_fastqset where analysis_id={analisis_id}
        '''
        
        return [i[0] for i in self.cursor.execute(sql,).fetchall()]
        

    def calculate_age(self, birth_date, prel_date):
        from dateutil.relativedelta import relativedelta
        rdelta = relativedelta(prel_date, birth_date)
        age_years = rdelta.years 
        if age_years == 0:
            age_months = rdelta.months
            if age_months == 0:
                age_weeks = rdelta.weeks
                return f'{age_weeks} weeks'
            else:
                return f'{age_months} months'
        else:
            return age_years


    def insert_molis_sample_metadata_from_xml(self, 
                                              df, 
                                              field_definition, 
                                              molis_alias="Numéro alias", 
                                              already_into_db="patient_sex"):
        import datetime
        samples_in_db = [] # set([str(i) for i in self.get_sample_metadata(already_into_db).keys()])
        
        metadata_list = []

        molis2sample_list = self.get_molis_id2sample_id_list()

        for n, row in df.iterrows():
            alias = row[molis_alias]
            try:
                sample_list = molis2sample_list[int(alias)]
            except KeyError:
                continue

            for sample_id in sample_list:
                if str(sample_id) not in samples_in_db:               
                    print(f"missing metadata: sample {sample_id} (alias {alias})")

                    for field in field_definition:
                        if field_definition[field]["type"] == 'simple':
                            val = row[field_definition[field]["fields"][0]]
                        elif field_definition[field]["type"] == 'date':
                            val = row[field_definition[field]["fields"][0]]
                            if val == '':
                                continue 
                            else:
                                val = datetime.datetime.strptime(val, '%Y-%m-%d').strftime('%Y-%m-%d')
                        elif field_definition[field]["type"] == 'split':
                            index = field_definition[field]["fields"][2]
                            delimiter = field_definition[field]["fields"][1]
                            val = row[field_definition[field]["fields"][0]].split(delimiter)[index]
                        elif field_definition[field]["type"] == 'concatenate':
                            field_list = field_definition[field]["fields"]
                            values = row[field_list].to_list()
                            val = ''.join(values)
                        elif field_definition[field]["type"] == 'age':
                            date_a = datetime.datetime.strptime(row[field_definition[field]["fields"][0]], '%Y-%m-%d')
                            date_b = datetime.datetime.strptime(row[field_definition[field]["fields"][1]], '%Y-%m-%d')
                            val = self.calculate_age(date_a, date_b).strftime('%Y-%m-%d')
                        else:
                            field_type = field_definition[field]["type"]
                            print(f'Unknown type: {field_type}')
                            raise IOError(f'Unknown type: {field_type}')
                    
                        # add value if not empty
                        if val != '':
                             metadata_list.append({'term_name':field, 'value': val, 'sample_id': sample_id})

        # insert data
        for metadata in metadata_list:
            # duplicates will raise error...
            print(metadata)
            try:
                self.add_sample_metadata(sample_id=metadata["sample_id"],
                                         term_name=metadata["term_name"],
                                         value=metadata["value"]) 
            except IntegrityError as e:
                if 'UNIQUE' in str(e):
                    print(f'UNIQUE constraint failed: {sample_id} -- {metadata["term_name"]} -- {metadata["value"]}')
                    continue
                else:
                    print("other error")
                    raise



def update_analysis_status(analysis_id, status):
    from GEN.models import Analysis
    m = Analysis.objects.filter(id=analysis_id)[0]
    m.status = status
    m.save()
    
def add_analysis_metadata(analysis_id, term, value, update=False):
    from GEN.models import Term
    from GEN.models import AnalysisMetadata
    term = Term.objects.get_or_create(name=term)[0]
    if not update:
        m = AnalysisMetadata(term=term, analysis_id=analysis_id, value=value)
        m.save()
    else:
        # update value
        m = AnalysisMetadata.objects.filter(term=term, analysis_id=analysis_id)[0]
        m.value = value
        m.save()
        
def create_analysis(fastq_id_list,
                    analysis_description,
                    reference_fastq_id_list = [],
                    subproject_id=False,
                    workflow_name="Airflow_epidemiology"):
                
    from GEN.models import Workflow
    from GEN.models import WorkflowSteps
    from GEN.models import FastqFiles
    from GEN.models import FastqSet
    from GEN.models import Term
    from GEN.models import AnalysisStatus
    from GEN.models import ProjectAnalysis
    from GEN.models import Analysis
    from GEN.models import AnalysisMetadata
    from datetime import datetime
    
    # create new analysis
    workflow_instance = Workflow.objects.filter(workflow_name=workflow_name)
    workflow_instance = workflow_instance[0]
    analysis = Analysis(workflow=workflow_instance, start_date=datetime.today(), status='started')
    analysis.save()
    
    # associated to project if a project id was provided
    if subproject_id:
        project_analysis = ProjectAnalysis(analysis=analysis, subproject_id=subproject_id)
        project_analysis.save()
    
    # insert description as metadata
    term = Term.objects.get_or_create(name="description")[0]
    desc = AnalysisMetadata(term=term, analysis=analysis, value=analysis_description)
    desc.save()
    
    if len(reference_fastq_id_list) > 0:
        term_ref_genome = Term.objects.get_or_create(name="reference_genome")[0]
        for ref_genome in reference_fastq_id_list:
            m = AnalysisMetadata(term=term_ref_genome, analysis=analysis, value=ref_genome)
            m.save()  
    
    # add each fastq to fastq set
    for fastq_id in fastq_id_list:
        print("fastq_id", fastq_id)
        fastq_instance = FastqFiles.objects.filter(id=fastq_id)[0]
        new_set_fastq = FastqSet(analysis=analysis, fastq=fastq_instance)
        new_set_fastq.save()

    # create "AnalysisStatus" entry for each step of the workflow
    # ==> all marked as not DONE
    all_steps = WorkflowSteps.objects.filter(workflow=workflow_instance)
    for one_step in all_steps:
        new_project_analysis_status_instance = AnalysisStatus(analysis=analysis, step=one_step.step, status=0)
        new_project_analysis_status_instance.save()
        
    
    return analysis.id

