


class DB:
    def __init__(self, db_path):
        import sqlite3
        
        self.db_path = db_path
        print("connecting to ", self.db_path)
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        
    def get_fastq_qc_metric(self, metric_name):
        
        sql = f'''select fastq_prefix,metric_value from GEN_fastqfiles t1
                  inner join GEN_runsqc t2 on t1.id=t2.fastq_id
                  inner join GEN_qcmetrics t3 on t2.metric_id=t3.id
                  where metric_name="{metric_name}";'''
                  
        print(sql)

        return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}
    
    def get_fastq(self, run_name=False):
        
        sql = 'select run_name,date_run,qc,fastq_prefix,xlsx_sample_ID,species_name,date_received,read_length from GEN_fastqfiles t1 ' \
            ' inner join GEN_runs t2 on t1.run_id=t2.id ' \
            ' left join GEN_fastqtosample t3 on t1.id=t3.fastq_id' \
            ' left join GEN_sample t4 on t3.sample_id=t4.id'

        sql = '''
        select run_name,qc,fastq_prefix,read_length,t3.id,species_name,date_received,t4.run_date from GEN_fastqfiles t1 
        left join GEN_fastqtosample t2 on t1.id=t2.fastq_id 
        left join GEN_sample t3 on t2.sample_id=t3.id 
        left join GEN_runs t4 on t1.run_id=t4.id
        '''

        if run_name:
            sql += f'\nwhere run_name="{run_name}"'

        print(sql)
    
        return self.cursor.execute(sql,).fetchall()
    
    def get_run_table(self,):
        sql = '''select run_date,run_name,read_length,filearc_folder,qc,qc_path,count(*) as n_fastq from GEN_runs t1
        inner join GEN_fastqfiles t2 on t1.id=t2.run_id group by run_date,run_name,read_length,filearc_folder,qc,qc_path
        '''       
        return [list(i) for i in self.cursor.execute(sql,).fetchall()]

    def get_run_samples(self, run_name):
        
        sql = f'''select fastq_prefix,R1,R2 from GEN_fastqfiles t1 
                  inner join GEN_runs t2 on t1.run_id=t2.id where run_name="{run_name}"
               '''
        return self.cursor.execute(sql,).fetchall()
    
    
    def get_run_sample2species(self, run_name):
        # left join because some fastq won't have match in the sample table
        sql = f'''select fastq_prefix,species_name from GEN_fastqfiles t1 
                inner join GEN_runs t2 on t1.run_id=t2.id
                left join GEN_fastqtosample t3 on t1.id=t3.fastq_id
                left join GEN_sample t4 on t3.sample_id=t4.id where run_name="{run_name}";
            '''
        return {i[0]: i[1] for i in self.cursor.execute(sql,).fetchall()}
    
    def add_metrics(self, metrics_name):
        
        sql = f'insert into GEN_qcmetrics(metric_name) values(?)'
        self.cursor.execute(sql, [metrics_name])
        self.conn.commit()
        
        return self.cursor.lastrowid
        
    
    def get_metrics_name2metrics_id(self, metrics_list):
        
        metrics_filter = '","'.join(metrics_list)
        sql = f'select metric_name,id from GEN_qcmetrics where metric_name in ("{metrics_filter}")'   
        dico = {i[0]:i[1] for i in self.cursor.execute(sql,)}
                
        for metric in metrics_list:
            if metric not in dico:
                dico[metric] = self.add_metrics(metric)
        
        return dico
        
    def fastq_prefix2fastq_id(self,fastq_prefix_list):
        
        fasta_filter = '","'.join(fastq_prefix_list)
        sql = f'select fastq_prefix,id from GEN_fastqfiles where fastq_prefix in ("{fasta_filter}")'
        
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}
    
    def add_metrics_values(self, lst):
        '''
        [{"sample": sample,
        "metrics_name": metric.lower(),
        "metrics_value": metric_value,
        "pipeline_version" version}]
        '''
        
        sample_name_list = list(set([i["sample"] for i in lst]))
        metric_list = list(set([i["metrics_name"] for i in lst]))
        
        metric2id = self.get_metrics_name2metrics_id(metric_list)
        sample2id = self.fastq_prefix2fastq_id(sample_name_list)
        
        sql_insert = 'INSERT OR REPLACE into GEN_runsqc (fastq_id,pipeline_version,metric_id,metric_value) values (?,?,?,?)'
       
        for entry in lst:
            self.cursor.execute(sql_insert,[
                                            sample2id[entry["sample"]],
                                            entry["pipeline_version"],
                                            metric2id[entry["metrics_name"]],
                                            entry["metrics_value"],
                                           ])
        self.conn.commit()
        
            
    def add_QC_report(self, run_name, run_path):
        
        sql = 'update GEN_runs set qc=1 where run_name=?'
        self.cursor.execute(sql, [run_name]) 
        
        sql = 'update GEN_runs set qc_path=? where run_name=?'
        self.cursor.execute(sql, (run_path, run_name))
        self.conn.commit()
    
    
    def get_sample2species(self, sample_list):
        
        sample_list_filter = '","'.join(sample_list)
        sql = f'''select distinct fastq_prefix,species_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where fastq_prefix in ("{sample_list_filter}");
           '''
    
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}
    
    
    def get_sample2species(self, sample_list):
        
        sample_list_filter = '","'.join(sample_list)
        sql = f'''select distinct fastq_prefix,species_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where fastq_prefix in ("{sample_list_filter}");
           '''
    
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}