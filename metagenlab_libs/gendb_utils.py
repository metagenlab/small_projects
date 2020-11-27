import pandas
from django.conf import settings
try:
    import GEN_database.settings as GEN_settings

    settings.configure(INSTALLED_APPS=GEN_settings.INSTALLED_APPS,
                    DATABASES=GEN_settings.DATABASES)

    import django
    django.setup()
except:
    pass

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
        
        sql = 'select run_name,date_run,qc,fastq_prefix,xlsx_sample_ID,species_name,date_received,read_length,t1.id from GEN_fastqfiles t1 ' \
            ' inner join GEN_runs t2 on t1.run_id=t2.id ' \
            ' left join GEN_fastqtosample t3 on t1.id=t3.fastq_id' \
            ' left join GEN_sample t4 on t3.sample_id=t4.id'

        sql = '''
        select run_name,qc,fastq_prefix,read_length,t3.id,species_name,date_received,t4.run_date,t1.id from GEN_fastqfiles t1 
        left join GEN_fastqtosample t2 on t1.id=t2.fastq_id 
        left join GEN_sample t3 on t2.sample_id=t3.id 
        left join GEN_runs t4 on t1.run_id=t4.id
        '''

        if run_name:
            sql += f'\nwhere run_name="{run_name}"'

        print(sql)
    
        return self.cursor.execute(sql,).fetchall()
    
    def get_fastq_metadata(self, fastq_id_list):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_id_list])
        # left join because some fastq won't have match in the sample table
        sql = f'''select t1.id as fastq_id,fastq_prefix,R1,R2,species_name from GEN_fastqfiles t1 
                left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
                left join GEN_sample t3 on t2.sample_id=t3.id 
                where t1.id in ("{fastq_list_filter}");
            '''
        
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
    
    def get_fastq_id2species(self, fastq_list):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_list])
        sql = f'''select distinct fastq_prefix,species_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where t1.id in ("{fastq_list_filter}");
           '''
        print(sql)
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}
    

def update_analysis_status(analysis_id, status):
    from GEN.models import Analysis
    m = Analysis.objects.filter(id=analysis_id)[0]
    m.status = status
    m.save()
    
def add_analysis_metadata(analysis, term, value, update=False):
    from GEN.models import Term
    from GEN.models import AnalysisMetadata
    if not update:
        term = Term.objects.get_or_create(name=term)[0]
        m = AnalysisMetadata(term=term, analysis_id=analysis, value=value)
        m.save()
    else:
        m = AnalysisMetadata.objects.filter(term=term_status, analysis_id=analysis, value=value)[0]
        m.value = status
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
