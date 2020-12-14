import pandas
from django.conf import settings
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
    def __init__(self, db_path):
        import sqlite3
        
        self.db_path = db_path
        print("connecting to ", self.db_path)
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        
    def get_fastq_metadata(self, metric_name):
        
        sql = f'''select t1.id,value from GEN_fastqfiles t1
                  inner join GEN_fastqfilesmetadata t2 on t1.id=t2.fastq_id
                  inner join GEN_term t3 on t2.term_id=t3.id
                  where t3.name="{metric_name}";'''
                  
        print(sql)

        return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}

    
    def get_fastq_and_sample_data(self, fastq_id_list):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_id_list])
        
        # left join because some fastq won't have match in the sample table
        # possible problem: fastq prefix match with multiple samples from different species
        # in that case: remove species name
        sql = f'''select distinct t1.id as fastq_id,fastq_prefix,R1,R2,species_name from GEN_fastqfiles t1 
                left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
                left join GEN_sample t3 on t2.sample_id=t3.id 
                where t1.id in ("{fastq_list_filter}");
            '''
        
        return pandas.read_sql(sql, self.conn)


    def get_analysis_metadata(self, term_name):
        
        sql = f"""select t1.analysis_id,value from GEN_analysismetadata t1 
                  inner join GEN_term t2 on t1.term_id =t2.id 
                  where t2.name like '{term_name}';"""

        return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}


    def get_run_name2run_id(self,):
        sql = 'select run_name,id from GEN_runs'
        return {i[0]:i[1] for i in self.cursor.execute(sql,).fetchall()}

    def match_fastq_to_sample(fastq_prefix):
        sql = f'select id from GEN_sample where sample_name="{fastq_prefix}"' 
        try:
            sample_id = self.cursor.execute(sql,).fetchall()[0][0]
        except:
            sample_id = None 
        return sample_id

    def get_sample_id(self, sample_xls_id):
        sql = 'select id from GEN_sample where xlsx_sample_ID=?'

        return self.cursor.execute(sql,(sample_xls_id,)).fetchall()[0][0]


    def add_sample_to_fastq_relation(self, fastq_id, sample_id):
        
            # ignore if relation already known
            sql2 = 'insert or ignore into GEN_fastqtosample(fastq_id, sample_id) values(?,?)'
            
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
        
        sql = '''INSERT into GEN_runs (run_name, run_date, assay, read_length, paired, filearc_folder) values(?,?,?,?,?,?) 
            ON CONFLICT(GEN_runs.run_name) DO UPDATE SET run_date=?, assay=?, read_length=?, paired=?, filearc_folder=?;
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

    def insert_sample(self,
                      col_names,
                      values_list,
                      sample_xls_id):
        
        update_str = '%s=?'

        update_str_comb = ','.join([update_str % colname for colname in col_names])
        # INSERT into GEN_sample(xlsx_sample_ID,species_name,date_received,sample_name,sample_type,analysis_type,description,molis_id,myseq_passage,run_date,date_registered,date_sample_modification,user_creation_id,user_modification_id) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(GEN_sample.xlsx_sample_ID) DO UPDATE SET xlsx_sample_ID=?,species_name=?,date_received=?,sample_name=?,sample_type=?,analysis_type=?,description=?,molis_id=?,myseq_passage=?,run_date=?,date_registered=?,date_sample_modification=?,user_creation_id=?,user_modification_id=?;
        # [29, 'Staphylococcus aureus', False, nan, 'strain', 'research', nan, 1306182530, nan, False, '2020-09-02', '2020-09-02', 2, 2, 29, 'Staphylococcus aureus', False, nan, 'strain', 'research', nan, 1306182530, nan, False, '2020-09-02', '2020-09-02', 2, 2]            
        # update all columns in case of conflict with xlsx_sample_id
        
        # NOTE: xlsx_sample_ID used as reference: if a row is updated in the xlsx table, the corresponding row is updated in the sql table
        sql_template = 'INSERT into GEN_sample(%s) values(%s)' \
                       ' ON CONFLICT(GEN_sample.xlsx_sample_ID) DO UPDATE SET %s;' % (','.join(col_names),
                                                                                      ','.join(['?']*len(col_names)),
                                                                                      update_str_comb)
               
        self.cursor.execute(sql_template, values_list + values_list)
        self.conn.commit()

        return self.get_sample_id(sample_xls_id)
  

    def match_sample_to_fastq(self, sample_prefix):
        sql = 'select id from GEN_fastqfiles where fastq_prefix=?' 
        try:
            fastq_id = self.cursor.execute(sql,(sample_prefix,)).fetchall()[0][0]
        except:
            fastq_id = None 
        return fastq_id


    def get_fastq(self, run_name=False):
        
        sql = 'select run_name,date_run,qc,fastq_prefix,xlsx_sample_ID,species_name,date_received,read_length,t1.id from GEN_fastqfiles t1 ' \
            ' inner join GEN_runs t2 on t1.run_id=t2.id ' \
            ' left join GEN_fastqtosample t3 on t1.id=t3.fastq_id' \
            ' left join GEN_sample t4 on t3.sample_id=t4.id'

        sql = '''
        select run_name,t5.status,fastq_prefix,read_length,t3.id,species_name,date_received,t4.run_date,t1.id from GEN_fastqfiles t1 
        left join GEN_fastqtosample t2 on t1.id=t2.fastq_id 
        left join GEN_sample t3 on t2.sample_id=t3.id 
        left join GEN_runs t4 on t1.run_id=t4.id
        left join GEN_analysis t5 on t4.qc_id=t5.id
        '''

        if run_name:
            sql += f'\nwhere run_name="{run_name}"'

        print(sql)
    
        return self.cursor.execute(sql,).fetchall()
      
    
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
            
    def add_QC_report(self, run_name, run_path):
        
        sql = 'update GEN_runs set qc=1 where run_name=?'
        self.cursor.execute(sql, [run_name]) 
        
        sql = 'update GEN_runs set qc_path=? where run_name=?'
        self.cursor.execute(sql, (run_path, run_name))
        self.conn.commit()
    
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
    
    def get_sample_table(self,):
        
        sql = '''
        select distinct A.id,A.sample_type,A.sample_name,A.species_name,A.date_registered, A.date_received,A.n_fastq, count(t3.subproject_id ) as n_projects from (
        select distinct t1.id,t1.sample_type,t1.sample_name,t1.species_name,t1.date_registered, t1.date_received, count(t2.fastq_id) as n_fastq from GEN_sample t1 
        left join GEN_fastqtosample t2 on t1.id=t2.sample_id 
        group by t1.id) A  
        left join GEN_subprojectsample t3 on A.id=t3.sample_id
        group by A.id 
        '''

        return self.cursor.execute(sql,).fetchall()
    

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
