import pandas
from django.conf import settings
from django.db import IntegrityError
import os 
import yaml 
import networkx
import numpy

# setup django do be able to access django db models 
import GEN_database.settings as GEN_settings

print("DB_DRIVER", GEN_settings.DB_DRIVER)

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
        
        
        self.db_type = GEN_settings.DB_DRIVER

        self.AIRFLOW_CONFIG =  yaml.safe_load(open(GEN_settings.AIRFLOW_CONF, 'r'))

        if self.db_type != "sqlite":
            '''
            import MySQLdb
            sqlpsw = os.environ['SQLPSW']
            self.conn = MySQLdb.connect(passwd=sqlpsw,
                                        user="root",
                                        host="127.0.0.1",
                                        db="GEN_LIMS",
                                        charset='utf8')
            self.cursor = self.conn.cursor()
            '''
            sqlpsw = os.environ['SQLPSW']
            from sqlalchemy import create_engine
            
            engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/GEN_LIMS")
            self.engine_conn = engine.connect()
            self.conn = engine.raw_connection()
            self.cursor = self.conn.cursor()

            # placeholder for sql querries (differ between sqlite and mysql)
            self.spl = '%s'
        else:
            import sqlite3
            self.db_path = GEN_settings.SQLITE_DB
            self.conn = sqlite3.connect(self.db_path)
            self.engine_conn = self.conn
            self.cursor = self.conn.cursor()
            # placeholder for sql querries (differ between sqlite and mysql)
            self.spl = '?'

    def close_cursor(self,):
        self.cursor.close() 
        self.cursor = self.conn.cursor()    


    def get_fastq_metadata(self, 
                           metric_name, 
                           index_str=True, 
                           analysis_id=False):
        
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

        self.cursor.execute(sql,)
        if index_str:
            return {str(i[0]):i[1] for i in self.cursor.fetchall()}
        else:
            return {int(i[0]):i[1] for i in self.cursor.fetchall()}

    def get_fastq_metadata_v2(self, 
                              metric_name,
                              fastq_id_list=False,
                              analysis_id_list=False):
        fastq_filter = ''
        sample_filter = ''
        if fastq_id_list:
            fastq_filter_str = '","'.join([str(i) for i in fastq_id_list])
            fastq_filter += f'and fastq_id in ("{fastq_filter_str}") '
            sample_filter += f'and fastq_id in ("{fastq_filter_str}") '
        if analysis_id_list:
            filter_str = ','.join([str(i) for i in analysis_id_list])
            fastq_filter += f' and analysis_id in ({filter_str})'
        else:
            analysis_filter = ''
        sql = f'''select analysis_id,fastq_id,value from  GEN_fastqfilesmetadata t1
                  inner join GEN_term t2 on t1.term_id=t2.id
                  where t2.name="{metric_name}" {fastq_filter}       
                  union
                  select NULL as analysis_id, t3.fastq_id,t1.value from  GEN_samplemetadata t1
                  inner join GEN_term t2 on t1.term_id=t2.id
                  inner join GEN_fastqtosample t3 on t1.sample_id=t3.sample_id
                  where t2.name="{metric_name}" {sample_filter}
                  '''
        print(sql)
        return pandas.read_sql(sql, self.conn)
           
        

    def get_sample_metadata(self, metric_name, index_str=True):
        sql = f'''select t1.sample_id,t1.value from GEN_samplemetadata t1
                  inner join GEN_term t2 on t1.term_id=t2.id
                  where t2.name="{metric_name}";'''
        self.cursor.execute(sql,)   
        if index_str:
            return {str(i[0]):i[1] for i in self.cursor.fetchall()}
        else:
            return {int(i[0]):i[1] for i in self.cursor.fetchall()}

    def get_fastq_id2run_name(self,):

        sql = '''select t1.id,run_name from GEN_fastqfiles t1
                 inner join GEN_runs t2 on t1.run_id=t2.id
        '''
        self.cursor.execute(sql,)
        return {str(i[0]):i[1] for i in self.cursor.fetchall()}

    def get_metadata_labels(self,):

        sql = '''select distinct t2.id,t2.name from GEN_fastqfilesmetadata t1
                 inner join GEN_term t2 on t1.term_id=t2.id
                 union 
                 select distinct t2.id,t2.name from GEN_samplemetadata t1
                 inner join GEN_term t2 on t1.term_id=t2.id
        '''
        self.cursor.execute(sql,) 
        return {str(i[0]):i[1] for i in self.cursor.fetchall()}


    def get_molis_id2fastq_id(self,):

        sql = ''' 
        select fastq_id, molis_id from GEN_fastqtosample t1
        inner join GEN_sample t2 on t1.sample_id=t2.id
        '''
        self.cursor.execute(sql,) 
        return {str(i[0]):i[1] for i in self.cursor.fetchall()}

    def get_molis_id2sample_id_list(self,):

        sql = ''' 
        select molis_id, id from GEN_sample
        '''

        molis_id2sample_list = {}
        self.cursor.execute(sql,) 
        for row in self.cursor.fetchall():
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

    def get_workflow_config(self, workflow_id):
        from GEN.models import Workflow
        
        workflow_name = Workflow.objects.filter(id=workflow_id)[0].workflow_name
        base = self.AIRFLOW_CONFIG["WORKFLOW"][workflow_name]["BASE_METRICS"]
        add = self.AIRFLOW_CONFIG["WORKFLOW"][workflow_name]["ADDITIONAL_METRICS"]
        color = self.AIRFLOW_CONFIG["WORKFLOW"][workflow_name]["APP_COLOR_FACTOR"]
        # APP_COLOR_FACTOR
        # ADDITIONAL_METRICS
        # BASE_METRICS
        return base, add, color

    def get_workflows(self,):
        from GEN.models import Workflow
        return {i.id:i.workflow_name for i in Workflow.objects.all() }


    def get_fastq_metadata_list(self, 
                                term_list=False, 
                                fastq_filter=None, 
                                run_name_list=False,
                                metadata_value_list=False,
                                add_molis=False,
                                analysis_id_list=False,
                                subproject_id_list=False,
                                workflow_id_list=False,
                                exclude_analysis_id_list=False,
                                fastq_intersection=False):
        '''
        retrieve metadata from both sample and fastq metadata table
        TODO will curently combine dta from multiple anlayses => should retun analysis_id as well
        '''

        print("add molis:", add_molis)

        res_filter_fastq = ''
        res_filter_sample = ''
        workflow_filter = ''
        if term_list:
            print("term_list", term_list)
            term_filter = '","'.join(term_list)
            res_filter_fastq += f' and t2.name in ("{term_filter}")\n' 
            res_filter_sample += f' and t2.name in ("{term_filter}")\n' 
        if fastq_filter:
            fastq_filter_str = ','.join([str(i) for i in fastq_filter])
            res_filter_fastq += f' and t1.fastq_id in ({fastq_filter_str})\n'
            res_filter_sample += f' and t3.fastq_id in ({fastq_filter_str})\n' 
        if run_name_list:
            run_filter = '","'.join(run_name_list)
            res_filter_fastq += f' and run_name in ("{run_filter}")'
            res_filter_sample += f' and run_name in ("{run_filter}")'
        if metadata_value_list:
            metadata_filter = '","'.join(metadata_value_list)
            res_filter_fastq += f' and t1.value in ("{metadata_filter}")'
            res_filter_sample += f' and t1.value in ("{metadata_filter}")'
        if analysis_id_list:
            metadata_filter = '","'.join(analysis_id_list)
            res_filter_fastq += f' and t1.analysis_id in ("{metadata_filter}")'
            res_filter_sample += f' and t6.analysis_id in ("{metadata_filter}")'    
        if subproject_id_list:
            metadata_filter = '","'.join(subproject_id_list)
            res_filter_fastq += f' and t5.subproject_id in ("{metadata_filter}")'  
        if workflow_id_list:
            workflow_filter = 'inner join GEN_analysis t6 on t1.analysis_id=t6.id'
            metadata_filter = '","'.join([str(i) for i in workflow_id_list])
            res_filter_fastq += f' and t6.workflow_id in ({metadata_filter})'
        if exclude_analysis_id_list:
            metadata_filter = ','.join(exclude_analysis_id_list)
            res_filter_fastq += f' and t1.analysis_id not in ({metadata_filter})'
            res_filter_sample += f' and t6.analysis_id not in ({metadata_filter})' 
        
        sql_fastq = f'''
            select distinct t1.fastq_id, t2.name,t1.value, run_name from GEN_fastqfilesmetadata t1 
            inner join GEN_term t2 on t1.term_id=t2.id 
            inner join GEN_fastqfiles t3 on t1.fastq_id=t3.id 
            inner join GEN_runs t4 on t3.run_id=t4.id 
            left join GEN_projectanalysis t5 on t1.analysis_id=t5.analysis_id
            {workflow_filter}
            where t3.fastq_prefix not like "Undetermined%"
            {res_filter_fastq}
        '''
        sql_sample = f'''
            select distinct t3.fastq_id, t2.name,t1.value,run_name from GEN_samplemetadata t1 
            inner join GEN_term t2 on t1.term_id=t2.id 
            inner join GEN_fastqtosample t3 on t1.sample_id=t3.sample_id
            inner join GEN_fastqfiles t4 on t3.fastq_id=t4.id 
            inner join GEN_runs t5 on t4.run_id=t5.id 
            left join GEN_fastqset t6 on t3.fastq_id=t6.fastq_id
            where t4.fastq_prefix not like "Undetermined%"
            {res_filter_sample}
            '''
            
        df_fastq = pandas.read_sql(sql_fastq, self.conn)
        df_samples = pandas.read_sql(sql_sample, self.conn).set_index("fastq_id")
        print("HEAD", df_samples.head())
        
        if fastq_intersection:
            intersection = set(df_fastq["fastq_id"].to_list()).intersection(set(df_samples.index.to_list()))
            df_samples = df_samples.loc[intersection,]

        df = pandas.concat([df_fastq, df_samples.reset_index()])
        
        if add_molis:
            df_molis = self.get_fastq_and_sample_data(df["fastq_id"].to_list()).set_index("fastq_id")
            df = df.set_index("fastq_id").join(df_molis, on="fastq_id", rsuffix='_other', how="left")
            df = df[["molis_id", "name", "value", "run_name"]]

        return df

    def get_snp_matrix(self, analysis_id): 

        sql = f'select fastq1_id, fastq2_id,n_diffs from GEN_snp_distance where analysis_id={analysis_id} union select fastq2_id as fastq1_id,fastq1_id as fastq2_id,n_diffs  from GEN_snp_distance where analysis_id={analysis_id}'
        snp_long = pandas.read_sql(sql, self.conn)
        snp_matrix = pandas.pivot(snp_long, index="fastq1_id", columns="fastq2_id", values="n_diffs").fillna(0)

        return snp_matrix



    def plot_phylogeny_from_parwise_snps_with_context(self, analysis_id):

        # retrieve data 

        backup_folder = self.AIRFLOW_CONFIG["OUTPUT_FOLDER_BASE"]
    
    
        # retrieve snp table and metadata table location 
        pairsnp_filtered_path = os.path.join(backup_folder, self.get_analysis_metadata("pairsnp_filtered",analysis_id_list=[analysis_id])[analysis_id])
        pairsnp_metadata_path = os.path.join(backup_folder, self.get_analysis_metadata("pairsnp_metadata",analysis_id_list=[analysis_id])[analysis_id])

        print("path", pairsnp_filtered_path)
        print("path", pairsnp_metadata_path)

        pairsnp_filtered_df = pandas.read_csv(pairsnp_filtered_path, sep="\t", header=0)[["genome_1", "genome_2", "SNPs"]]
        pairsnp_filtered_df_reverse = pandas.read_csv(pairsnp_filtered_path, sep="\t", header=0)[["genome_2", "genome_1", "SNPs"]]
        pairsnp_filtered_df_reverse.columns = ["genome_1", "genome_2", "SNPs"]
        pairsnp_filtered_df = pairsnp_filtered_df.append(pairsnp_filtered_df_reverse)
        pairsnp_metadata_df = pandas.read_csv(pairsnp_metadata_path, sep="\t", header=0, index_col="strain")

        print(pairsnp_filtered_df.head())

        nr_set = set(pairsnp_filtered_df["genome_1"].to_list() + pairsnp_filtered_df["genome_2"].to_list())
        
        print("nr_set", len(nr_set))

        df_self_comp = pandas.DataFrame([[i, i, 0] for i in nr_set], columns=["genome_1", "genome_2", "SNPs"])

        pairsnp_filtered_df = pairsnp_filtered_df.append(df_self_comp)

        print(pairsnp_filtered_df.tail())

        snp_matrix = pandas.pivot(pairsnp_filtered_df, index="genome_1", columns="genome_2", values="SNPs").fillna(20)

        snp_matrix.to_csv("/media/IMU/GEN/PROJECTS/156_NOSOCOV/Analysis/phylogeny_test_LIMS_with_context_v2.csv", sep="\t")

        cols = set(snp_matrix.columns)
        rows = set(snp_matrix.index)

        print("unique cols",len(cols),cols.difference(rows) )
        print("unique rows",len(rows),rows.difference(cols) )

        print("snp_matrix", snp_matrix.shape)
        print(snp_matrix)

        # nj
        from skbio import DistanceMatrix
        from skbio.tree import nj
        
        print(snp_matrix.shape)
        print()
        dist_mat = snp_matrix.astype('int64').values.tolist()

        dm = DistanceMatrix(dist_mat, [str(i) for i in snp_matrix.index])

        newick_str = nj(dm, result_constructor=str)

        from metagenlab_libs import ete_phylo
        ete_tree = ete_phylo.EteTool(newick_str)

        '''
        # clustering 
        from scipy.cluster import hierarchy

        Z = hierarchy.linkage(snp_matrix, 'single')
        tree = hierarchy.to_tree(Z,False)
        
        # convert to newick
        from metagenlab_libs import ete_phylo

        newick = ete_phylo.get_newick(tree, "", tree.dist, snp_matrix.index)
        ete_tree = ete_phylo.EteTool(newick)
        '''
        # generate plot
        ################

        # retrieve metadata 
        fastq_id2lineage = pairsnp_metadata_df["pango_lineage"].to_dict()
        
        fastq_id2sample_date = pairsnp_metadata_df["date"].to_dict()
        fastq_id2name = {key:value for key,value in pairsnp_metadata_df["patient_name"].to_dict().items() if not pandas.isna(value)}
        fastq_id2IPP = {key:int(value) for key,value in pairsnp_metadata_df["ipp"].to_dict().items() if not pandas.isna(value)}
        fastq_id2qc_status = {key:value for key,value in pairsnp_metadata_df["qc_status"].to_dict().items() if not pandas.isna(value)}
        try:
            fastq_id2Cluster = {key:value for key,value in pairsnp_metadata_df["Cluster"].to_dict().items() if not pandas.isna(value)}
            fastq_id2Service_test = {key:value for key,value in pairsnp_metadata_df["Service_test"].to_dict().items() if not pandas.isna(value)}
            fastq_idCode_Service = {key:value for key,value in pairsnp_metadata_df["Code_Service"].to_dict().items() if not pandas.isna(value)}
        except:
            pass
        fastq_id2origin = {key:("CHUV" if value=='local' else "Switzerland") for key,value in pairsnp_metadata_df["origin"].to_dict().items()}
        
        # plot
        taxon2new_taxon = {i:i for i in snp_matrix.index}

        fastq_id2year = {i:fastq_id2sample_date[i].split("-")[0] for i in fastq_id2sample_date}
        fastq_id2month = {i:fastq_id2sample_date[i].split("-")[1] for i in fastq_id2sample_date}

        # 
        fastq_id2week = {i:pandas.to_datetime(fastq_id2sample_date[i]).week for i in fastq_id2sample_date}

        fastq_id2fastq_id = {i:i for i in fastq_id2sample_date}

        ete_tree.add_text_face(fastq_id2origin, 
                            header_name="Origin",
                            color_scale=True)

        ete_tree.add_text_face(fastq_id2IPP, 
                            header_name="IPP",
                            color_scale=False)

        ete_tree.add_text_face(fastq_id2fastq_id, 
                            header_name="fastq_id",
                            color_scale=False)
        try:
            ete_tree.add_text_face(fastq_id2qc_status, 
                                header_name="QC",
                                color_scale=True)
        except:
            pass

        ete_tree.add_text_face(fastq_id2lineage, 
                            header_name="lineage",
                            color_scale=True)
        try:
            ete_tree.add_text_face(fastq_id2Cluster, 
                                header_name="Cluster",
                                color_scale=True)

            # Service_test	Service_acquisition	Code_Service
            ete_tree.add_text_face(fastq_id2Service_test, 
                                header_name="Service_test",
                                color_scale=True)

            ete_tree.add_text_face(fastq_idCode_Service, 
                                header_name="Code_Service",
                                color_scale=True)
        except:
            pass

        ete_tree.add_text_face(fastq_id2sample_date, 
                            header_name="date",
                            color_scale=False)

        ete_tree.add_text_face(fastq_id2week, 
                            header_name="week",
                            color_scale=True)

        #ete_tree.add_text_face(fastq_id2month, 
        #                    header_name="month",
        #                    color_scale=True)
        

        ete_tree.add_text_face(fastq_id2name, 
                            header_name="patient",
                            color_scale=False)
        
        '''
        # ADD GROUPING DATA
        for group in group_list:
            fastq_id2group = {}
            group_data, group_label = group
            for i, group in enumerate(group_data):
                for member in group:
                    fastq_id2group[str(member)] = i
                    
            # ADD TO PLOT
            ete_tree.add_text_face(fastq_id2group, 
                                   header_name=group_label,
                                   color_scale=True)
        '''

        ete_tree.rename_leaves(taxon2new_taxon,
                               keep_original=False,
                               add_face=False)
        

        print("plotting...")
        ete_tree.remove_dots()
        os.environ['QT_QPA_PLATFORM']='offscreen'
        ete_tree.tree.render("/media/IMU/GEN/PROJECTS/156_NOSOCOV/Analysis/phylogeny_test_LIMS_external_v2.svg",tree_style=ete_tree.tss, w=183, units="mm")



    def plot_phylogeny_from_parwise_snps(self, analysis_id, group_list):

        # retrieve data 

        snp_matrix = self.get_snp_matrix(analysis_id)

        #snp_matrix.to_csv("/media/IMU/GEN/PROJECTS/135_COVIDGEN_20_019/NOSOCOV/results/phylogeny_test_LIMS.csv", sep="\t")

        # nj
        from skbio import DistanceMatrix
        from skbio.tree import nj
        
        print(snp_matrix.shape)
        print()
        dist_mat = snp_matrix.astype('int64').values.tolist()

        dm = DistanceMatrix(dist_mat, [str(i) for i in snp_matrix.index])

        newick_str = nj(dm, result_constructor=str)

        from metagenlab_libs import ete_phylo
        ete_tree = ete_phylo.EteTool(newick_str)

        '''
        # clustering 
        from scipy.cluster import hierarchy

        Z = hierarchy.linkage(snp_matrix, 'single')
        tree = hierarchy.to_tree(Z,False)
        
        # convert to newick
        from metagenlab_libs import ete_phylo

        newick = ete_phylo.get_newick(tree, "", tree.dist, snp_matrix.index)
        ete_tree = ete_phylo.EteTool(newick)
        '''
        # generate plot
        ################

        # retrieve metadata 
        fastq_id2lineage = self.get_fastq_metadata("pangolin_lineage", index_str=True)
        fastq_id2sample_date = self.get_fastq_metadata("sample_date", index_str=True)
        fastq_id2first_name = self.get_fastq_metadata("first_name", index_str=True)
        fastq_id2last_name = self.get_fastq_metadata("last_name", index_str=True)
        fastq_id2IPP = self.get_fastq_metadata("patient_id", index_str=True)
        fastq_id2unit = self.get_fastq_metadata("NOSOCOV_unit", index_str=True)
        fastq_id2n_SNPs = self.get_fastq_metadata("freebayes_n_SNPs", index_str=True)
        fastq_id2qc_status = self.get_fastq_metadata("qc_status", index_str=True)
        
        fastq_id2name = {str(i):f"{fastq_id2first_name[str(i)][0]}. {fastq_id2last_name[str(i)]}" for i in snp_matrix.index}
        # plot
        taxon2new_taxon = {i:i for i in snp_matrix.index}

        fastq_id2year = {i:fastq_id2sample_date[i].split("-")[0] for i in fastq_id2sample_date}
        fastq_id2month = {i:fastq_id2sample_date[i].split("-")[1] for i in fastq_id2sample_date}

        fastq_id2fastq_id = {i:i for i in fastq_id2sample_date}

        ete_tree.add_text_face(fastq_id2IPP, 
                            header_name="IPP",
                            color_scale=False)

        ete_tree.add_text_face(fastq_id2fastq_id, 
                            header_name="fastq_id",
                            color_scale=False)

        ete_tree.add_text_face(fastq_id2qc_status, 
                            header_name="QC",
                            color_scale=True)

        ete_tree.add_text_face(fastq_id2lineage, 
                            header_name="lineage",
                            color_scale=True)

        ete_tree.add_text_face(fastq_id2n_SNPs, 
                            header_name="SNPs",
                            color_scale=False)

        ete_tree.add_text_face(fastq_id2unit, 
                            header_name="unit",
                            color_scale=True)

        ete_tree.add_text_face(fastq_id2sample_date, 
                            header_name="date",
                            color_scale=False)

        ete_tree.add_text_face(fastq_id2year, 
                            header_name="year",
                            color_scale=True)

        ete_tree.add_text_face(fastq_id2month, 
                            header_name="month",
                            color_scale=True)

        '''
        ete_tree.add_text_face(ipp2group, 
                            header_name="cluster",
                            color_scale=True)

        ete_tree.add_text_face(ipp2group_1, 
                            header_name="cluster1",
                            color_scale=True)
                            
        ete_tree.add_text_face(ipp2group_3, 
                            header_name="cluster2",
                            color_scale=True)
        '''

        ete_tree.add_text_face(fastq_id2name, 
                            header_name="patient",
                            color_scale=False)
        
        # ADD GROUPING DATA
        for group in group_list:
            fastq_id2group = {}
            group_data, group_label = group
            for i, group in enumerate(group_data):
                for member in group:
                    fastq_id2group[str(member)] = i
                    
            # ADD TO PLOT
            ete_tree.add_text_face(fastq_id2group, 
                                   header_name=group_label,
                                   color_scale=True)
                 
        ete_tree.rename_leaves(taxon2new_taxon,
                               keep_original=False,
                               add_face=False)
        

        print("plotting...")
        ete_tree.remove_dots()

        #ete_tree.tree.render("/media/IMU/GEN/PROJECTS/135_COVIDGEN_20_019/NOSOCOV/results/phylogeny_test_LIMS.svg",tree_style=ete_tree.tss, w=183, units="mm")



    def find_clusters(self, G, node_list):
        print("search clusters in", type(G))
        import itertools
        
        comb = itertools.combinations(node_list, 2)
        groups = []
        for i in comb:
            if i[0] == i[1]:
                continue
            try:
                print (G[i[0]][i[1]])
            except:
                no_match = True
                for n, group in enumerate(groups):
                    if i[0] in groups[n] and i[1] in groups[n]:
                        no_match = False
                    elif i[0] in groups[n] and i[1] not in groups[n]:
                        groups[n].append(i[1])
                        no_match = False
                    elif i[1] in groups[n] and i[0] not in groups[n]:
                        groups[n].append(i[0])
                        no_match = False
                if no_match:
                    groups.append([i[0], i[1]])
        return groups

    def merge_group_nodes(self, G, groups, node_list):

        node2merged_group_label = {}
        for group in groups:
            #print("one group", group)
            median_dico = {}
            for node in node_list:
                if node in group:
                    continue
                data = []
                for member in group:
                    if node in node2merged_group_label:
                        data.append(G[member][node2merged_group_label[node]]['weight'])
                    else:
                        data.append(G[member][node]['weight'])
                m = numpy.median(data)
                #print("adding median;", node)
                median_dico[node] = m
            mapping = {group[0]: '\n'.join([str(i) for i in group])}
            for node in group:
                node2merged_group_label[node] = '\n'.join([str(i) for i in group])
            G = networkx.relabel_nodes(G, mapping)
            for i in group[1:len(group)]:
                G.remove_node(i)
            #print("median_dico", median_dico)
            for i in node_list:
                #print("i, group", i, group, i in group)
                if i in group:
                    continue
                else:
                    med = median_dico[i]
                    if i in node2merged_group_label:
                        G['\n'.join([str(i) for i in group])][node2merged_group_label[i]]['weight'] = med
                    else:
                        G['\n'.join([str(i) for i in group])][i]['weight'] = med
        return G

    # this function is used to convert networkx to Cytoscape.js JSON format
    # returns string of JSON
    def convert2cytoscapeJSON(self, G):
        # load all nodes into nodes array
        final = {}
        final["nodes"] = []
        final["edges"] = []
        for node in G.nodes():
            nx = {}
            nx["data"] = {}
            nx["data"]["id"] = str(node)
            nx["data"]["label"] = str(node)
            final["nodes"].append(nx.copy())
        #load all edges to edges array
        for edge in G.edges(data=True):
            print(edge)
            nx = {}
            nx["data"]={}
            nx["data"]["id"]=str(edge[0])+str(edge[1])
            nx["data"]["strength"] = edge[2]["weight"]
            nx["data"]["source"]=str(edge[0])
            nx["data"]["target"]=str(edge[1])
            final["edges"].append(nx)
        return final

    def group_closely_related(self, G, group_list, cutoff):
        '''
        Input network 
        Extract closely related samples
        '''
        groups_cutoff_1 = []
        for group in group_list:
            group_updated = group.copy()
            group_label = '\n'.join([str(i) for i in group])
            for alt in G[group_label]:
                if G[group_label][alt]['weight'] <= cutoff:
                    lst = str(alt).split("\n")
                    group_updated = group_updated+lst
            groups_cutoff_1.append(group_updated)
        return groups_cutoff_1


    def get_MS_tree(self, analysis_id):

        #m = pandas.read_csv(dist_matrix, delimiter='\t', header=0, index_col=0)
        dist_matrix_df = self.get_snp_matrix(analysis_id)

        nodes = dist_matrix_df.index

        G = networkx.from_pandas_adjacency(dist_matrix_df)  # networkx.from_numpy_matrix(A)

        groups = self.find_clusters(G, nodes)

        G = self.merge_group_nodes(G, groups, nodes)

        T=networkx.minimum_spanning_tree(G)

        return G, T, groups


    def get_fastq_metadata_list_v2(self, 
                                   term_list=False, 
                                   fastq_filter=None, 
                                   run_name_list=False,
                                   metadata_value_list=False,
                                   add_molis=False,
                                   analysis_id_list=False,
                                   subproject_id_list=False,
                                   fastq_intersection=False):
        '''
        retrieve metadata from both sample and fastq metadata table
        TODO will curently combine dta from multiple anlayses => should return analysis_id as well
        '''

        res_filter_fastq = ''
        res_filter_sample = ''
        if term_list:
            print("term_list", term_list)
            term_filter = '","'.join(term_list)
            res_filter_fastq += f'and t2.name in ("{term_filter}")\n' 
            res_filter_sample += f'and t2.name in ("{term_filter}")\n' 
        if fastq_filter:
            fastq_filter_str = ','.join([str(i) for i in fastq_filter])
            res_filter_fastq += f'and fastq_id in ({fastq_filter_str})\n'
            res_filter_sample += f'and fastq_id in ({fastq_filter_str})\n' 
        if run_name_list:
            run_filter = '","'.join(run_name_list)
            res_filter_fastq += f'and run_name in ("{run_filter}")'
            res_filter_sample += f'and run_name in ("{run_filter}")'
        if metadata_value_list:
            metadata_filter = '","'.join(metadata_value_list)
            res_filter_fastq += f'and t1.value in ("{metadata_filter}")'
            res_filter_sample += f'and t1.value in ("{metadata_filter}")'
        if analysis_id_list:
            metadata_filter = '","'.join(analysis_id_list)
            res_filter_fastq += f'and t1.analysis_id in ("{metadata_filter}")'        
        if subproject_id_list:
            metadata_filter = '","'.join(subproject_id_list)
            res_filter_fastq += f'and t5.subproject_id in ("{metadata_filter}")'   
        
        sql_fastq = f'''
            select distinct t5.analysis_id,fastq_id, t2.name,t1.value, run_name from GEN_fastqfilesmetadata t1 
            inner join GEN_term t2 on t1.term_id=t2.id 
            inner join GEN_fastqfiles t3 on t1.fastq_id=t3.id 
            inner join GEN_runs t4 on t3.run_id=t4.id 
            left join GEN_projectanalysis t5 on t1.analysis_id=t5.analysis_id
            where t3.fastq_prefix not like "Undetermined%"
            {res_filter_fastq}
            group by fastq_id, t2.name,t3.fastq_prefix,t4.run_date,t4.run_name,t4.read_length
        '''
        
        sql_sample = f'''
            select distinct NULL as analysis_id, t3.fastq_id, t2.name,t1.value,run_name from GEN_samplemetadata t1 
            inner join GEN_term t2 on t1.term_id=t2.id 
            inner join GEN_fastqtosample t3 on t1.sample_id=t3.sample_id
            inner join GEN_fastqfiles t4 on t3.fastq_id=t4.id 
            inner join GEN_runs t5 on t4.run_id=t5.id 
            where t4.fastq_prefix not like "Undetermined%"
            {res_filter_sample}
            group by t3.fastq_id, t2.name,t4.fastq_prefix,t5.run_date,t5.run_name,t5.read_length
            '''

        df_fastq = pandas.read_sql(sql_fastq, self.conn)
        df_samples = pandas.read_sql(sql_sample, self.conn).set_index("fastq_id")
        
        if fastq_intersection:
            intersection = set(df_fastq["fastq_id"].to_list()).intersection(set(df_samples.index.to_list()))
            df_samples = df_samples.loc[intersection,]

        df = pandas.concat([df_fastq, df_samples.reset_index()])
        
        if add_molis:
            df_molis = self.get_fastq_and_sample_data(df["fastq_id"].to_list()).set_index("fastq_id")
            df = df.set_index("fastq_id").join(df_molis, on="fastq_id", rsuffix='_other')
            df = df[["molis_id", "name", "value", "run_name"]]

        return df

    def get_xslx_id2fastq_id_list(self,):
        '''
        TODO: deal with samples with multiple fastq ids 
        '''
        sql = '''select t2.xlsx_sample_ID,t1.fastq_id from GEN_fastqtosample t1 
                inner join GEN_sample t2 on t1.sample_id=t2.id '''
        self.cursor.execute(sql,) 

        xslx_id2fastq_id = {}
        for row in self.cursor.fetchall():
            if row[0] not in xslx_id2fastq_id:
                xslx_id2fastq_id[row[0]] = [row[1]]
            else:
                xslx_id2fastq_id[row[0]].append(row[1])

        return xslx_id2fastq_id

    def get_xslx_id2fastq_id(self,):
        '''
        TODO: deal with samples with multiple fastq ids 
        '''
        sql = '''select t2.xlsx_sample_ID,t1.fastq_id from GEN_fastqtosample t1 
                inner join GEN_sample t2 on t1.sample_id=t2.id '''
        self.cursor.execute(sql,) 

        return {int(i[0]):i[1] for i in self.cursor.fetchall()}

    def get_xslx_id2sample_id(self,):
        
        sql = '''select xlsx_sample_ID,id from GEN_sample'''
        self.cursor.execute(sql,) 
        return {int(i[0]):i[1] for i in self.cursor.fetchall()}


    def count_qc_warning(self, key_str=False):
        
        df = self.get_fastq_metadata_list()[["fastq_id","value"]]

        term2n_warnings = df.query('value == "WARN"').groupby(["fastq_id"])['value'].count().to_dict()

        if key_str:
            return {str(i):term2n_warnings[i] for i in term2n_warnings}
        else:
            return term2n_warnings




    def get_fastq_metadata_stats(self, term_list, analysis_id_list=False, workflow_id=False):
        
        if workflow_id:
            workflow_id_list = [workflow_id]
        else:
            workflow_id_list = False
        df = self.get_fastq_metadata_list(term_list, analysis_id_list=analysis_id_list, workflow_id_list=workflow_id_list)[["fastq_id","name","value"]]
        
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
        sql = f'''select distinct t1.id as fastq_id,fastq_prefix,R1,R2,species_name,molis_id, sample_name,t4.run_name,t2.sample_id from GEN_fastqfiles t1 
                left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
                left join GEN_sample t3 on t2.sample_id=t3.id
                inner join GEN_runs t4 on t1.run_id=t4.id
                where t1.id in ("{fastq_list_filter}");
            '''
        print(sql,)
        return pandas.read_sql(sql, self.conn)


    def get_analysis_metadata(self, term_name, analysis_id_list=False):
        

        add_filter = ''
        if analysis_id_list:
            filter_str = ','.join([str(i) for i in analysis_id_list])
            add_filter = f' and t1.analysis_id in ({filter_str})'

        sql = f"""select t1.analysis_id,value from GEN_analysismetadata t1 
                  inner join GEN_term t2 on t1.term_id =t2.id 
                  where t2.name like '{term_name}'
                  {add_filter}
                  """
        self.cursor.execute(sql,) 
        return {int(i[0]):i[1] for i in self.cursor.fetchall()}

    def parse_CT_scoring_table(self, table_path):
        
        df = pandas.read_csv(table_path, 
                               index_col=None,
                               header=0,
                               sep=",")

        df = df.fillna(value=False)
        term_name2data = {}
        for n, row in df.iterrows():
            term_name2data[row["name"]] = {"LCL_failed": row["LCL_failed"], "UCL_failed":row["UCL_failed"], "LCL_warning":row["LCL"], "UCL_warning":row["UCL"]}

        return term_name2data

    def qc_score(self, 
                 fastq_id_list, 
                 config,
                 workflow_id=False,
                 analysis_id_list=False,
                 index_int=False):

        from GEN.models import Workflow

        if workflow_id:
            workflow_name = Workflow.objects.filter(id=workflow_id)[0].workflow_name
        else:
            # use base metrics if not workflow specified
            workflow_name = 'BASE'

        print("workflow_name", workflow_name)

        metric2scoring = self.parse_CT_scoring_table(config["WORKFLOW"][workflow_name]["SCORING_TABLE"])

        print("metric2scoring", metric2scoring)

        if workflow_id:
            workflow_id_list = [workflow_id]
        else:
            workflow_id_list = False
        df = self.get_fastq_metadata_list(term_list=list(metric2scoring.keys()), 
                                          fastq_filter=fastq_id_list,
                                          analysis_id_list=analysis_id_list,
                                          workflow_id_list=workflow_id_list)

        # fastq_id, t2.name,t1.value, run_name
        # default score to 0
        fastq_id2n_fail = {str(i):0 for i in df["fastq_id"].to_list()}
        fastq_id2n_warn = {str(i):0 for i in df["fastq_id"].to_list()}
        fastq_id2metric2score = {str(i):{} for i in df["fastq_id"].to_list()}
                
        for n, row in df.iterrows(): 
            fastq_id = str(row["fastq_id"])

            if row["name"] in metric2scoring:
                
                # default green
                LCL_failed = metric2scoring[row["name"]]["LCL_failed"]
                UCL_failed = metric2scoring[row["name"]]["UCL_failed"]
                LCL_warning = metric2scoring[row["name"]]["LCL_warning"]
                UCL_warning = metric2scoring[row["name"]]["UCL_warning"]
                
                if LCL_failed:
                    if float(float(row["value"])) < float(LCL_failed):
                        print(row["name"])
                        print("LCL_failed", LCL_failed)
                        fastq_id2n_fail[fastq_id] += 1
                        fastq_id2metric2score[fastq_id][row["name"]] = 'FAIL'
                if UCL_failed:
                    if float(float(row["value"])) > float(UCL_failed):
                        print(row["name"])
                        print("UCL_failed", UCL_failed)
                        fastq_id2n_fail[fastq_id] += 1
                        fastq_id2metric2score[fastq_id][row["name"]] = 'FAIL'
                if row["name"] not in fastq_id2metric2score[fastq_id]:
                    if LCL_warning:
                        if float(float(row["value"])) < float(LCL_warning):
                            #print("LCL_warning", LCL_warning)
                            fastq_id2n_warn[fastq_id] += 1
                            fastq_id2metric2score[fastq_id][row["name"]] = 'WARN'
                    if UCL_warning:
                        
                        if float(float(row["value"])) > float(UCL_warning):
                            #print("UCL_warning", UCL_warning)
                            fastq_id2n_warn[fastq_id] += 1
                            fastq_id2metric2score[fastq_id][row["name"]] = 'WARN' 

        print("fastq_id2n_warn", fastq_id2n_warn)
        print("fastq_id2n_fail", fastq_id2n_fail)

        if index_int:
            fastq_id2n_fail = {int(k):v for k,v in fastq_id2n_fail.items()}
            fastq_id2n_warn = {int(k):v for k,v in fastq_id2n_warn.items()}
            fastq_id2metric2score = {int(k):v for k,v in fastq_id2metric2score.items()}
        return fastq_id2n_fail, fastq_id2n_warn, fastq_id2metric2score

    def get_run_name2run_id(self,):
        sql = 'select run_name,id from GEN_runs'
        self.cursor.execute(sql,) 
        return {i[0]:i[1] for i in self.cursor.fetchall()}

    def match_fastq_to_sample(self, fastq_prefix):

        # mapping based on xlsx sample id
        # conflict possible with several projects with sample with numeric labels (1,2,3, 232,...)
        # use date filter to exclude mapping older than X on the various columns
        sql = f'select id from GEN_sample where xlsx_sample_ID="{fastq_prefix}"' 
        try:
            self.cursor.execute(sql,) 
            sample_id = self.cursor.fetchall()[0][0]
        except:
            sql = f'select id from GEN_sample where sample_name="{fastq_prefix}"' 
            try:
                self.cursor.execute(sql,) 
                sample_id = self.cursor.fetchall()[0][0]
            except:
                sql = f'select id from GEN_sample where alias="{fastq_prefix}"' 
                try:
                    self.cursor.execute(sql,) 
                    sample_id = self.cursor.fetchall()[0][0]
                except:
                    sample_id = None 
        return sample_id

    def get_sample_id(self, sample_xls_id):
        sql = f'select id from GEN_sample where xlsx_sample_ID={self.spl}'
        self.cursor.execute(sql,(sample_xls_id,)) 
        return self.cursor.fetchall()[0][0]

    def add_sample_to_fastq_relation(self, fastq_id, sample_id):
        
        if GEN_settings.DB_DRIVER == 'sqlite':
            # ignore if relation already known
            sql2 = f'insert or ignore into GEN_fastqtosample(fastq_id, sample_id) values({self.spl},{self.spl})'
        elif GEN_settings.DB_DRIVER == 'mysql':
            # ignore if relation already known
            sql2 = f'insert ignore into GEN_fastqtosample(fastq_id, sample_id) values({self.spl},{self.spl})'
        else:                                                        
            raise IOError(f"Unknown db driver: {GEN_settings.DB_DRIVER}")
           
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
        


        if GEN_settings.DB_DRIVER == 'sqlite':
            sql = f'''INSERT into GEN_runs (run_name, run_date, assay, read_length, paired, filearc_folder) values({self.spl},{self.spl},{self.spl},{self.spl},{self.spl},{self.spl}) 
            ON CONFLICT(GEN_runs.run_name) DO UPDATE SET run_date={self.spl}, assay={self.spl}, read_length={self.spl}, paired={self.spl}, filearc_folder={self.spl};
            '''
        elif GEN_settings.DB_DRIVER == 'mysql':
            sql = f'''INSERT into GEN_runs (run_name, run_date, assay, read_length, paired, filearc_folder) values({self.spl},{self.spl},{self.spl},{self.spl},{self.spl},{self.spl}) 
            ON DUPLICATE KEY UPDATE run_date={self.spl}, assay={self.spl}, read_length={self.spl}, paired={self.spl}, filearc_folder={self.spl};
            '''
        else:                                                        
            raise IOError(f"Unknown db driver: {GEN_settings.DB_DRIVER}")



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

    def compare_samples(self, 
                        fastq_id_1, 
                        fastq_id2,
                        alt_freq_cutoff=70):
        
        sql = f'select * from GEN_snps where fastq_id={fastq_id_1} and alt_percent>{alt_freq_cutoff}'
        df1 = pandas.read_sql(sql, self.conn)
        sql2 = f'select * from GEN_snps where fastq_id={fastq_id2} and alt_percent>{alt_freq_cutoff}'
        df2 = pandas.read_sql(sql2, self.conn)
        df1["VAR"] = df1["ref"] + df1["position"].astype(str) + df1["alt"]
        df2["VAR"] = df2["ref"] + df2["position"].astype(str) + df2["alt"]

        df1_vars = set(df1["VAR"].to_list())
        df2_vars = set(df2["VAR"].to_list())
        
        shared = df1_vars.intersection(df2_vars)
        unique_df1 = df1_vars.difference(df2_vars)
        unique_df2= df2_vars.difference(df1_vars)
        union = df1_vars.union(df2_vars)

        n_shared = len(shared)
        n_unique_df1 = len(unique_df1)
        n_unique_df2 = len(unique_df2)
        n_union = len(union)
        n_diffs = n_unique_df1 +  n_unique_df2
        return [fastq_id_1, fastq_id2, n_union, n_shared, n_diffs, n_unique_df1, n_unique_df2]

    def get_analysis_from_description(self, description):
        sql = f'''select analysis_id from GEN_analysismetadata t1 inner join GEN_term t2 on t1.term_id=t2.id 
                 where t2.name='description' and t1.value="{description}";'''
                 
        self.cursor.execute(sql,)
        res = self.cursor.fetchall()
        
        if len(res) > 1:
            raise(f"Multiple analyses matche the description: {description}")
        elif len(res) == 0:
            return None
        else:
            return res[0][0]

    def create_search_index(self,):

        # sample table
        # sqlite: fts5 mysql: FULLTEXT
        # metadata 
        pass

    def get_analysis_overview(self, analysis_id_list=False):
        print(analysis_id_list)
        if isinstance(analysis_id_list, list):
            analysis_filter = '","'.join([str(i) for i in analysis_id_list])
            filter_str = f'and t1.analysis_id in ("{analysis_filter}") '
        else:
            filter_str = ''
        print("filter_str", filter_str)
        sql_analyse_metadata = f'''select t1.analysis_id,t4.value as description, start_date,workflow_name,count(*) as n_samples,t2.status from GEN_fastqset t1 
                                    inner join GEN_analysis t2 on t1.analysis_id=t2.id 
                                    inner join GEN_workflow t3 on t2.workflow_id=t3.id 
                                    inner join GEN_analysismetadata t4 on t1.analysis_id=t4.analysis_id
                                    inner join GEN_term t5 on t4.term_id=t5.id
                                    where t5.name='description' {filter_str} group by t1.analysis_id,t4.value,start_date,workflow_name'''
        print(sql_analyse_metadata)
        analyses = pandas.read_sql(sql_analyse_metadata, self.conn)

        return analyses


    def searchdb(self, term, search_type="exact_match"):

        if search_type=="exact_match":
            # search sample table 
            # sample_name, xlsx_sample_id, molis_id
            sql_sample = f'''select t1.id,t2.fastq_id,sample_name, xlsx_sample_id, molis_id from GEN_sample t1 
                            inner join GEN_fastqtosample t2 on t1.id=t2.sample_id where t1.id={self.spl} or sample_name={self.spl} or xlsx_sample_id={self.spl} or molis_id={self.spl} limit 1000'''

            # search metadata fastq 
            sql_fastq_metadata = f'''select t2.sample_id,t1.fastq_id, value,t3.name from GEN_fastqfilesmetadata t1 
                                    inner join GEN_fastqtosample t2 on t1.fastq_id=t2.fastq_id 
                                    inner join GEN_term t3 on t1.term_id=t3.id
                                    where t1.fastq_id={self.spl} or value={self.spl} limit 1000
                                '''

            # search metadata sample
            
            sql_sample_metadata = f'''select t1.sample_id,t2.fastq_id, value,t3.name from GEN_samplemetadata t1 
                                     inner join GEN_fastqtosample t2 on t1.sample_id=t2.sample_id 
                                     inner join GEN_term t3 on t1.term_id=t3.id
                                     where value={self.spl} limit 1000
                                   '''

            hits_sample_df = pandas.read_sql(sql_sample, self.conn, params=[term,term,term,term])
            hits_fastq_metadata_df = pandas.read_sql(sql_fastq_metadata, self.conn, params=[term,term])
            hits_sample_metadata_df = pandas.read_sql(sql_sample_metadata, self.conn, params=[term])

            nr_fastq_list = [] 
            if not hits_sample_df.empty:
                nr_fastq_list += hits_sample_df["fastq_id"].to_list()
            if not hits_sample_metadata_df.empty:
                nr_fastq_list += hits_sample_metadata_df["fastq_id"].to_list()
            if not hits_fastq_metadata_df.empty:
                nr_fastq_list += hits_fastq_metadata_df["fastq_id"].to_list()

            nr_fastq_list =  [str(i) for i in set(nr_fastq_list)]

            detail_df = self.get_fastq_and_sample_data(nr_fastq_list)[["fastq_id","fastq_prefix","run_name","species_name","molis_id", "sample_name"]]

            if not hits_sample_metadata_df.empty:
                # add metadata to detail_df
                detail_df = detail_df.set_index("fastq_id").join(hits_sample_metadata_df.set_index("fastq_id"))
            if not hits_fastq_metadata_df.empty:
                if 'fastq_id' in detail_df.columns:
                    detail_df = detail_df.set_index("fastq_id").join(hits_fastq_metadata_df[["fastq_id", "value", "name"]].set_index("fastq_id"))
                else:
                    detail_df = detail_df.join(hits_fastq_metadata_df[["fastq_id", "value", "name"]].set_index("fastq_id"), rsuffix='r')
            

            nr_fastq_list_filter = ','.join(nr_fastq_list)

            sql_analyse = f'''select distinct fastq_id,analysis_id from GEN_fastqset t1 
                              where fastq_id in ({nr_fastq_list_filter})'''
            
            analyses = pandas.read_sql(sql_analyse, self.conn)


            # add qc results                       
            df_qc = self.get_fastq_metadata_v2("qc_status", fastq_id_list=nr_fastq_list)
            analyses = pandas.merge(analyses, df_qc,  how='left', left_on=['fastq_id','analysis_id'], right_on = ['fastq_id','analysis_id'])

            # add n samples, description
            analyse_metadata = self.get_analysis_overview(analyses["analysis_id"].to_list())
            analyses = pandas.merge(analyses, analyse_metadata,  how='left', left_on=['analysis_id'], right_on = ['analysis_id'])

            return detail_df, analyses
    
    def compare_snp_tables(self, fasq_id_list, min_alt_freq):

        df_list = [self.fastq_snp_table(fastq_id, min_alt_freq) for fastq_id in fasq_id_list]
        
        df_1 = df_list[0].set_index("nucl_change")[['depth', 'alt_percent']]
        
        print("df_1", df_1.head())

        df_1.columns = [f"depth_{fasq_id_list[0]}", f"alt_percent_{fasq_id_list[0]}"]
        
        for n,df in enumerate(df_list[1:]):
            #print(n)
            df_2 = df.set_index("nucl_change")[['depth', 'alt_percent']]
            df_2.columns = [f"depth_{fasq_id_list[n+1]}", f"alt_percent_{fasq_id_list[n+1]}"]
            df_1 = df_1.join(df_2, how='outer')
        
        print("combined", df_1.head())

        df_1["diff"] = 0
        #print(df_1)
        #print(df_1.isna().any(axis=1))
        df_1.loc[df_1.isna().any(axis=1), "diff"] = 1
        #print(df_1[df_1.isna().any(axis=1)])
    
        snp2frequency = self.snp_frequency(min_alt_freq=min_alt_freq)

        summary_data, matrix = self.parwise_snps_comp(fasq_id_list, min_alt_freq)

        df_1["frequency"] = [snp2frequency[i] for i in df_1.index]

        return df_1, summary_data

    def parwise_snps_comp(self, fastq_list, alt_freq_cutoff=70):
        import itertools
        comb = itertools.combinations(fastq_list, 2)

        res = []
        complete_mat = []
        for one_pair in comb:
            print("one_pair", one_pair)
            vals = self.compare_samples(one_pair[0], one_pair[1], alt_freq_cutoff=alt_freq_cutoff)
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
                      sample_xls_id,
                      update=False):
        from MySQLdb._exceptions import ProgrammingError
        

        # INSERT into GEN_sample(xlsx_sample_ID,species_name,date_received,sample_name,sample_type,analysis_type,description,molis_id,myseq_passage,run_date,date_registered,date_sample_modification,user_creation_id,user_modification_id) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(GEN_sample.xlsx_sample_ID) DO UPDATE SET xlsx_sample_ID=?,species_name=?,date_received=?,sample_name=?,sample_type=?,analysis_type=?,description=?,molis_id=?,myseq_passage=?,run_date=?,date_registered=?,date_sample_modification=?,user_creation_id=?,user_modification_id=?;
        # [29, 'Staphylococcus aureus', False, nan, 'strain', 'research', nan, 1306182530, nan, False, '2020-09-02', '2020-09-02', 2, 2, 29, 'Staphylococcus aureus', False, nan, 'strain', 'research', nan, 1306182530, nan, False, '2020-09-02', '2020-09-02', 2, 2]            
        # update all columns in case of conflict with xlsx_sample_id
        

        # NOTE: xlsx_sample_ID used as reference: if a row is updated in the xlsx table, the corresponding row is updated in the sql table
        print(GEN_settings.DB_DRIVER)
        if GEN_settings.DB_DRIVER == 'sqlite':
            update_str = f'%s={self.spl}'
            update_str_comb = ','.join([update_str % colname for colname in col_names])

            sql_template = f'''

            INSERT into GEN_sample(%s) values(%s)
                            ON CONFLICT(GEN_sample.xlsx_sample_ID) DO UPDATE SET %s;''' % (','.join(col_names),
                                                                                           ','.join([f'{self.spl}']*len(col_names)),
                                                                                           update_str_comb)
            print(sql_template, values_list + values_list)
            self.cursor.execute(sql_template, values_list + values_list)
            
        elif GEN_settings.DB_DRIVER == 'mysql':
            update_str = f'%s=%{self.spl}'
            update_str_comb = ','.join([update_str % colname for colname in col_names])
            if not update:
                sql_template = '''INSERT into GEN_sample(%s) values(%s);''' % (','.join(col_names),
                                                                               ','.join([f'{self.spl}']*len(col_names)))
            else:
                sql_template = '''UPDATE GEN_sample SET
                                 %s where xlsx_sample_ID=%s;''' % (update_str_comb,
                                                                   sample_xls_id)
            print(sql_template, values_list)
            #print(values_list)
            try:
                self.cursor.execute(sql_template, values_list)
            except IntegrityError as e:
                if 'UNIQUE' in str(e) or 'Duplicate' in str(e):
                    print(f'UNIQUE constraint failed for sample ID: {sample_xls_id} -- skipping row')
                    return None
                else:
                    print(f"Problem with sample ID: {sample_xls_id}")
                    raise                                               

        else:                                                        
            raise IOError(f"Unknown db driver: {GEN_settings.DB_DRIVER}")
        



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
        sql = f'select id from GEN_fastqfiles where fastq_prefix={self.spl}'
        
        try:
            self.cursor.execute(sql,(sample_prefix,)) 
            fastq_id_list = [i[0] for i in self.cursor.fetchall()]
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
        self.cursor.execute(sql,) 
        data = [list(i) for i in self.cursor.fetchall()]
        for n, row in enumerate(data):
            if row[4]:
                data[n][4] = Analysis.objects.filter(id=row[4])[0]
        return data

    def get_run_samples(self, run_name):
        
        sql = f'''select fastq_prefix,R1,R2 from GEN_fastqfiles t1 
                  inner join GEN_runs t2 on t1.run_id=t2.id where run_name="{run_name}"
               '''
        self.cursor.execute(sql,) 
        return self.cursor.fetchall()

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
        self.cursor.execute(sql,) 
        return {i[0]: i[1] for i in self.cursor.fetchall()}
    
        
    def fastq_prefix2fastq_id(self,fastq_prefix_list):
        
        fasta_filter = '","'.join(fastq_prefix_list)
        sql = f'select fastq_prefix,id from GEN_fastqfiles where fastq_prefix in ("{fasta_filter}")'
        self.cursor.execute(sql,)
        return {i[0]:i[1] for i in self.cursor.fetchall()}
    
    def add_fastq_metadata(self, 
                           fastq_id,
                           term_name,
                           value,
                           analysis_id=False):
        '''
        skip duplicates by default
        [{"fastq_id": fastq_id,
        "term_name": <name>,
        "value": <value>,
        "analysis_id" analysis_id}]
        '''
        from GEN.models import Term
        from GEN.models import FastqFilesMetadata
        term = Term.objects.get_or_create(name=term_name)[0]
        try:
            if not analysis_id:
                m = FastqFilesMetadata(term=term, fastq_id=fastq_id, value=value)
                m.save()
            else:
                m = FastqFilesMetadata(term=term, fastq_id=fastq_id, value=value, analysis_id=analysis_id)
                m.save()
        except IntegrityError as e:
            if 'UNIQUE' in str(e) or 'Duplicate' in str(e):
                  #print(f'UNIQUE constraint failed: {sample_id} -- {key}')
                  pass
            else:
                  print("other error")
                  raise

    def add_sample_metadata(self, 
                           sample_id,
                           term_name,
                           value,
                           update=False):
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


        if not update:
            m = SampleMetadata(term=term, sample_id=sample_id, value=value)
            m.save()
        else:
            # update value
            try:
                m = SampleMetadata.objects.filter(term=term, sample_id=sample_id)[0]
                # update entry only if necessary
                if m.value == value:
                    return
                m.value = value
            except:
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

    def fastq_mutation(self, 
                       aa_change_list, 
                       analysis_id_list=False, 
                       fastq_id_list=False,
                       min_alt_freq=70,
                       return_df=False):
        
        change_filter = '","'.join(aa_change_list)
        if analysis_id_list:
            analysis_list = ','.join([str(i) for i in analysis_id_list])
            analysis_filter = f'and analysis_id in ({analysis_list})'
        else:
            analysis_filter = ''
        if fastq_id_list:
            fastq_list = ','.join([str(i) for i in fastq_id_list])
            fastq_filter = f'and fastq_id in ({fastq_list})'
        else:
            fastq_filter = ''
        sql = f'''select fastq_id,aa_change from GEN_snps where aa_change in ("{change_filter}") 
                  {analysis_filter} 
                  {fastq_filter}
                  and alt_percent>{min_alt_freq};'''
        df = pandas.read_sql(sql, self.conn)
        if return_df:
            return df
        else:
            fastq2changes = {i:[] for i in df["fastq_id"].to_list()}
            for n, row in df.iterrows():
                fastq2changes[row["fastq_id"]].append(row["aa_change"])
            return fastq2changes


    def snp_to_fastq(self, nucl_change, min_alt_freq=70, type='nucl_change'):
        
        sql = f'select distinct analysis_id,fastq_id,position,ref,alt,gene,aa_change,nucl_change,depth,ref_count,alt_count,alt_percent from GEN_snps where {type}="{nucl_change}" and alt_percent>{min_alt_freq};'

        df = pandas.read_sql(sql, self.conn)

        sql2 = f'select distinct fastq_id,position,ref,alt,gene,aa_change,nucl_change,depth,ref_count,alt_count,alt_percent from GEN_snps where alt_percent>{min_alt_freq};'

        df2 = pandas.read_sql(sql2, self.conn)

        # add number of snp present in each sample
        fastq_id2n_snps = df2.groupby(["fastq_id"])["fastq_id"].count().to_dict()
        df["n_snps"] = [fastq_id2n_snps[i] for i in df["fastq_id"]]

        # add run name
        fastq2run_name = self.get_fastq_id2run_name()
        df["run_name"] = [fastq2run_name[str(i)] for i in df["fastq_id"]]

        # add patient id
        fastq_id2patient_id = self.get_fastq_metadata("patient_id", index_str=False)
        df["patient_id"] = [fastq_id2patient_id[i] if i in fastq_id2patient_id else '-' for i in df["fastq_id"]]

        # add patient id
        fastq_id2pangolin_lineage = self.get_fastq_metadata("pangolin_lineage", index_str=False)
        df["pangolin_lineage"] = [fastq_id2pangolin_lineage[i] if i in fastq_id2pangolin_lineage else '-' for i in df["fastq_id"]]

        return df


    def fastq_snp_table(self, 
                        fastq_id, 
                        min_alt_freq=70,
                        analysis_id_list=False):
        
        add_filter = ''
        if analysis_id_list:
            analysis_filter = ','.join([str(i) for i in analysis_id_list])
            add_filter += f' and analysis_id in ({analysis_filter})'

        sql = f'''select distinct fastq_id,reference,position,ref,alt,effect_id,gene,aa_position,aa_change,nucl_change,depth,ref_count,alt_count,alt_percent from GEN_snps 
                  where fastq_id={fastq_id} 
                  and alt_percent>{min_alt_freq} {add_filter};'''

        return pandas.read_sql(sql, self.conn)

    def snp_frequency(self, 
                      analysis_id=False, 
                      percent=False,
                      min_alt_freq=70):
        if analysis_id:
            sql = f'select distinct fastq_id,nucl_change from GEN_snps where analysis_id={analysis_id} and alt_percent>{min_alt_freq} ;'
        else:
            sql = f'select distinct fastq_id,nucl_change from GEN_snps where alt_percent > {min_alt_freq}'

        df_snps = pandas.read_sql(sql, self.conn)

        if not percent:
            return df_snps.groupby("nucl_change")["nucl_change"].count().to_dict()
        else:
            mut_freq = pd.DataFrame(df_snps.groupby("nucl_change")["nucl_change"].count())
            n_genomes = len(df_snps["fastq_id"].unique())
            mut_freq.columns = ["mutation_frequency_count"]
            mut_freq["mutation_frequency_percent"] = [(i/float(n_genomes))*100 for i in mut_freq["mutation_frequency_count"]]
            return mut_freq.set_index("nucl_change")["mutation_frequency_percent"].to_dict()


    def add_QC_report(self, run_name, run_path):
        
        sql = f'update GEN_runs set qc=1 where run_name={self.spl}'
        self.cursor.execute(sql, [run_name]) 
        
        sql = f'update GEN_runs set qc_path={self.spl} where run_name={self.spl}'
        self.cursor.execute(sql, (run_path, run_name))
        self.conn.commit()
    
    
    def fastq_qc_filter(self, analysis_id, value, fastq_id_list=False):
        
        if fastq_id_list:
            fastq_id_filter = ','.join([str(i) for i in fastq_id_list])
            filter_str = f' and fastq_id in ({fastq_id_filter})'
        else:
            filter_str = ''
                    
        sql = f'''select fastq_id from GEN_fastqfilesmetadata t1 
                 inner join GEN_term t2 on t1.term_id=t2.id 
                 where t1.analysis_id={analysis_id} and t2.name="qc_status" and t1.value="{value}" {filter_str};
        '''
        self.cursor.execute(sql,)
        
        return [str(i[0]) for i in self.cursor.fetchall()]
        
    def format_snps(self, fastq_id_list, alt_freq_cutoff=70):
        
        fastq_filter = ','.join([str(i) for i in fastq_id_list])
        sql = f'select fastq_id,nucl_change, aa_change, gene from GEN_snps where fastq_id in ({fastq_filter}) and alt_percent>{alt_freq_cutoff}'

        df = pandas.read_sql(sql, self.conn).set_index("fastq_id")

        fastq2data = {str(fastq_id): {} for fastq_id in df.index.unique()}
        
        for fastq_id in list(fastq2data.keys()):
            target_fastq = df.loc[fastq_id]
            aa_changes = [f'{pos["aa_change"]} ({pos["gene"]})' for n, pos in target_fastq.iterrows() if not pandas.isna(pos["aa_change"])]
            nucl_changes = [str(pos["nucl_change"]) for n, pos in target_fastq.iterrows()]
            fastq2data[str(fastq_id)]["aa_changes"] = ';'.join(aa_changes)
            fastq2data[str(fastq_id)]["nucl_changes"] = ';'.join(nucl_changes)
            
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
        self.cursor.execute(sql,)
        return {i[0]:i[1] for i in self.cursor.fetchall()}
    
    def get_fastq_id2species(self, fastq_list):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_list])
        sql = f'''select distinct fastq_prefix,species_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where t1.id in ("{fastq_list_filter}");
           '''
        self.cursor.execute(sql,)
        return {i[0]:i[1] for i in self.cursor.fetchall()}
    
    def get_fastq_id2sample_name(self, 
                                 fastq_list, 
                                 key_str=True):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_list])
        sql = f'''select distinct fastq_id,sample_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where t1.id in ("{fastq_list_filter}");
           '''
        self.cursor.execute(sql,)
        if key_str:
            return {str(i[0]):i[1] for i in self.cursor.fetchall()}
        else:
            return {int(i[0]):i[1] for i in self.cursor.fetchall()}


    def generate_sample_file_GENCOV(self,fastq_id_list):

        '''
        Yaml of the form 
        '1016073361_1546':
            alt_id: Sample_1
            read1: /scratch/hdd/tpillone/projects/airflow_test/epidemiology/2021_01_05-1501_Z1UBOZ/cov_minipipe/1016073361_1546_R1.fastq.gz
            read2: /scratch/hdd/tpillone/projects/airflow_test/epidemiology/2021_01_05-1501_Z1UBOZ/cov_minipipe/1016073361_1546_R2.fastq.gz
        '''
        # id,fastq_prefix,R1,R2,species_name
        fastq_df = self.get_fastq_and_sample_data(fastq_id_list)
        
        dict_file = {}

        for n, row in fastq_df.iterrows():
            
            # deal with multiple fastq with same name
            R1 = row["R1"]
            R2 = row["R2"]
            #if '20210201' in R1:
            #    R1 = re.sub("20210201", "20210201_CleanPlex", row["R1"])
            #    R2 = re.sub("20210201", "20210201_CleanPlex", row["R2"])
            #if '20210208' in R1:
            #    R1 = re.sub("20210208", "20210208_CleanPlex", row["R1"])
            #    R2 = re.sub("20210208", "20210208_CleanPlex", row["R2"])
            fastq_id = row["fastq_id"]
            sample_name = f'{row["fastq_prefix"]}_{fastq_id}'
            dict_file[sample_name] = {'alt_id': f'Sample_{n+1}',
                                    'read1': R1,
                                    'read2': R2,
                                    }
            
        return dict_file


    def get_sample_table(self, suproject_id_list=False):
        

        project_filter = ''
        if suproject_id_list:

            project_filter = 'where t3.subproject_id in (%s)' % ','.join([str(i) for i in suproject_id_list])

        sql = f'''
        select distinct A.id,A.sample_type,A.sample_name,A.taxonomy,A.date_registered, A.date_received,A.n_fastq, count(t3.subproject_id ) as n_projects from (
        select distinct t1.id,t1.sample_type,t1.sample_name,t1.taxonomy,t1.date_registered, t1.date_received, count(t2.fastq_id) as n_fastq from GEN_sample t1 
        left join GEN_fastqtosample t2 on t1.id=t2.sample_id 
        group by t1.id) A  
        left join GEN_subprojectsample t3 on A.id=t3.sample_id
        {project_filter}
        group by A.id 
        '''
        self.cursor.execute(sql,) 
        return self.cursor.fetchall()
    
    def get_analysis_fastq_list(self,analisis_id):
        
        sql = f'''
        select fastq_id from GEN_fastqset where analysis_id={analisis_id}
        '''
        self.cursor.execute(sql,) 
        return [i[0] for i in self.cursor.fetchall()]
        

    def calculate_age(self, birth_date, prel_date, fraction_year=False):
        from dateutil.relativedelta import relativedelta
        rdelta = relativedelta(prel_date, birth_date)
        age_years = rdelta.years 
        if age_years == 0:
            age_months = rdelta.months
            if age_months == 0:
                age_weeks = rdelta.weeks
                if fraction_year:
                    return round(age_weeks/52.1429, 2)
                else:
                    return f'{age_weeks} weeks'
            else:
                if fraction_year:
                    return round(age_months/12, 2)
                else:
                    return f'{age_months} months'
        else:
            return age_years


    def insert_molis_sample_metadata_from_xml(self, 
                                              df, 
                                              field_definition, 
                                              molis_alias="Numéro alias"):
        import datetime

        print("field_definition", field_definition)
        
        metadata_list = []

        molis2sample_list = self.get_molis_id2sample_id_list()

        print("df shape", df.shape)
        for n, row in df.iterrows():
            alias = row[molis_alias]

            print("alias", alias)
            try:
                sample_list = molis2sample_list[str(alias)]
            except KeyError:
                print(f"No samples for alias: {alias}")
                continue
            print("number of samples", len(sample_list))
            for sample_id in sample_list:          
                print(f"Sample {sample_id} (alias {alias})")

                for field in field_definition:
                    # deal with incomplete tables
                    if field_definition[field]["fields"][0] not in df.columns:
                        print(f'field {field_definition[field]["fields"][0]} not in table')
                        continue

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
                        if 'zfill' in field_definition[field]:
                            print('ZFILL')
                            print('zipping1', values)
                            print('zipping2', field_definition[field]['zfill'])
                            print(zip(values, field_definition[field]['zfill']))
                            val = ''.join([str(x).zfill(y) for x,y in zip(values, field_definition[field]['zfill'])])
                        else:
                            val = ''.join(values)
                    elif field_definition[field]["type"] == 'age':
                        date_a = datetime.datetime.strptime(row[field_definition[field]["fields"][0]], '%Y-%m-%d')
                        date_b = datetime.datetime.strptime(row[field_definition[field]["fields"][1]], '%Y-%m-%d')
                        val = self.calculate_age(date_a, date_b).strftime('%Y-%m-%d')
                    else:
                        field_type = field_definition[field]["type"]
                        print(f'Unknown type: {field_type}')
                        raise IOError(f'Unknown type: {field_type}')
                    print({'term_name':field, 'value': val, 'sample_id': sample_id})
                    # add value if not empty
                    if val != '':
                            metadata_list.append({'term_name':field, 'value': val, 'sample_id': sample_id})

        # insert data
        for metadata in metadata_list:
            # duplicates will raise error...
            # print(metadata)
            #try:
            self.add_sample_metadata(sample_id=metadata["sample_id"],
                                     term_name=metadata["term_name"],
                                     value=metadata["value"],
                                     update=True) 
            '''
            except IntegrityError as e:
                if 'UNIQUE' in str(e) or 'Duplicate' in str(e):
                    print(f'UNIQUE constraint failed: {sample_id} -- {metadata["term_name"]} -- {metadata["value"]}')
                    continue
                else:
                    print("other error")
                    raise
            '''



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
        try:
            m = AnalysisMetadata.objects.filter(term=term, analysis_id=analysis_id)[0]
            m.value = value
        except:
            m = AnalysisMetadata(term=term, analysis_id=analysis_id, value=value)
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

