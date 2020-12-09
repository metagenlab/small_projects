

import os
import yaml
from datetime import datetime
from metagenlab_libs import gendb_utils
from airflow.hooks.base_hook import BaseHook
from airflow.contrib.operators.slack_webhook_operator import SlackWebhookOperator
import string 
import random
    
SLACK_CONN_ID = 'slack'

def clean_species(species_string):
    import re 
    # case when fastq file could not be matched to sample table
    if species_string is None:
        return "Unspecified species"
    if len(species_string) == 0:
        return "Unspecified species"
    species = re.sub(",.*", "",species_string)
    species = " ".join(species.split(" ")[0:2])
    return species 


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def make_run_dir(execution_folder, 
                 analysis_id):
    import shutil
    import errno
    
    folder_name = datetime.now().strftime("%Y_%m_%d-%H%M_" + id_generator(6))
    
    run_execution_folder = os.path.join(execution_folder, folder_name)
    
    if analysis_id:
        # update status
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "running", update=True)
        # save execution folder
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_folder", run_execution_folder)
    
    # remove existing dir
    shutil.rmtree(run_execution_folder, ignore_errors=True)
    # make dir
    os.makedirs(run_execution_folder)
    
    return folder_name

def write_sample_file(gen_db,
                      fastq_list,
                      analysis_name, 
                      execution_folder,
                      analysis_id):
    
    # update status
    gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "running", update=True)
    
    if not isinstance(fastq_list, list):
        fastq_list = fastq_list.split(",")
    
    # id,fastq_prefix,R1,R2,species_name
    fastq_df = gen_db.get_fastq_and_sample_data(fastq_list)
    
    run_execution_folder = os.path.join(execution_folder, analysis_name)
    
    header = ["SampleName",
              "ScientificName",
              "R1",
              "R2", 
              "fastq_id"]
    
    with open(os.path.join(run_execution_folder, f'{analysis_name}.tsv'), 'w') as f:
        f.write("\t".join(header) + '\n')
        for n, row in fastq_df.iterrows():
            
            # deal with multiple fastq with same name
            R1 = row["R1"]
            R2 = row["R2"]
            fastq_id = row["fastq_id"]
            species = row["species_name"]
            sample_name = f'{row["fastq_prefix"]}_{fastq_id}'
            f.write(f"{sample_name}\t{species}\t{R1}\t{R2}\t{fastq_id}\n")
            
            
            
def write_snakemake_config_file(analysis_name,
                                fastq_list,
                                execution_folder,
                                snakemake_config,
                                gen_db,
                                analysis_id,
                                reference_list=False,
                                scientific_name=False,
                                check_single_species=False,
                                reference_docx=False):

    # update status
    gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "running", update=True)

    run_execution_folder = os.path.join(execution_folder, analysis_name)
    
    if check_single_species and not scientific_name:
        species_list = list(set(gen_db.get_fastq_id2species(fastq_list.split(",")).values()))
    if check_single_species:
        if len(species_list) > 1:
            raise IOError("More than one different species in the dataset: %s" % ','.join(species_list))
        else:
            scientific_name = species_list[0]
    
    # if references, prepare list
    if reference_list:
        reference_fastq_list = reference_list.split(",")
        fastq_df = gen_db.get_fastq_and_sample_data(reference_fastq_list)
        
        ref_list = []
        for n, row in fastq_df.iterrows():
            fastq_id = row["fastq_id"]
            sample_name = f'{row["fastq_prefix"]}_{fastq_id}'
            ref_list.append(sample_name)
    
    with open(os.path.join(run_execution_folder, f'{analysis_name}.config'), 'w') as f:
        # update sample table name
        snakemake_config["local_samples"] = f'{analysis_name}.tsv'
        if reference_list:
            snakemake_config["reference"] = f'{",".join(ref_list)}'
        if check_single_species:
            snakemake_config["species"] = f'{scientific_name}'
        if reference_docx:
            snakemake_config["reference_docx"] = f'{reference_docx}'

        documents = yaml.dump(snakemake_config, f)
        

def backup(execution_folder, 
           backup_folder,
           file_or_folder_list,
           analysis_id=False,
           analysis_metadata={}):
    
    '''
    Analysis name: folder within execution_folder (generally execution date)
    Analysis_metadata: dictionnary of status to add to LIMS: 
        {"analysis_id" : <id>,
         "value": <value>}
    '''
    
    import shutil
    
    if analysis_id:
        # update status
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "running", update=True)

    
    for output in file_or_folder_list:
        original = os.path.join(execution_folder, output)
        target = os.path.join(backup_folder, output)
        print("original", original)
        print("target", target)
        if os.path.isdir(original):
            shutil.rmtree(target, ignore_errors=True)
            shutil.copytree(original, target)
        if os.path.isfile(original): 
            if os.path.exists(target):
                os.remove(target)
            shutil.rmtree(target, ignore_errors=True)
            shutil.copy(original, target)
    
    if analysis_metadata:
        for status_entry in analysis_metadata:
            gendb_utils.add_analysis_metadata(analysis_metadata[status_entry]["analysis_id"], 
                                              status_entry, 
                                              analysis_metadata[status_entry]["value"], 
                                              update=False)


def task_fail_slack_alert(context):

    print("task_fail_slack_alert")
    
    slack_webhook_token = BaseHook.get_connection(SLACK_CONN_ID).password

    dag_run = context.get('dag_run')
    analysis_id = dag_run.conf.get('analysis_id')
    
    if analysis_id:
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "failed", update=True)
        gendb_utils.update_analysis_status(analysis_id, "failed")
    else:
        print("NO ANALYSE NUMBER")

    slack_msg = """
            :red_circle: Task Failed. 
            *Task*: {task}  
            *Dag*: {dag} 
            *Execution Time*: {exec_date}  
            *Log Url*: {log_url} 
            *Analysis*: {analysis_id}
            """.format(
            task=context.get('task_instance').task_id,
            dag=context.get('task_instance').dag_id,
            exec_date=context.get('execution_date'),
            log_url=context.get('task_instance').log_url,
            analysis_id=analysis_id,
        )
          
    print("message:\n", slack_msg)
            
    failed_alert = SlackWebhookOperator(
        task_id='slack_test',
        http_conn_id='slack',
        webhook_token=slack_webhook_token,
        message=slack_msg,
        username='airflow')

    return failed_alert.execute(context=context)


def task_success_slack_alert(context):
    
    '''
    Slack message + update analysis status
    '''
    
    slack_webhook_token = BaseHook.get_connection(SLACK_CONN_ID).password

    dag_run = context.get('dag_run')
    analysis_id = dag_run.conf.get('analysis_id')

    if analysis_id:
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "success", update=True)
        gendb_utils.update_analysis_status(analysis_id, "success")

    slack_msg = """
            :heavy_check_mark: Dag Success. 
            *Dag*: {dag} 
            *Execution Time*: {exec_date}  
            *Log Url*: {log_url}
            *Analysis ID*: {analysis_id}
            """.format(
            dag=context.get('task_instance').dag_id,
            ti=context.get('task_instance'),
            exec_date=context.get('execution_date'),
            log_url=context.get('task_instance').log_url,
            analysis_id=analysis_id
        )
    
    success_alert = SlackWebhookOperator(
        task_id='slack_test',
        http_conn_id='slack',
        webhook_token=slack_webhook_token,
        message=slack_msg,
        username='airflow')

    return success_alert.execute(context=context)